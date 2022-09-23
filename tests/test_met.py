# coding: utf-8

_all_ = [ 'test_met' ]

import os
import sys
parent_dir = os.path.abspath(__file__ + 2 * '/..')
sys.path.insert(0, parent_dir)
import argparse
import glob
import multiprocessing
import itertools as it

import inclusion
from inclusion.config import main
from inclusion.utils import utils

import ROOT

def check_bit(bitpos, bit):
    bitdigit = 1
    res = bool(bit&(bitdigit<<bitpos))
    return res

def get_trigger_bit(trigger_name):
    """
    Returns the trigger bit corresponding to 'main.trig_map'
    """
    res = main.trig_map[trigger_name]
    try:
        res = res['mc']
    except KeyError:
        print('You likely forgot to add your custom trigger to `main.trig_custom`.')
        raise
    return res

def pass_triggers(trigs, bit):
    """
    Internal only function.
    Checks at least one trigger was fired.
    """
    flag = False
    for trig in trigs:
        if trig in main.trig_custom:
            flag = set_custom_trigger_bit(trig, bit)
        else:
            flag = check_bit(get_trigger_bit(trig), bit)
        if flag:
            return True
    return False    

def set_custom_trigger_bit(trigger, bit):
    """
    The VBF trigger was updated during data taking, adding HPS
    https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauTrigger
    """
    if trigger not in main.trig_custom:
        import inspect
        currentFunction = inspect.getframeinfo(frame).function
        raise ValueError('[{}] option not supported.'.format(currentFunction))

    if trigger == 'VBFTauCustom':
        bits = check_bit(main.trig_map[trigger]['VBFTauHPS']['mc'], bit)
    elif trigger == 'IsoDoubleTauCustom':
        bits = check_bit(main.trig_map[trigger]['IsoDoubleTauHPS']['mc'], bit)
    elif trigger == 'IsoMuIsoTauCustom':
        bits = check_bit(main.trig_map[trigger]['IsoMuIsoTauHPS']['mc'], bit)
    elif trigger == 'EleIsoTauCustom':
        bits = check_bit(main.trig_map[trigger]['EleIsoTauHPS']['mc'], bit)

    return bits

def sel_category(entries, category):
    btagLL = entries.bjet1_bID_deepFlavor > 0.0490 and entries.bjet2_bID_deepFlavor > 0.0490
    btagM  = ((entries.bjet1_bID_deepFlavor > 0.2783 and entries.bjet2_bID_deepFlavor < 0.2783) or
              (entries.bjet1_bID_deepFlavor < 0.2783 and entries.bjet2_bID_deepFlavor > 0.2783))
    btagMM = entries.bjet1_bID_deepFlavor > 0.2783 and entries.bjet2_bID_deepFlavor > 0.2783

    common = not (entries.isVBF == 1 and entries.VBFjj_mass > 500 and entries.VBFjj_deltaEta > 3 and
                  (entries.bjet1_bID_deepFlavor > 0.2783 or entries.bjet2_bID_deepFlavor > 0.2783))

    if category == 's1b1jresolvedMcut':
        specific = entries.isBoosted != 1 and btagM
    elif category == 's2b0jresolvedMcut':
        specific = entries.isBoosted != 1 and btagMM
    elif category == 'sboostedLLMcut':
        specific = entries.isBoosted == 1 and btagLL

    return common and specific

def sel_cuts(entries, iso_cuts=dict(), lepton_veto=True,
             bjets_cut=True, invert_mass_cut=True):
        """
        Applies selection cut to one event.
        Returns `True` only if all selection cuts pass.
        """
        mhh = entries['HHKin_mass']

        # When one only has 0 or 1 bjet th HH mass is not well defined,
        # and a value of -1 is assigned. One thus has to remove the cut below
        # when considering events with less than 2 b-jets.
        if mhh < 1 and bjets_cut:
            return False

        pairtype    = entries['pairType']
        dau1_eleiso = entries['dau1_eleMVAiso']
        dau1_muiso  = entries['dau1_iso']
        dau1_tauiso = entries['dau1_deepTauVsJet']
        dau2_tauiso = entries['dau2_deepTauVsJet']

        # third lepton veto
        nleps = entries['nleps']
        if nleps > 0 and lepton_veto:
            return False

        # require at least two b jet candidates
        nbjetscand = entries['nbjetscand']
        if nbjetscand <= 1 and bjets_cut:
            return False

        # Loose / Medium / Tight
        iso_allowed = { 'dau1_ele': 1., 'dau1_mu': 0.15,
                        'dau1_tau': 5., 'dau2_tau': 5. }
        if any(x not in iso_allowed for x in iso_cuts.keys()):
            mes = 'At least one of the keys is not allowed. '
            mes += 'Keys introduced: {}.'.format(iso_cuts.keys())
            raise ValueError(mes)

        # setting to the defaults in case the user did not specify the values
        for k, v in iso_allowed.items():
            if k not in iso_cuts: iso_cuts[k] = v
        
        bool0 = pairtype==0 and (dau1_muiso >= iso_cuts['dau1_mu'] or
                                 dau2_tauiso < iso_cuts['dau2_tau'])
        bool1 = pairtype==1 and (dau1_eleiso != iso_cuts['dau1_ele'] or
                                 dau2_tauiso < iso_cuts['dau2_tau'])
        bool2 = pairtype==2 and (dau1_tauiso < iso_cuts['dau1_tau'] or
                                 dau2_tauiso < iso_cuts['dau2_tau'])
        if bool0 or bool1 or bool2:
            return False

        #((tauH_SVFIT_mass-116.)*(tauH_SVFIT_mass-116.))/(35.*35.) + ((bH_mass_raw-111.)*(bH_mass_raw-111.))/(45.*45.) <  1.0
        svfit_mass = entries['tauH_SVFIT_mass']
        bh_mass    = entries['bH_mass_raw']

        # mcut_outside = ((svfit_mass-129.)*(svfit_mass-129.) / (53.*53.) +
        #                 (bh_mass-169.)*(bh_mass-169.) / (145.*145.) ) > 1.0
        # if mcut_outside: # inverted elliptical mass cut (-> ttCR)
        #     return False

        return True

def plot(hmet, hnomet, var, channel, sample, category, directory):
    ROOT.gROOT.SetBatch(ROOT.kTRUE)
    ROOT.gStyle.SetOptStat(ROOT.kFALSE)

    c = ROOT.TCanvas('c', '', 600, 400)
    c.Divide(2,1)
    pad1 = ROOT.TPad('pad1', 'pad1', 0, 0.0, 1.0, 1.0)
    pad1.SetFrameLineWidth(3)
    pad1.SetLeftMargin(0.12);
    pad1.SetBottomMargin(0.12);
    pad1.SetTopMargin(0.055);
    pad1.Draw()
    pad1.cd()

    cmsTextSize = 0.05 
    t = ROOT.gPad.GetTopMargin()
    b = ROOT.gPad.GetBottomMargin()
    l = ROOT.gPad.GetLeftMargin()
    r = ROOT.gPad.GetRightMargin()
    
    selBox = ROOT.TLatex(l + 0.03 , 1 - t - 0.02, category)
    selBox.SetNDC()
    selBox.SetTextSize(cmsTextSize+20)
    selBox.SetTextFont(43)
    selBox.SetTextColor(ROOT.kBlack)
    selBox.SetTextAlign(13)
    selBox.Draw('same')
    
    #hmet.SetAxisRange(0,1.1,'Y');
    hmet.GetXaxis().SetTitle(var + ' [GeV]');
    hmet.GetYaxis().SetTitle('a.u.');
    hmet.GetXaxis().SetTitleSize(0.045);
    hmet.GetYaxis().SetTitleSize(0.045);
    hmet.SetLineWidth(2);
    hmet.SetLineColor(8);
    #hmet.SetFillColor(8);

    hnomet.SetLineWidth(2);
    hnomet.SetLineColor(4);
    #hnomet.SetFillColor(4);

    # hmet.Add(hnomet)
    hmet.Scale(1/hmet.Integral())
    hmet.Draw('hist');
    hnomet.Scale(1/hnomet.Integral())
    hnomet.Draw('histsame');

    leg = ROOT.TLegend(0.7, 0.8, 0.90, 0.93)
    leg.SetNColumns(1)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetTextFont(43)
    leg.AddEntry(hmet, 'MET Only')
    leg.AddEntry(hnomet, '+\n'.join(triggers[channel]))
    leg.Draw('same')

    outdir = os.path.join(directory, channel, sample, category)
    utils.create_single_dir(outdir)
    c.Update();
    c.SaveAs( os.path.join(outdir, 'met_' + var + '.png') )

def plot2D(hmet, hnomet, channel, sample, category, directory):
    ROOT.gROOT.SetBatch(ROOT.kTRUE)
    ROOT.gStyle.SetOptStat(ROOT.kFALSE)

    c = ROOT.TCanvas('c', '', 600, 400)
    pad1 = ROOT.TPad('pad1', 'pad1', 0.0, 0.0, .5, 1.0)
    pad1.SetFrameLineWidth(3)
    pad1.SetLeftMargin(0.12);
    pad1.SetBottomMargin(0.12);
    pad1.SetTopMargin(0.055);
    pad1.Draw()
    pad1.cd()

    cmsTextSize = 0.05 
    t = ROOT.gPad.GetTopMargin()
    b = ROOT.gPad.GetBottomMargin()
    l = ROOT.gPad.GetLeftMargin()
    r = ROOT.gPad.GetRightMargin()
    
    selBox = ROOT.TLatex(l + 0.03 , 1 - t - 0.02, category)
    selBox.SetNDC()
    selBox.SetTextSize(cmsTextSize+20)
    selBox.SetTextFont(43)
    selBox.SetTextColor(ROOT.kBlack)
    selBox.SetTextAlign(13)
    selBox.Draw('same')
    
    #hmet.SetAxisRange(0,1.1,'Y');
    hmet.GetXaxis().SetTitle('dau1_pt [GeV]');
    hmet.GetYaxis().SetTitle('dau2_pt [GeV]');
    hmet.GetXaxis().SetTitleSize(0.045);
    hmet.GetYaxis().SetTitleSize(0.045);
    hmet.SetLineWidth(2);
    hmet.SetLineColor(8);
    #hmet.SetFillColor(8);
    
    # hmet.Add(hnomet)
    hmet.Draw('hist colz');

    c.cd()
    pad2 = ROOT.TPad('pad2', 'pad2', 0.5, 0.0, 1.0, 1.0)
    pad2.SetFrameLineWidth(3)
    pad2.SetLeftMargin(0.12);
    pad2.SetBottomMargin(0.12);
    pad2.SetTopMargin(0.055);
    pad2.Draw()
    pad2.cd()

    hnomet.GetXaxis().SetTitle('dau1_pt [GeV]');
    hnomet.GetYaxis().SetTitle('dau2_pt [GeV]');
    hnomet.GetXaxis().SetTitleSize(0.045);
    hnomet.GetYaxis().SetTitleSize(0.045);
    hnomet.SetLineWidth(2);
    hnomet.SetLineColor(8);
    hnomet.Draw('hist colz same');
    #hnomet.SetFillColor(4);

    outdir = os.path.join(directory, channel, sample, category)
    utils.create_single_dir(outdir)
    c.Update();
    c.SaveAs( os.path.join(outdir, 'met_2D.png') )

def test_met(indir, sample, channel, plot_only):
    binning.update({'HHKin_mass': (20, float(sample)-300, float(sample)+300)})
    
    outname = 'met_{}_{}.root'.format(sample, channel)
    full_sample = 'GluGluToBulkGravitonToHHTo2B2Tau_M-' + sample + '_'
    
    if not plot_only:
        t_in = ROOT.TChain('HTauTauTree');
        glob_files = glob.glob( os.path.join(indir, full_sample, 'output_*.root') )
        for f in glob_files:
            t_in.Add(f)
     
        hMET, hNoMET = ({} for _ in range(2))
        for v in tuple(variables):
            hMET[v], hNoMET[v] = ({} for _ in range(2))
            for cat in categories:
                hMET[v][cat] = ROOT.TH1D('hMET'+v+'_'+cat, '', *binning[v])
                hNoMET[v][cat] = ROOT.TH1D('hNoMET'+v+'_'+cat, '', *binning[v])

        hMET_2D, hNoMET_2D = ({} for _ in range(2))
        for cat in categories:
            hMET_2D[cat] = ROOT.TH2D('hMET_2D_'+cat, '', *binning['dau1_pt'], *binning['dau2_pt'])
            hNoMET_2D[cat] = ROOT.TH2D('hNoMET_2D_'+cat, '', *binning['dau1_pt'], *binning['dau2_pt'])        
     
        t_in.SetBranchStatus('*', 0)
        _entries = ('triggerbit', 'bjet1_bID_deepFlavor', 'bjet2_bID_deepFlavor', 'isBoosted',
                    'isVBF', 'VBFjj_mass', 'VBFjj_deltaEta', 'PUReweight', 'lumi', 'IdAndIsoSF_deep_pt',
                    'HHKin_mass', 'pairType', 'dau1_eleMVAiso', 'dau1_iso', 'dau1_deepTauVsJet', 'dau2_deepTauVsJet',
                    'nleps', 'nbjetscand', 'tauH_SVFIT_mass', 'bH_mass_raw',)
        _entries += tuple(variables)
        for ientry in _entries:
            t_in.SetBranchStatus(ientry, 1)
     
        for entry in t_in:
            # this is slow: do it once only
            entries = utils.dot_dict({x: getattr(entry, x) for x in _entries})
            
            # mcweight   = entries.MC_weight
            pureweight = entries.PUReweight
            lumi       = entries.lumi
            idandiso   = entries.IdAndIsoSF_deep_pt
            
            #if utils.is_nan(mcweight)  : mcweight=1
            if utils.is_nan(pureweight) : pureweight=1
            if utils.is_nan(lumi)       : lumi=1
            if utils.is_nan(idandiso)   : idandiso=1
     
            evt_weight = pureweight*lumi*idandiso
            if utils.is_nan(evt_weight):
                evt_weight = 1

            for v in variables:
                
                if utils.is_channel_consistent(channel, entries.pairType):
     
                    if not sel_cuts(entries, lepton_veto=True, bjets_cut=True):
                        continue
     
                    for cat in categories:
                        if sel_category(entries, cat):

                            if pass_triggers(triggers[channel], entries.triggerbit):
                                hNoMET[v][cat].Fill(entries[v], evt_weight)
                                if v == variables[0]:
                                    hNoMET_2D[cat].Fill(entries['dau1_pt'], entries['dau2_pt'], evt_weight)
     
                            if (pass_triggers(('METNoMu120',), entries.triggerbit) and
                                not pass_triggers(triggers[channel], entries.triggerbit)):
                                hMET[v][cat].Fill(entries[v], evt_weight)
                                if v == variables[0]:
                                    hMET_2D[cat].Fill(entries['dau1_pt'], entries['dau2_pt'], evt_weight)

     
        f_out = ROOT.TFile(outname, 'RECREATE')
        f_out.cd()
        for v in variables:
            for cat in categories:
                hMET[v][cat].Write('hMET_' + v + '_' + cat)
                hNoMET[v][cat].Write('hNoMET_' + v + '_' + cat)
                if v == variables[0]:
                    hMET_2D[cat].Write('hMET_2D_' + cat)
                    hNoMET_2D[cat].Write('hNoMET_2D_' + cat)
        f_out.Close()
        print('Raw histograms saved in {}.'.format(outname) )


    f_in = ROOT.TFile(outname, 'READ')
    f_in.cd()
    from_directory = 'MET_Histograms'
    for v in variables:
        for cat in categories:
            hMET = f_in.Get('hMET_' + v + '_' + cat)
            hNoMET = f_in.Get('hNoMET_' + v + '_' + cat)
            plot(hMET, hNoMET, v, channel, sample, cat, from_directory)
            if v == variables[0]:
                hMET_2D = f_in.Get('hMET_2D_' + cat)
                hNoMET_2D = f_in.Get('hNoMET_2D_' + cat)
                plot2D(hMET_2D, hNoMET_2D, channel, sample, cat, from_directory)
                
    f_in.Close()

    from distutils.dir_util import copy_tree
    to_directory = '/eos/user/b/bfontana/www/TriggerScaleFactors/{}'.format(from_directory)
    copy_tree(from_directory, to_directory)

if __name__ == '__main__':
    triggers = {'etau': ('Ele32', 'EleIsoTauCustom'),
                'mutau': ('IsoMu24', 'IsoMuIsoTauCustom'),
                'tautau': ('IsoDoubleTauCustom',)}
    binning = {'metnomu_et': (20, 0, 450),
               'dau1_pt': (30, 0, 500),
               'dau1_eta': (20, -2.5, 2.5),
               'dau1_iso': (24, 0, 8),
               'dau2_iso': (20, 0.8, 1.),
               'dau2_pt': (30, 0, 500),
               'dau2_eta': (20, -2.5, 2.5),
               'ditau_deltaR': (25, 0, 5),
               'dib_deltaR': (25, 0, 5),
               'bH_pt': (20, 70, 400),
               'bH_mass': (30, 0, 300),
               'tauH_mass': (30, 0, 200),
               'tauH_pt': (30, 0, 400),
               'tauH_SVFIT_mass': (30, 0, 300),
               'tauH_SVFIT_pt': (30, 0, 500),
               'bjet1_pt': (20, 20, 200),
               }
    variables = tuple(binning.keys()) + ('HHKin_mass',)
    categories = ('s1b1jresolvedMcut', 's2b0jresolvedMcut', 'sboostedLLMcut')
    # massLimits = {'250': (0,500),
    #               '260': (0,500),
    #               '270': (0,500),
    #               '280': (0,500),
    #               '300': (0,500),
    #               '320': (0,500),
    #               '350': (0,500),
    #               '400': (0,500),
    #               '450': (0,500),
    #               '500': (0,500),
    #               '550': (0,500),
    #               '600': (0,500),
    #               '650': (0,500),
    #               '700': (0,500),
    #               '750': (0,500),
    #               '800': (0,500),
    #               '850': (300,1200),
    #               '900': (0,500),
    #               '1000': (0,500),
    #               '1250': (0,500),
    #               '1500': (0,500),
    #               '1750': (0,500),
    #               '2000': (0,500),
    #               '2500': (0,500),
    #               '3000': (0,500)}
    
    # Parse input arguments
    parser = argparse.ArgumentParser(description='Producer trigger histograms.')

    parser.add_argument('--indir', dest='indir', required=True, type=str,
                        help='Full path of ROOT input file')
    parser.add_argument('--samples', dest='samples', required=True, nargs='+', type=str,
                        help='Full path of ROOT input file')
    parser.add_argument('--channel', dest='channel', required=True, type=str,  
                        help='Select the channel over which the workflow will be run.' )
    parser.add_argument('--plot_only', dest='plot_only', help='use log scale', action='store_true')
    args = utils.parse_args(parser)

    # pool = multiprocessing.Pool(processes=4)    
    # pool.starmap(test_met, zip(it.repeat(args.indir), args.samples,
    #                            it.repeat(args.channel), it.repeat(args.plot_only)))
    
    for sample in args.samples:
        print('Processing sample {}'.format(sample))
        test_met(args.indir, sample, args.channel, args.plot_only)
