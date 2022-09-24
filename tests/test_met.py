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

    if category == 'baseline':
        specific = True
    elif category == 's1b1jresolvedMcut':
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

def set_plot_definitions():
    ROOT.gROOT.SetBatch(ROOT.kTRUE)
    ROOT.gStyle.SetOptStat(ROOT.kFALSE)
    ret = {'BoxTextSize': 50,
           'BoxTextFont': 43,
           'BoxTextColor': ROOT.kBlack,
           'XTitleSize': 0.045,
           'YTitleSize': 0.045,
           'LineWidth': 2,
           }
    return ret
    
def plot(hmet, hnomet, hmetwithcut, var, channel, sample, category, directory):
    defs = set_plot_definitions()
    c = ROOT.TCanvas('c', '', 600, 400)
    c.cd()
    
    pad = ROOT.TPad('pad', 'pad', 0., 0., 1., 1.)
    pad.SetFrameLineWidth(3)

    selBox = ROOT.TLatex(0.5, 0.5, category)
    selBox.SetNDC()
    selBox.SetTextSize(defs['BoxTextSize'])
    selBox.SetTextFont(defs['BoxTextFont'])
    selBox.SetTextColor(defs['BoxTextColor'])
    selBox.SetTextAlign(13)
    selBox.Draw('same')
    
    hmet.GetXaxis().SetTitleSize(defs['XTitleSize']);
    hmet.GetXaxis().SetTitle(var + ' [GeV]');
    hmet.GetYaxis().SetTitleSize(defs['YTitleSize']);
    hmet.GetYaxis().SetTitle('a.u.');
    hmet.SetLineWidth(defs['LineWidth']);
    hmet.SetLineColor(8);

    hnomet.SetLineWidth(2);
    hnomet.SetLineColor(4);
    hmetwithcut.SetLineWidth(2);
    hmetwithcut.SetLineColor(2);

    hmet.Scale(1/hmet.Integral())
    hnomet.Scale(1/hnomet.Integral())
    hmetwithcut.Scale(1/hmetwithcut.Integral())

    max_met   = hmet.GetMaximum() + (hmet.GetMaximum()-hmet.GetMinimum())/5.
    max_nomet = hnomet.GetMaximum() + (hnomet.GetMaximum()-hnomet.GetMinimum())/5.
    max_withcut = hmetwithcut.GetMaximum() + (hmetwithcut.GetMaximum()-hmetwithcut.GetMinimum())/5.
    hmet.SetMaximum( max(max_met, max_nomet, max_withcut) )

    hmet.Draw('hist')
    hnomet.Draw('histsame')
    hmetwithcut.Draw('histsame')

    leg = ROOT.TLegend(0.65, 0.75, 0.90, 0.9)
    leg.SetNColumns(1)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetTextFont(43)
    leg.SetTextSize(10)
    leg.AddEntry(hmet, 'MET')
    leg.AddEntry(hmetwithcut, 'Met + Cut')
    leg.AddEntry(hnomet, '+\n'.join(triggers[channel]))
    leg.Draw('same')

    cat_folder = os.path.join(directory, category)
    utils.create_single_dir(cat_folder)
    c.Update();
    c.SaveAs( os.path.join(cat_folder, 'met_' + var + '.png') )
    c.Close()

def plot2D(hmet, hnomet, hmetwithcut, two_vars, channel, sample, category, directory):
    defs = set_plot_definitions()    
    c = ROOT.TCanvas('c', '', 600, 400)
    c.cd()

    pad = ROOT.TPad('pad1', 'pad1', 0., 0., 0.33, 1.)
    pad.SetFrameLineWidth(3)
    pad.SetLeftMargin(0.12);
    pad.SetBottomMargin(0.12);
    pad.SetTopMargin(0.055);
    pad.Draw()
    pad.cd()

    hmet.GetXaxis().SetTitle(two_vars[0])
    hmet.GetYaxis().SetTitle(two_vars[1])
    hmet.GetXaxis().SetTitleSize(0.045);
    hmet.GetYaxis().SetTitleSize(0.045);
    hmet.SetLineWidth(2);
    hmet.SetLineColor(8);
    #hmet.SetFillColor(8);
    
    # hmet.Add(hnomet)
    hmet.Draw('hist colz');

    c.cd()
    pad2 = ROOT.TPad('pad2', 'pad2', 0.33, 0.0, 0.67, 1.0)
    pad2.SetFrameLineWidth(3)
    pad2.SetLeftMargin(0.12);
    pad2.SetBottomMargin(0.12);
    pad2.SetTopMargin(0.055);
    pad2.Draw()
    pad2.cd()

    hnomet.GetXaxis().SetTitle(two_vars[0]);
    hnomet.GetYaxis().SetTitle(two_vars[1]);
    hnomet.GetXaxis().SetTitleSize(0.045);
    hnomet.GetYaxis().SetTitleSize(0.045);
    hnomet.SetLineWidth(2);
    hnomet.SetLineColor(8);
    hnomet.Draw('hist colz same');
    #hnomet.SetFillColor(4);

    c.cd()
    pad3 = ROOT.TPad('pad3', 'pad3', 0.67, 0.0, 1.0, 1.0)
    pad3.SetFrameLineWidth(3)
    pad3.SetLeftMargin(0.12);
    pad3.SetBottomMargin(0.12);
    pad3.SetTopMargin(0.055);
    pad3.Draw()
    pad3.cd()

    hmetwithcut.GetXaxis().SetTitle(two_vars[0]);
    hmetwithcut.GetYaxis().SetTitle(two_vars[1]);
    hmetwithcut.GetXaxis().SetTitleSize(0.045);
    hmetwithcut.GetYaxis().SetTitleSize(0.045);
    hmetwithcut.SetLineWidth(2);
    hmetwithcut.SetLineColor(8);
    hmetwithcut.Draw('hist colz same');
    #hmetwithcut.SetFillColor(4);

    outdir = os.path.join(directory, channel, sample)
    utils.create_single_dir(outdir)
    c.Update();
    c.SaveAs( os.path.join(outdir, 'met_' + '_VS_'.join(two_vars) + '.png') )
    c.Close()

def test_met(indir, sample, channel, plot_only):
    if channel == 'etau' or 'mutau':
        iso1 = (24, 0, 8)
    elif channel == 'tautau':
        iso1 = (20, 0.86, 1.01)
    binning.update({'HHKin_mass': (20, float(sample)-300, float(sample)+300),
                    'dau1_iso': iso1})
    
    outname = 'met_{}_{}.root'.format(sample, channel)
    full_sample = 'GluGluToBulkGravitonToHHTo2B2Tau_M-' + sample + '_'
    
    if not plot_only:
        t_in = ROOT.TChain('HTauTauTree');
        glob_files = glob.glob( os.path.join(indir, full_sample, 'output_*.root') )
        for f in glob_files:
            t_in.Add(f)
     
        hMET, hNoMET, hMETWithCut = ({} for _ in range(3))
        for v in tuple(variables):
            hMET[v], hNoMET[v], hMETWithCut[v] = ({} for _ in range(3))
            for cat in categories:
                hMET[v][cat] = ROOT.TH1D('hMET_'+v+'_'+cat, '', *binning[v])
                hNoMET[v][cat] = ROOT.TH1D('hNoMET_'+v+'_'+cat, '', *binning[v])
                hMETWithCut[v][cat] = ROOT.TH1D('hMETWithCut_'+v+'_'+cat, '', *binning[v])

        hMET_2D, hNoMET_2D, hMETWithCut_2D = ({} for _ in range(3))
        for v in variables_2D:
            hMET_2D[v], hNoMET_2D[v], hMETWithCut_2D[v] = ({} for _ in range(3))
            for cat in categories:
                hMET_2D[v][cat] = ROOT.TH2D('hMET_2D_'+'_'.join(v)+'_'+cat, '', *binning[v[0]], *binning[v[1]])
                hNoMET_2D[v][cat] = ROOT.TH2D('hNoMET_2D_'+'_'.join(v)+'_'+cat, '', *binning[v[0]], *binning[v[1]])
                hMETWithCut_2D[v][cat] = ROOT.TH2D('hMETWithCut_2D_'+'_'.join(v)+'_'+cat, '', *binning[v[0]], *binning[v[1]])        

        t_in.SetBranchStatus('*', 0)
        _entries = ('triggerbit', 'bjet1_bID_deepFlavor', 'bjet2_bID_deepFlavor', 'isBoosted',
                    'isVBF', 'VBFjj_mass', 'VBFjj_deltaEta', 'PUReweight', 'lumi', 'IdAndIsoSF_deep_pt',
                    'pairType', 'dau1_eleMVAiso', 'dau1_iso', 'dau1_deepTauVsJet', 'dau2_deepTauVsJet',
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

            if utils.is_channel_consistent(channel, entries.pairType):
                if not sel_cuts(entries, lepton_veto=True, bjets_cut=True):
                    continue

                for v in variables:
                    for cat in categories:
                        if sel_category(entries, cat):

                            # passes the OR of the trigger baseline (not including METNoMu120 trigger)
                            if pass_triggers(triggers[channel], entries.triggerbit):
                                hNoMET[v][cat].Fill(entries[v], evt_weight)

                            # passes the METNoMu120 trigger and does *not* pass the OR of the baseline
                            if (pass_triggers(('METNoMu120',), entries.triggerbit) and
                                not pass_triggers(triggers[channel], entries.triggerbit)):
                                hMET[v][cat].Fill(entries[v], evt_weight)
                                if entries.metnomu_et > met_cut:
                                    hMETWithCut[v][cat].Fill(entries[v], evt_weight)

                for v in variables_2D:
                    for cat in categories:
                        if sel_category(entries, cat):

                            # passes the OR of the trigger baseline (not including METNoMu120 trigger)
                            if pass_triggers(triggers[channel], entries.triggerbit):
                                hNoMET_2D[v][cat].Fill(entries[v[0]], entries[v[1]], evt_weight)

                            # passes the METNoMu120 trigger and does *not* pass the OR of the baseline
                            if (pass_triggers(('METNoMu120',), entries.triggerbit) and
                                not pass_triggers(triggers[channel], entries.triggerbit)):
                                hMET_2D[v][cat].Fill(entries[v[0]], entries[v[1]], evt_weight)
                                if entries.metnomu_et > met_cut:
                                    hMETWithCut_2D[v][cat].Fill(entries[v[0]], entries[v[1]], evt_weight)


        f_out = ROOT.TFile(outname, 'RECREATE')
        f_out.cd()
        for cat in categories:
            for v in variables:
                hMET[v][cat].Write('hMET_' + v + '_' + cat)
                hNoMET[v][cat].Write('hNoMET_' + v + '_' + cat)
                hMETWithCut[v][cat].Write('hMETWithCut_' + v + '_' + cat)
            for v in variables_2D:
                hMET_2D[v][cat].Write('hMET_2D_' + '_'.join(v)+'_'+ cat)
                hNoMET_2D[v][cat].Write('hNoMET_2D_' + '_'.join(v)+'_'+ cat)
                hMETWithCut_2D[v][cat].Write('hMETWithCut_2D_' + '_'.join(v)+'_'+ cat)
        f_out.Close()
        print('Raw histograms saved in {}.'.format(outname), flush=True)


    f_in = ROOT.TFile(outname, 'READ')
    f_in.cd()
    from_directory = os.path.join('MET_Histograms', channel, sample)
    for cat in categories:
        for v in variables:
            hMET = f_in.Get('hMET_' + v + '_' + cat)
            hNoMET = f_in.Get('hNoMET_' + v + '_' + cat)
            hMETWithCut = f_in.Get('hMETWithCut_' + v + '_' + cat)
            plot(hMET, hNoMET, hMETWithCut, v, channel, sample, cat, from_directory)
        for v in variables_2D:
            hMET_2D = f_in.Get('hMET_2D_' + '_'.join(v)+'_'+ cat)
            hNoMET_2D = f_in.Get('hNoMET_2D_' + '_'.join(v)+'_'+ cat)
            hMETWithCut_2D = f_in.Get('hMETWithCut_2D_' + '_'.join(v)+'_'+ cat)
            plot2D(hMET_2D, hNoMET_2D, hMETWithCut_2D, v, channel, sample, cat, from_directory)
                
    f_in.Close()

    from distutils.dir_util import copy_tree
    to_directory = '/eos/user/b/bfontana/www/TriggerScaleFactors/{}'.format(from_directory)
    copy_tree(from_directory, to_directory)
    print(from_directory)
    print(to_directory)
    print('Pictures copied to {}.'.format(to_directory), flush=True)

if __name__ == '__main__':
    triggers = {'etau': ('Ele32', 'EleIsoTauCustom'),
                'mutau': ('IsoMu24', 'IsoMuIsoTauCustom'),
                'tautau': ('IsoDoubleTauCustom',)}
    binning = {'metnomu_et': (20, 0, 450),
               'dau1_pt': (30, 0, 450),
               'dau1_eta': (20, -2.5, 2.5),
               'dau2_iso': (20, 0.87, 1.01),
               'dau2_pt': (30, 0, 350),
               'dau2_eta': (20, -2.5, 2.5),
               'ditau_deltaR': (25, 0, 1.7),
               'dib_deltaR': (25, 0, 2.5),
               'bH_pt': (20, 70, 600),
               'bH_mass': (30, 0, 280),
               'tauH_mass': (30, 0, 170),
               'tauH_pt': (30, 0, 500),
               'tauH_SVFIT_mass': (50, 0, 250),
               'tauH_SVFIT_pt': (30, 100, 600),
               'bjet1_pt': (16, 20, 500),
               'bjet2_pt': (16, 20, 500),
               'bjet1_eta': (20, -2.5, 2.5),
               'bjet2_eta': (20, -2.5, 2.5),
               }
    variables = tuple(binning.keys()) + ('HHKin_mass', 'dau1_iso')
    variables_2D = (('dau1_pt', 'dau2_pt'), ('dau1_iso', 'dau2_iso'))
    categories = ('baseline', 's1b1jresolvedMcut', 's2b0jresolvedMcut', 'sboostedLLMcut')
    met_cut = 200
    
    # Parse input arguments
    parser = argparse.ArgumentParser(description='Producer trigger histograms.')

    parser.add_argument('--indir', dest='indir', required=True, type=str,
                        help='Full path of ROOT input file')
    parser.add_argument('--samples', dest='samples', required=True, nargs='+', type=str,
                        help='Full path of ROOT input file')
    parser.add_argument('--channel', dest='channel', required=True, type=str,  
                        help='Select the channel over which the workflow will be run.' )
    parser.add_argument('--plot_only', dest='plot_only', action='store_true',
                        help='Reuse previously produced data for quick plot changes.')
    args = utils.parse_args(parser)

    # pool = multiprocessing.Pool(processes=4)    
    # pool.starmap(test_met, zip(it.repeat(args.indir), args.samples,
    #                            it.repeat(args.channel), it.repeat(args.plot_only)))
    
    for sample in args.samples:
        print('Processing sample {}'.format(sample))
        test_met(args.indir, sample, args.channel, args.plot_only)
