# coding: utf-8

_all_ = [ 'test_trigger_regions' ]

import os
import sys
parent_dir = os.path.abspath(__file__ + 2 * '/..')
sys.path.insert(0, parent_dir)
import argparse
import glob
import multiprocessing
import itertools as it
import csv

import inclusion
from inclusion import selection
from inclusion.config import main
from inclusion.utils import utils

import ROOT

def get_outname(suffix, mode, cut='', ext=''):
    pref = {'met': 'MET', 'tau': 'Tau', 'met_tau': 'MET_and_Tau'}
    assert mode in list(pref.keys())
    if ext == 'csv':
        s = 'counts_{}_{}'.format(suffix, pref[mode])
        if cut:
            s += '_{}'.format(cut)
        s += '/table.csv'
    elif ext == 'root':
        utils.create_single_dir('data')
        s = 'data/data_{}_{}'.format(suffix, pref[mode])
        if cut:
            s += '_{}'
        s += '.{}'.format(ext)
    else:
        raise ValueError('The {} extension is not supported.'.format(ext))
    return s

def set_plot_definitions():
    ROOT.gROOT.SetBatch(ROOT.kTRUE)
    ROOT.gStyle.SetOptStat(ROOT.kFALSE)
    ret = {'BoxTextSize'    : 50,
           'BoxTextFont'    : 43,
           'BoxTextColor'   : ROOT.kBlack,
           'XTitleSize'     : 0.045,
           'YTitleSize'     : 0.045,
           'LineWidth'      : 2,
           'FrameLineWidth' : 1,
           }
    return ret
        
def plot2D(histo, two_vars, channel, sample, category, directory, region):
    histo = histo.Clone('histo')
    defs = set_plot_definitions()    

    c = ROOT.TCanvas('c', '', 800, 800)
    c.cd()
    
    # pad1 = ROOT.TPad('pad1', 'pad1', 0., 0., 1., 1.)
    # pad1.cd()
    # pad1.SetFrameLineWidth(defs['FrameLineWidth'])
    # pad1.SetLeftMargin(0.15);
    # pad1.SetRightMargin(0.0);
    # pad1.SetBottomMargin(0.08);
    # pad1.SetTopMargin(0.055);
    # pad1.Draw()

    histo.GetXaxis().SetTitle(two_vars[0])
    histo.GetXaxis().SetTitleSize(0.045)
    histo.GetYaxis().SetTitle(two_vars[1])
    histo.GetYaxis().SetTitleSize(0.045)
    histo.Draw('colz')

    cat_folder = os.path.join(directory, sample, category)
    utils.create_single_dir(cat_folder)
    c.Update();
    for ext in ('png', 'pdf'):
        c.SaveAs( os.path.join(cat_folder, 'reg' + region + '_' + '_VS_'.join(two_vars) + '.' + ext) )
    c.Close()

def test_trigger_regions(indir, sample, channel):
    outname = get_outname(suffix=sample+'_'+channel, mode='met', ext='root')

    if channel == 'etau' or channel == 'mutau':
        iso1 = (24, 0, 8)
    elif channel == 'tautau':
        iso1 = binning['dau2_iso']
    binning.update({'HHKin_mass': (20, float(sample)-300, float(sample)+300),
                    'dau1_iso': iso1})
    
    full_sample = 'GluGluToBulkGravitonToHHTo2B2Tau_M-' + sample + '_'
    
    t_in = ROOT.TChain('HTauTauTree')
    glob_files = glob.glob( os.path.join(indir, full_sample, 'output_*.root') )
    for f in glob_files:
        t_in.Add(f)

    hBase, hMET, hTau = ({} for _ in range(3)) # one trigger
    hBase_MET, hBase_Tau, hMET_Tau = ({} for _ in range(3)) # intersection of two triggers
    hBase_MET_Tau = {} # intersection of three triggers
    for ireg in range(nregions):
        hBase[ireg], hMET[ireg], hTau[ireg] = ({} for _ in range(3))
        hBase_MET[ireg], hBase_Tau[ireg], hMET_Tau[ireg] = ({} for _ in range(3))
        hBase_MET_Tau[ireg] = {}
        
        for v in tuple(variables_2D):
            hBase[ireg][v], hMET[ireg][v], hTau[ireg][v] = ({} for _ in range(3))
            hBase_MET[ireg][v], hBase_Tau[ireg][v], hMET_Tau[ireg][v] = ({} for _ in range(3))
            hBase_MET_Tau[ireg][v] = {}

            for cat in categories:
                suff = lambda x: x + '_' + str(ireg) + '_' + v[0] + '_VS_' + v[1] + '_' + cat
                hopt = (*binning[v[0]], *binning[v[1]])
                hBase[ireg][v][cat] = ROOT.TH2D(suff('hBase'), '', *hopt)
                hMET[ireg][v][cat]  = ROOT.TH2D(suff('hMET'), '', *hopt)
                hTau[ireg][v][cat]  = ROOT.TH2D(suff('hTau'), '', *hopt)
                
                hBase_MET[ireg][v][cat] = ROOT.TH2D(suff('hBase_MET'), '', *hopt)
                hBase_Tau[ireg][v][cat] = ROOT.TH2D(suff('hBase_Tau'), '', *hopt)
                hMET_Tau[ireg][v][cat]  = ROOT.TH2D(suff('hMET_Tau') , '', *hopt)
                
                hBase_MET_Tau[ireg][v][cat] = ROOT.TH2D(suff('hBase_MET_Tau'), '', *hopt)

    t_in.SetBranchStatus('*', 0)
    _entries = ('triggerbit', 'RunNumber', 'isLeptrigger',
                'bjet1_bID_deepFlavor', 'bjet2_bID_deepFlavor', 'isBoosted',
                'isVBF', 'VBFjj_mass', 'VBFjj_deltaEta', 'PUReweight', 'lumi', 'IdAndIsoSF_deep_pt',
                'pairType', 'dau1_eleMVAiso', 'dau1_iso', 'dau1_deepTauVsJet', 'dau2_deepTauVsJet',
                'nleps', 'nbjetscand', 'tauH_SVFIT_mass', 'bH_mass_raw',)
    _entries += tuple(variables)
    for ientry in _entries:
        t_in.SetBranchStatus(ientry, 1)
  
    for entry in t_in:
        # this is slow: do it once only
        entries = utils.dot_dict({x: getattr(entry, x) for x in _entries})
        sel = selection.EventSelection(entries, isdata=False, configuration=None)
        
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
            if not sel.selection_cuts(lepton_veto=True, bjets_cut=True,
                                      standard_mass_cut=True, invert_mass_cut=False):
                continue

            for vx, vy in variables_2D:
                for cat in categories:
                    if sel.sel_category(cat) and entries.ditau_deltaR > 0.5:

                        if eval(met_region):
                            ireg = 1
                        elif eval(tau_region):
                            ireg = 2
                        else: # trigger baseline
                            ireg = 0

                        met_cut_expr = entries.metnomu_et > met_cut
                        tau_cut_expr = ((entries.dau1_pt > tau_cut and args.channel=='tautau') or
                                        (entries.dau2_pt > tau_cut and args.channel!='tautau'))
                        pass_trg = sel.pass_triggers(triggers[channel]) and entries.isLeptrigger
                        pass_met = sel.pass_triggers(('METNoMu120',))
                        pass_tau = sel.pass_triggers(('IsoTau180',))
                        
                        if pass_trg:
                            hBase[ireg][v][cat].Fill(entries[vx], entries[vy], evt_weight)

                        if pass_met:
                            hMET[ireg][v][cat].Fill(entries[vx], entries[vy], evt_weight)
                            
                        if pass_tau:
                            hTau[ireg][v][cat].Fill(entries[vx], entries[vy], evt_weight)

                        if pass_trg and pass_met:
                            hBase_MET[ireg][v][cat].Fill(entries[vx], entries[vy], evt_weight)

                        if pass_trg and pass_tau:
                            hBase_Tau[ireg][v][cat].Fill(entries[vx], entries[vy], evt_weight)

                        if pass_met and pass_tau:
                            hMET_Tau[ireg][v][cat].Fill(entries[vx], entries[vy], evt_weight)

                        if pass_trg and pass_met and pass_tau:
                            hBase_MET_Tau[ireg][v][cat].Fill(entries[vx], entries[vy], evt_weight)

    f_out = ROOT.TFile(outname, 'RECREATE')
    f_out.cd()
    for ireg in range(nregions):
        for cat in categories:
            for v in variables_2D:
                suff = lambda x: x + '_' + str(ireg) + '_' + v[0] + '_VS_' + v[1] + '_' + cat
                hBase[ireg][v][cat].Write(suff('hBase'))
                hMET[ireg][v][cat].Write(suff('hMET'))
                hTau[ireg][v][cat].Write(suff('hTau'))

                hBase_MET[ireg][v][cat].Write(suff('hBase_MET'))
                hBase_Tau[ireg][v][cat].Write(suff('hBase_Tau'))
                hMET_Tau[ireg][v][cat].Write(suff('hMET_Tau'))

                hBase_MET_Tau[ireg][v][cat].Write(suff('hBase_MET_Tau'))
                
    f_out.Close()
    print('Raw histograms saved in {}.'.format(outname), flush=True)

if __name__ == '__main__':
    triggers = {'etau': ('Ele32', 'EleIsoTauCustom'),
                'mutau': ('IsoMu24', 'IsoMuIsoTauCustom'),
                'tautau': ('IsoDoubleTauCustom',)}
    binning = {'metnomu_et': (20, 0, 450),
               'dau1_pt': (30, 0, 450),
               'dau1_eta': (20, -2.5, 2.5),
               'dau2_iso': (20, 0.88, 1.01),
               'dau2_pt': (30, 0, 350),
               'dau2_eta': (20, -2.5, 2.5),
               'ditau_deltaR': (30, 0.3, 1.3),
               'dib_deltaR': (25, 0, 2.5),
               'bH_pt': (20, 70, 600),
               'bH_mass': (30, 0, 280),
               'tauH_mass': (30, 0, 170),
               'tauH_pt': (30, 0, 500),
               'tauH_SVFIT_mass': (30, 0, 250),
               'tauH_SVFIT_pt': (20, 200, 650),
               'bjet1_pt': (25, 10, 600),
               'bjet2_pt': (25, 10, 550),
               'bjet1_eta': (20, -2.5, 2.5),
               'bjet2_eta': (20, -2.5, 2.5),
               }
    variables = tuple(binning.keys()) + ('HHKin_mass', 'dau1_iso')
    #variables = tuple(('dau1_pt', 'dau2_pt', 'metnomu_et', 'HHKin_mass', 'dau1_iso', 'dau2_iso'))
    variables_2D = (('dau1_pt', 'dau2_pt'),)#('dau1_iso', 'dau2_iso'))

    #categories = ('baseline', 's1b1jresolvedMcut', 's2b0jresolvedMcut', 'sboostedLLMcut')
    categories = ('baseline',)
    met_cut = 200
    tau_cut = 190
    
    # Parse input arguments
    parser = argparse.ArgumentParser(description='Producer trigger histograms.')

    parser.add_argument('--indir', required=True, type=str,
                        help='Full path of ROOT input file')
    parser.add_argument('--samples', required=True, nargs='+', type=str,
                        help='Full path of ROOT input file')
    parser.add_argument('--channel', required=True, type=str,  
                        help='Select the channel over which the workflow will be run.' )
    parser.add_argument('--plot', action='store_true',
                        help='Reuse previously produced data for quick plot changes.')
    parser.add_argument('--copy', action='store_true',
                        help='Do not copy the outputs to EOS at the end.')
    parser.add_argument('--sequential', action='store_true',
                        help='Do not use the multiprocess package.')
    args = utils.parse_args(parser)

    region_cut = '160'
    main_dir = 'Region_' + region_cut

    nregions = 3
    met_region = 'entries.dau1_pt < 40 and entries.dau2_pt < {}'.format(region_cut)
    tau_region = 'entries.dau1_pt < 40 and entries.dau2_pt >= {}'.format(region_cut)

    #### run major function ###
    if args.sequential:
        if not args.plot:
            for sample in args.samples:
                test_trigger_regions(args.indir, sample, args.channel)
    else:
        if not args.plot:
            pool = multiprocessing.Pool(processes=6)
            pool.starmap(test_trigger_regions,
                         zip(it.repeat(args.indir), args.samples, it.repeat(args.channel)))
    ###########################

    from_directory = os.path.join(main_dir, args.channel)
    for sample in args.samples:
        outname = get_outname(suffix=sample+'_'+args.channel, mode='met', ext='root')
        f_in = ROOT.TFile(outname, 'READ')
        f_in.cd()

        # write csv header, one per category
        out_counts = []
        for cat in categories:
            out_counts.append( os.path.join(from_directory, sample, cat, 'counts') )
            utils.create_single_dir(out_counts[-1])
            with open(os.path.join(out_counts[-1], 'table.csv'), 'w') as f:
                reader = csv.writer(f, delimiter=',', quotechar='|')
                header_row = ['Region', 'Base', 'MET', 'Tau', 'Base && MET',
                              'Base && Tau', 'MET && Tau', 'Base && MET && Tau']
                reader.writerow(header_row)

        # plot histograms and fill CSV with histogram integrals
        for ireg in range(nregions):
            for cat in categories:
                for v in variables_2D:
                    suff = lambda x : x+'_'+str(ireg)+'_'+v[0]+'_VS_'+v[1]+'_'+cat
                    
                    hBase = f_in.Get(suff('hBase'))
                    hMET = f_in.Get(suff('hMET'))
                    hTau = f_in.Get(suff('hTau'))
 
                    hBase_MET = f_in.Get(suff('hBase_MET'))
                    hBase_Tau = f_in.Get(suff('hBase_Tau'))
                    hMET_Tau = f_in.Get(suff('hMET_Tau'))
 
                    hBase_MET_Tau = f_in.Get(suff('hBase_MET_Tau'))
 
                    plot2D(hMET, v, args.channel, sample, cat, from_directory, str(ireg))
 
                nbinsx, nbinsy = hBase.GetNbinsX(), hBase.GetNbinsY()
                cBase         = round(hBase.Integral(0, nbinsx+1, 0, nbinsy+1), 3)
                cMET          = round(hMET.Integral(0, nbinsx+1, 0, nbinsy+1), 3)
                cTau          = round(hTau.Integral(0, nbinsx+1, 0, nbinsy+1), 3)
                cBase_MET     = round(hBase_MET.Integral(0, nbinsx+1, 0, nbinsy+1), 3)
                cBase_Tau     = round(hBase_Tau.Integral(0, nbinsx+1, 0, nbinsy+1), 3)
                cMET_Tau      = round(hMET_Tau.Integral(0, nbinsx+1, 0, nbinsy+1), 3)
                cBase_MET_Tau = round(hBase_MET_Tau.Integral(0, nbinsx+1, 0, nbinsy+1), 3)

                with open(os.path.join(out_counts[categories.index(cat)], 'table.csv'), 'a') as f:
                    reader = csv.writer(f, delimiter=',', quotechar='|')
                    row = [ireg, cBase, cMET, cTau, cBase_MET, cBase_Tau, cMET_Tau, cBase_MET_Tau]
                    reader.writerow(row)
                        
        f_in.Close()

    if args.copy:
        import subprocess
        to_directory = os.path.join('/eos/user/b/bfontana/www/TriggerScaleFactors', main_dir)
        to_directory = os.path.join(to_directory, args.channel)
     
        for ireg in range(nregions):
            for cat in categories:
                folder_name = 'counts_total_' + cat + '_' + m
                folder_to = os.path.join(to_directory, folder_name)
                utils.create_single_dir(folder_to)
                counts_files = os.path.join(from_directory, folder_name, 'table.csv')
                print('Copying: {}\t\t--->\t{}'.format(counts_files, folder_to), flush=True)
                subprocess.run(['rsync', '-ah', counts_files, os.path.join(folder_to, 'table.csv')])
     
        for sample in args.samples:
            sample_from = os.path.join(from_directory, sample)
            print('Copying: {}\t\t--->\t{}'.format(sample_from, to_directory), flush=True)
            subprocess.run(['rsync', '-ah', sample_from, to_directory])
     
    print('Done.')
