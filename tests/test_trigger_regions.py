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
from collections import defaultdict

import inclusion
from inclusion import selection
from inclusion.config import main
from inclusion.utils import utils

import ROOT

import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn3_circles

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
        
def plot2D(histo, two_vars, channel, sample, suffix, category, directory, region):
    histo = histo.Clone('histo_clone')
    defs = set_plot_definitions()    

    c = ROOT.TCanvas('c', '', 800, 800)
    c.cd()
    
    pad1 = ROOT.TPad('pad1', 'pad1', 0., 0., 1., 1.)
    pad1.SetFrameLineWidth(defs['FrameLineWidth'])
    pad1.SetLeftMargin(0.12);
    pad1.SetRightMargin(0.1);
    pad1.SetBottomMargin(0.12);
    pad1.SetTopMargin(0.04);
    pad1.Draw()
    pad1.cd()

    histo.GetXaxis().SetTitle(two_vars[0])
    histo.GetXaxis().SetTitleSize(0.04)
    histo.GetYaxis().SetTitle(two_vars[1])
    histo.GetYaxis().SetTitleSize(0.04)
    histo.Draw('colz')

    c.cd()
    cat_folder = os.path.join(directory, sample, category)
    utils.create_single_dir(cat_folder)
    c.Update();
    for ext in ('png', 'pdf'):
        c.SaveAs( os.path.join(cat_folder, suffix + '_' + 'reg' + region + '_' + '_VS_'.join(two_vars) + '.' + ext) )
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

    hBase, hMET, hTau = (defaultdict(lambda: defaultdict(dict)) for _ in range(3)) # one trigger
    hBase_MET, hBase_Tau, hMET_Tau = (defaultdict(lambda: defaultdict(dict)) for _ in range(3)) # intersection of two triggers
    hNoBase_MET, hNoBase_NoMET_Tau = (defaultdict(lambda: defaultdict(dict)) for _ in range(2)) # negations
    hNoBase_Tau, hNoBase_MET_NoTau = (defaultdict(lambda: defaultdict(dict)) for _ in range(2)) # negations
    hBase_MET_Tau = defaultdict(lambda: defaultdict(dict)) # intersection of three triggers
    for ireg in range(nregions):
        for v in tuple(variables_2D):
            for cat in categories:
                suff = lambda x: x + '_' + str(ireg) + '_' + v[0] + '_VS_' + v[1] + '_' + cat
                hopt = (*binning[v[0]], *binning[v[1]])

                hBase[ireg][v][cat] = ROOT.TH2D(suff('hBase'), '', *hopt)
                hMET[ireg][v][cat]  = ROOT.TH2D(suff('hMET'), '', *hopt)
                hTau[ireg][v][cat]  = ROOT.TH2D(suff('hTau'), '', *hopt)
                
                hBase_MET[ireg][v][cat] = ROOT.TH2D(suff('hBase_MET'), '', *hopt)
                hBase_Tau[ireg][v][cat] = ROOT.TH2D(suff('hBase_Tau'), '', *hopt)
                hMET_Tau[ireg][v][cat]  = ROOT.TH2D(suff('hMET_Tau') , '', *hopt)

                hNoBase_MET[ireg][v][cat]       = ROOT.TH2D(suff('hNoBase_MET'), '', *hopt)
                hNoBase_NoMET_Tau[ireg][v][cat] = ROOT.TH2D(suff('hNoBase_NoMET_Tau'), '', *hopt)
                hNoBase_Tau[ireg][v][cat]       = ROOT.TH2D(suff('hNoBase_Tau'), '', *hopt)
                hNoBase_MET_NoTau[ireg][v][cat] = ROOT.TH2D(suff('hNoBase_MET_NoTau'), '', *hopt)

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
                        pass_met = sel.pass_triggers(('METNoMu120',)) and met_cut_expr
                        pass_tau = sel.pass_triggers(('IsoTau180',)) and tau_cut_expr
                        
                        if pass_trg:
                            hBase[ireg][v][cat].Fill(entries[vx], entries[vy], evt_weight)

                        if pass_met:
                            hMET[ireg][v][cat].Fill(entries[vx], entries[vy], evt_weight)
                            
                        if pass_tau:
                            hTau[ireg][v][cat].Fill(entries[vx], entries[vy], evt_weight)

                        # intersections with two triggers
                        if pass_trg and pass_met:
                            hBase_MET[ireg][v][cat].Fill(entries[vx], entries[vy], evt_weight)

                        if pass_trg and pass_tau:
                            hBase_Tau[ireg][v][cat].Fill(entries[vx], entries[vy], evt_weight)

                        if pass_met and pass_tau:
                            hMET_Tau[ireg][v][cat].Fill(entries[vx], entries[vy], evt_weight)

                        # negations
                        if not pass_trg and pass_met:
                            hNoBase_MET[ireg][v][cat].Fill(entries[vx], entries[vy], evt_weight)

                        if not pass_trg and not pass_met and pass_tau:
                            hNoBase_NoMET_Tau[ireg][v][cat].Fill(entries[vx], entries[vy], evt_weight)

                        if not pass_trg and pass_tau:
                            hNoBase_Tau[ireg][v][cat].Fill(entries[vx], entries[vy], evt_weight)

                        if not pass_trg and pass_met and not pass_tau:
                            hNoBase_MET_NoTau[ireg][v][cat].Fill(entries[vx], entries[vy], evt_weight)

                        # intersection with three triggers
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

                hNoBase_MET[ireg][v][cat].Write(suff('hNoBase_MET'))
                hNoBase_NoMET_Tau[ireg][v][cat].Write(suff('hNoBase_NoMET_Tau'))
                hNoBase_Tau[ireg][v][cat].Write(suff('hNoBase_Tau'))
                hNoBase_MET_NoTau[ireg][v][cat].Write(suff('hNoBase_MET_NoTau'))

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

    region_cut = '200'
    main_dir = 'Region_' + region_cut

    nregions = 3
    met_region = 'entries.dau2_pt < 40 and entries.dau1_pt < {}'.format(region_cut)
    tau_region = 'entries.dau2_pt < 40 and entries.dau1_pt >= {}'.format(region_cut)

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
                header_row = ['Region', 'Base', 'MET', 'Tau',
                              'Base && MET', 'Base && Tau', 'MET && Tau',
                              '!Base && MET', '!Base && !MET && Tau',
                              '!Base && Tau', '!Base && MET && !Tau',
                              'Base && MET && Tau']
                reader.writerow(header_row)

        # plot histograms and fill CSV with histogram integrals
        for ireg in range(nregions):
            for cat in categories:
                for v in variables_2D:
                    suff = lambda x : x+'_'+str(ireg)+'_'+v[0]+'_VS_'+v[1]+'_'+cat
                    
                    hBase = f_in.Get(suff('hBase'))
                    hMET  = f_in.Get(suff('hMET'))
                    hTau  = f_in.Get(suff('hTau'))
 
                    hBase_MET = f_in.Get(suff('hBase_MET'))
                    hBase_Tau = f_in.Get(suff('hBase_Tau'))
                    hMET_Tau  = f_in.Get(suff('hMET_Tau'))

                    hNoBase_MET       = f_in.Get(suff('hNoBase_MET'))
                    hNoBase_NoMET_Tau = f_in.Get(suff('hNoBase_NoMET_Tau'))
                    hNoBase_Tau       = f_in.Get(suff('hNoBase_Tau'))
                    hNoBase_MET_NoTau = f_in.Get(suff('hNoBase_MET_NoTau'))
 
                    hBase_MET_Tau = f_in.Get(suff('hBase_MET_Tau'))

                    plot2D(hBase, v, args.channel, sample, 'Base', cat, from_directory, str(ireg))
                    plot2D(hMET, v, args.channel, sample, 'MET', cat, from_directory, str(ireg))
                    plot2D(hTau, v, args.channel, sample, 'Tau', cat, from_directory, str(ireg))

                # count only for one of the variables, it does not matter which
                nbinsx, nbinsy = hBase.GetNbinsX(), hBase.GetNbinsY()
                cBase = round(hBase.Integral(0, nbinsx+1, 0, nbinsy+1), 2)
                cMET  = round(hMET.Integral(0, nbinsx+1, 0, nbinsy+1), 2)
                cTau  = round(hTau.Integral(0, nbinsx+1, 0, nbinsy+1), 2)
                
                cBase_MET = round(hBase_MET.Integral(0, nbinsx+1, 0, nbinsy+1), 2)
                cBase_Tau = round(hBase_Tau.Integral(0, nbinsx+1, 0, nbinsy+1), 2)
                cMET_Tau  = round(hMET_Tau.Integral(0, nbinsx+1, 0, nbinsy+1), 2)
                
                cNoBase_MET       = round(hNoBase_MET.Integral(0, nbinsx+1, 0, nbinsy+1), 2)
                cNoBase_NoMET_Tau = round(hNoBase_NoMET_Tau.Integral(0, nbinsx+1, 0, nbinsy+1), 2)
                cNoBase_Tau       = round(hNoBase_Tau.Integral(0, nbinsx+1, 0, nbinsy+1), 2)
                cNoBase_MET_NoTau = round(hNoBase_MET_NoTau.Integral(0, nbinsx+1, 0, nbinsy+1), 2)
                
                cBase_MET_Tau = round(hBase_MET_Tau.Integral(0, nbinsx+1, 0, nbinsy+1), 2)

                if ireg==0:
                    subsets = (cBase, cMET, cBase_MET, cTau, cBase_Tau, cMET_Tau, cBase_MET_Tau)
                    plt.figure(figsize=(10,8))
                    v = venn3(subsets=subsets, set_labels = (r'$\tau\tau$', 'MET', 'Tau'))
                    v.get_patch_by_id('100').set_alpha(1.0)
                    #v.get_patch_by_id('100').set_color('white')
                    #v.get_label_by_id('100').set_text('Unknown')
                    c = venn3_circles(subsets=subsets, linestyle='dashed')
                    #c[0].set_lw(1.0)
                    #c[0].set_ls('dotted')
                    #plt.title('Baseline, MET, Tau')
                    for ext in ('png', 'pdf'):
                        plt.savefig(os.path.join(out_counts[categories.index(cat)], 'venn.' + ext))

                with open(os.path.join(out_counts[categories.index(cat)], 'table.csv'), 'a') as f:
                    reader = csv.writer(f, delimiter=',', quotechar='|')
                    row = [ireg, cBase, cMET, cTau,
                           cBase_MET, cBase_Tau, cMET_Tau,
                           cNoBase_MET, cNoBase_NoMET_Tau, cNoBase_Tau, cNoBase_MET_NoTau,
                           cBase_MET_Tau]
                    reader.writerow(row)
                        
        f_in.Close()

    if args.copy:
        import subprocess
        to_directory = os.path.join('/eos/user/b/bfontana/www/TriggerScaleFactors', main_dir)
        to_directory = os.path.join(to_directory, args.channel)
     
        for sample in args.samples:
            sample_from = os.path.join(from_directory, sample)
            print('Copying: {}\t\t--->\t{}'.format(sample_from, to_directory), flush=True)
            subprocess.run(['rsync', '-ah', sample_from, to_directory])
     
    print('Done.')
