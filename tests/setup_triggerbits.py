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
import numpy as np
import h5py
from collections import defaultdict as dd
import importlib
from array import array

import inclusion
from inclusion import selection
from inclusion.config import main
from inclusion.utils import utils

import ROOT
import hist
import pickle
import uproot

from bokeh.plotting import figure, output_file, save
from bokeh.models import Range1d, Label

tau = '\u03C4'
mu  = '\u03BC'
pm  = '\u00B1'
ditau = tau+tau
        
def rec_dd():
    return dd(rec_dd)
    
def get_outname(channel, bigtau):
    utils.create_single_dir('data')
    
    name = ""
    if bigtau:
        name += '_BIGTAU'
    name += '_all.root'

    s = 'data/regions_2023_postBPix_15p36-ditaujet{}'.format(name)
    return s


def trigger_regions(indir, channel, year, deltaR):
    outname = get_outname(channel,
                          args.bigtau)
    config_module = importlib.import_module(args.configuration)

    if channel == 'etau' or channel == 'mutau':
        iso1 = range(0, 8.1, 8/24)
    elif channel == 'tautau':
        iso1 = range(0, 401, 8)
    pNet_dist = [x/12. for x in range(13)]
    binning.update({
                    'genHH_mass': ([250, 300, 350, 400, 450, 500, 550, 600, 675, 800, 1000, 1600],),
                    'dau1_iso': (iso1,),
                    # 'dau1_pt': ([0, 20, 30, 40, 60, 80, 100, 125, 150, 200, 250],),
                    'dau2_iso': (iso1,),
                    # 'dau2_pt': ([0, 20, 30, 40, 60, 80, 100, 125, 150, 200, 250],),
                    'dau1_eta': ([-3, -2.5, -2.1, -1.8, -1.5, -1.2, -0.8, -0.4, 0, 0.4, 0.8, 1.2, 1.5, 1.8, 2.1, 2.5, 3],),
                    'dau2_eta': ([-3, -2.5, -2.1, -1.8, -1.5, -1.2, -0.8, -0.4, 0, 0.4, 0.8, 1.2, 1.5, 1.8, 2.1, 2.5, 3],),
                    'dau1_tauIdVSjet': ([0, 1, 2, 3, 4, 5, 6, 7], ),
                    'dau2_tauIdVSjet': ([0, 1, 2, 3, 4, 5, 6, 7], ),
                    'bjet1_pNet': (pNet_dist,),
                    'bjet2_pNet': (pNet_dist,),
                    'bjet1_pt': ([0, 20, 30, 40, 60, 80, 100, 125, 150, 200, 250],),
                    'bjet2_pt': ([0, 20, 30, 40, 60, 80, 100, 125, 150, 200, 250],),
                    'bjet1_eta': ([-3, -2.5, -2.1, -1.8, -1.5, -1.2, -0.8, -0.4, 0, 0.4, 0.8, 1.2, 1.5, 1.8, 2.1, 2.5, 3],),
                    'bjet2_eta': ([-3, -2.5, -2.1, -1.8, -1.5, -1.2, -0.8, -0.4, 0, 0.4, 0.8, 1.2, 1.5, 1.8, 2.1, 2.5, 3],),
                    # 'bjet1_triggerbits': ([i for i in range(32)],),
                    # 'bjet2_triggerbits': ([i for i in range(32)],),
                    # 'tau1_triggerbits': ([i for i in range(32)],),
                    # 'tau2_triggerbits': ([i for i in range(32)],),
                    # 'SoftActivityJetHT': ([0, 50, 100, 150, 200, 250, 300, 350, 400, 500, 750, 1000, 1500, 2000],)
                    })

    
    norphans, ntotal = ({k:0 for k in categories} for _ in range(2))
    
    ahistos = rec_dd()        
    for chn in ['mutau', 'etau', 'tautau']:
        for i in binning.keys():
            ahistos[chn][i] = (
                hist.Hist.new.Variable(*binning[i], name=i)
                .Weight()
            )

    t_in = ROOT.TChain('Events')
    glob_files = glob.glob( os.path.join(indir, 'data_*.root') )
    if len(glob_files) < 1:
        raise RuntimeError("No files!")
    for f in glob_files:
        t_in.Add(f)
    t_in.SetBranchStatus('*', 0)
    _entries, mc_corrections = utils.define_used_tree_variables(cut=config_module.custom_cut)

    _entries += tuple([ x + mc_corrections[x] for x in mc_corrections])
    for ientry in _entries:
        t_in.SetBranchStatus(ientry, 1)

    for entry in t_in:
        # this is slow: do it once only
        entries = utils.dot_dict({x: getattr(entry, x) for x in _entries})

        if entries.pairType == -1:
            continue

        if (entries.pairType == 0 or entries.pairType == 1) and (entries.isOS != 1 or entries.dau2_tauIdVSjet < 5):
            continue

        if entries.pairType == 2 and (entries.isOS != 1 or entries.dau2_tauIdVSjet < 5 or entries.dau1_tauIdVSjet < 5):
            continue
            
        if entries.pairType > 2:
            continue

     
        sel = selection.EventSelection(entries, year=year, isdata=False, configuration=config_module)
        # in_legacy, in_met, in_tau = which_region(entries, year, ptcuts, regcuts, channel, sel,
        #                                          bigtau=args.bigtau)
        
        w_mc     = entries.genWeight
        # w_pure   = entries.puWweight
        # # w_l1pref = entries.L1pref_weight
        # w_trig   = 1 # entries.trigSF
        # w_idiso  = entries.IdSF_deep_2d
        # w_jetpu  = entries.PUjetID_SF
        # w_btag   = entries.bTagweightReshape
        
        if utils.is_nan(w_mc)     : w_mc=1
        # if utils.is_nan(w_pure)   : w_pure=1
        # # if utils.is_nan(w_l1pref) : w_l1pref=1
        # if utils.is_nan(w_trig)   : w_trig=1
        # if utils.is_nan(w_idiso)  : w_idiso=1
        # if utils.is_nan(w_jetpu)  : w_jetpu=1
        # if utils.is_nan(w_btag)   : w_btag=1
        evt_weight = 1 # * w_btag
        # if evt_weight < 0.:
        #     print(w_mc, w_pure, w_trig, w_idiso, w_jetpu)
        tau_gen_cut = {"etau": None, "mutau": None,
                       "tautau": 'self.entries["isTau1real"] == 1 and self.entries["isTau2real"] == 1'}
        
        if entries.pairType == 0:
            chn = 'mutau'
        elif entries.pairType == 1:
            chn = 'etau'
        elif entries.pairType == 2:
            chn = 'tautau'

       

        for i in binning.keys():
            z = 0
            if 'bjet1' in i:
                if 'eta' in i:
                    z = entries['Jet_eta'][entries['bjet1_JetIdx']]
                if 'pt' in i:
                    z = entries['Jet_pt'][entries['bjet1_JetIdx']]
                if 'pNet' in i:
                    z = entries['Jet_btagPNetB'][entries['bjet1_JetIdx']]
            elif 'bjet2' in i:
                if 'eta' in i:
                    z = entries['Jet_eta'][entries['bjet2_JetIdx']]
                if 'pt' in i:
                    z = entries['Jet_pt'][entries['bjet2_JetIdx']]
                if 'pNet' in i:
                    z = entries['Jet_btagPNetB'][entries['bjet2_JetIdx']]
            else:
                if 'triggerbits' in i:
                    continue
                z = entries[i]

            ahistos[chn][i].fill(z, weight=evt_weight)
        
        # if chn == 'tautau' and entries['isQuadJetTrigger'] and entries['bjet1_filterbits'] >= 0:
        #     for i in range(31):
        #         print(i, entries["bjet1_filterbits"], bool(entries["bjet1_filterbits"] & (1 << i)))
        #         if bool(entries["bjet1_filterbits"] & (1 << i)):
        #             ahistos[chn]['bjet1_triggerbits'].fill(i)
        #         if bool(entries["bjet2_filterbits"] & (1 << i)):
        #             ahistos[chn]['bjet2_triggerbits'].fill(i)
        #         if bool(entries["tau1_filterbits"] & (1 << i)):
        #             ahistos[chn]['tau1_triggerbits'].fill(i)
        #         if bool(entries["tau2_filterbits"] & (1 << i)):
        #             ahistos[chn]['tau2_triggerbits'].fill(i)

   
    file = uproot.recreate(outname)
    for chn in ['mutau', 'etau', 'tautau']:
        for i in binning.keys():
            file[str(chn) + "_" + str(i)] = ahistos[chn][i]
   
    print('Raw histograms saved in {}.'.format(outname), flush=True)

if __name__ == '__main__':
    extensions = ('png',) #('png', 'pdf')
    triggers = {'etau': ("HLT_Ele30_WPTight_Gsf",
                    "HLT_Ele24_eta2p1_WPTight_Gsf_LooseDeepTauPFTauHPS30_eta2p1_CrossL1"),
                'mutau': ('HLT_IsoMu24', 'HLT_IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1'),
                'tautau': ('HLT_DoubleMediumDeepTauPFTauHPS35_L2NN_eta2p1',)
                }
    binning = {
        # 'metnomu_et': (20, 0, 450),
        #        'dau1_pt': (30, 0, 450),
        #        'dau1_eta': (20, -2.5, 2.5),
        #        'dau2_iso': (20, 0.88, 1.005),
        #        'dau2_pt': (30, 0, 400),
        #        'dau2_eta': (20, -2.5, 2.5),
        #        'ditau_deltaR': (30, 0.3, 1.3),
        #        'dib_deltaR': (25, 0, 2.5),
        #        'bH_pt': (20, 70, 600),
        #        'bH_mass': (30, 0, 280),
        #        'tauH_mass': (30, 0, 170),
        #        'tauH_pt': (30, 0, 500),
        #        'tauH_SVFIT_mass': (30, 0, 250),
        #        'tauH_SVFIT_pt': (20, 200, 650),
        #        'bjet1_pt': (25, 10, 600),
        #        'bjet2_pt': (25, 10, 550),
        #        'bjet1_eta': (20, -2.5, 2.5),
        #        'bjet2_eta': (20, -2.5, 2.5),
               }
    variables = tuple(binning.keys()) + ('HHKin_mass', 'dau1_iso')

    categories = ('baseline',) #('baseline', 's1b1jresolvedMcut', 's2b0jresolvedMcut', 'sboostedLLMcut')
    
    # Parse input arguments
    desc = 'Producer trigger histograms.\n'
    desc += "Run example: python tests/test_trigger_regions.py --indir /data_CMS/cms/alves/HHresonant_SKIMS/SKIMS_UL18_EOSv4_Signal/ --masses 400 500 600 700 800 900 1000 1250 1500 --channels ETau --met_turnon 180 --region_cuts 40 40 --copy"
    parser = argparse.ArgumentParser(description=desc, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('--indir', required=True, type=str,
                        help='Full path of ROOT input file')
    # parser.add_argument('--masses', required=True, nargs='+', type=str,
    #                     help='Resonance mass')
    parser.add_argument('--channel', required=True, type=str,  
                        help='Select the channel over which the workflow will be run.' )
    parser.add_argument('--year', required=True, type=str, choices=('2016', '2016APV', '2017', '2018', '2022'),
                        help='Select the year over which the workflow will be run.' )
    # parser.add_argument('--spin', required=True, type=int, choices=(0, 2), 
    #                     help='Select the spin hypothesis over which the workflow will be run.' )
    parser.add_argument('--deltaR', type=float, default=0.5, help='DeltaR between the two leptons.' )
    parser.add_argument('--plot', action='store_true',
                        help='Reuse previously produced data for quick plot changes.')
    parser.add_argument('--copy', action='store_true',
                        help='Copy the outputs to EOS at the end.')
    parser.add_argument('--notext', action='store_true', help='Square diagram without text.')
    parser.add_argument('--sequential', action='store_true',
                        help='Do not use the multiprocess package.')
    parser.add_argument('--bigtau', action='store_true',
                        help='Consider a larger single tau region, reducing the ditau one.')
    parser.add_argument('--met_turnon', type=float, default=180,
                        help='MET trigger turnon cut [GeV].' )
    parser.add_argument('--region_cuts', required=False, type=float, nargs=2, default=(190, 190),
                        help='High/low regions pT1 and pT2 selection cuts [GeV].' )
    parser.add_argument('--configuration', dest='configuration', required=True,
                        help='Name of the configuration module to use.')
    args = utils.parse_args(parser)

    met_turnon = args.met_turnon
    regcuts = args.region_cuts
    ptcuts = utils.get_ptcuts(args.channel, args.year)
        
    main_dir = os.path.join(os.path.join('/t3home/', os.environ['USER'], 'TriggerScaleFactors'),
                            '_'.join((args.channel, *[str(x) for x in regcuts],
                                      'DR', str(args.deltaR), 'PT', *[str(x) for x in ptcuts], 'TURNON',
                                      str(met_turnon))))
    if args.bigtau:
        main_dir += '_BIGTAU'

    regions = ('legacy', 'met', 'tau')
        
    #### run main function ###
    if not args.plot:
        # for sample in args.masses:
        trigger_regions(args.indir, args.channel, args.year, args.deltaR)
        # if args.sequential:
        #     pass
        # else:
        #     pool = multiprocessing.Pool(processes=6)
        #     pool.starmap(trigger_regions,
        #                  zip(it.repeat(args.indir), it.repeat(args.channel), it.repeat(args.year), it.repeat(args.deltaR)))

    ###########################

    sum_stats, err_sum_stats = ([] for _ in range(2))
    contam1, contam2, contam1_errors, contam2_errors = ([] for _ in range(4))
    from_directory = os.path.join(main_dir, args.channel)
    outname = get_outname(args.channel, args.bigtau)
    with open(outname, "rb") as f:
        ahistos = pickle.load(f)

    # write csv header, one per category
    out_counts = []
    # plot histograms and fill CSV with histogram integrals
    c_legacy_trg, c_met_trg, c_tau_trg = ({} for _ in range(3))

    acounts = rec_dd()
    for reg in regions:
        for key,cut in ahistos.items():
            acounts[key][reg] = round(ahistos[key][reg]["baseline"].values().sum(), 2)
        
        # append to table, one line per region
        with open(os.path.join(out_counts[categories.index(cat)], 'table.csv'), 'a') as f:
            reader = csv.writer(f, delimiter=',', quotechar='|')
            row = [reg]
            row.extend([acounts[k][reg] for k in acounts.keys()])
            reader.writerow(row)

        if reg=='legacy':
            c_legacy_trg[reg] = acounts["Base"][reg]
            c_met_trg[reg]    = acounts["NoBaseMET"][reg]
            c_tau_trg[reg]    = acounts["NoBaseNoMETTau"][reg]
        elif reg=='met':
            c_legacy_trg[reg] = acounts["BaseNoMET"][reg]
            c_met_trg[reg]    = acounts["MET"][reg]
            c_tau_trg[reg]    = acounts["NoBaseNoMETTau"][reg]
        elif reg=='tau':
            c_legacy_trg[reg] = acounts["BaseNoTau"][reg]
            c_met_trg[reg]    = acounts["NoBaseMETNoTau"][reg]
            c_tau_trg[reg]    = acounts["Tau"][reg]

    c1, c2, e1, e2 = sq_res
    contam1.append(c1)
    contam2.append(c2)
    contam1_errors.append(e1)
    contam2_errors.append(e2)

    stats_l = [c_legacy_trg['legacy'],c_met_trg['legacy'],c_tau_trg['tau'],c_legacy_trg['tau']]
    stats = sum(stats_l)
    estats = sum(np.sqrt(stats_l))
    sum_stats.append(stats)
    err_sum_stats.append(estats)

    contam1 = [float(x) for x in contam1]
    contam2 = [float(x) for x in contam2]
    contam1_errors  = [float(x) for x in contam1_errors]
    contam2_errors  = [float(x) for x in contam2_errors]
    # masses = [float(x) for x in args.masses]
    contamination_save('data', '_'.join([str(x) for x in regcuts]) + '_' + args.channel,
                       contam1, contam2, contam1_errors, contam2_errors)
    stats_save('data', '_'.join([str(x) for x in regcuts]) + '_' + args.channel,
               sum_stats, err_sum_stats, mode='a')

    if args.copy:
        import subprocess
        to_directory = os.path.join('/eos/home-b/bfontana/www/TriggerScaleFactors', main_dir)
        to_directory = os.path.join(to_directory, args.channel)
     
        for sample in args.masses:
            sample_from = os.path.join(from_directory, sample)
            print('Copying: {}\t\t--->\t{}'.format(sample_from, to_directory), flush=True)
            subprocess.run(['rsync', '-ah', sample_from, to_directory])
