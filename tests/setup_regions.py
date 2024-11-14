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
    
def get_outname(channel, bigtau):
    utils.create_single_dir('data')
    
    name = ""
    name += '.root'

    s = 'data/regions_{}'.format(name)
    return s


def trigger_regions(indir, channel, year, outname):
    outname = get_outname(outname)
    config_module = importlib.import_module(args.configuration)

    if channel == 'etau' or channel == 'mutau':
        iso1 = range(0, 8.1, 8/24)
    elif channel == 'tautau':
        iso1 = range(0, 401, 8)
    pNet_dist = [x/12. for x in range(13)]
    binning.update({
                    'genHH_mass': ([250, 300, 350, 400, 450, 500, 550, 600, 675, 800, 1000, 1600],),
                    'dau1_iso': (iso1,),
                    'dau1_pt': ([0, 20, 30, 40, 60, 80, 100, 125, 150, 200, 250],),
                    'dau2_iso': (iso1,),
                    'dau2_pt': ([0, 20, 30, 40, 60, 80, 100, 125, 150, 200, 250],),
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
                    # 'triggerbits': (range(31))
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
                z = entries[i]

            ahistos[chn][i].fill(z, weight=evt_weight)
        
        
    # all MC and signal must be rescaled to get the correct number of events
    # for key,_ in cuts.items():
    #     for reg in regions:
    #         for cat in categories:
    #             for cc in ['mutau', 'etau', 'tautau']:
    #                 print("lumi", utils.get_lumi(args.year), utils.total_sum_weights(glob_files[0].replace("PreprocessRDF", "PreCounter").replace("/cat_base_selection", "").replace(".root", ".json"), isdata=False))
    #                 ahistos[key][reg][cat][cc] *= (utils.get_lumi(args.year) /
    #                                     utils.total_sum_weights(glob_files[0].replace("PreprocessRDF", "PreCounter").replace("/cat_base_selection", "").replace(".root", ".json"), isdata=False))
    
    file = uproot.recreate(outname)
    print('Raw histograms saved in {}.'.format(outname), flush=True)

if __name__ == '__main__':
    extensions = ('png',) #('png', 'pdf')
    triggers = {'etau': ("HLT_Ele30_WPTight_Gsf",
                    "HLT_Ele24_eta2p1_WPTight_Gsf_LooseDeepTauPFTauHPS30_eta2p1_CrossL1"),
                'mutau': ('HLT_IsoMu24', 'HLT_IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1'),
                'tautau': ('HLT_DoubleMediumDeepTauPFTauHPS35_L2NN_eta2p1',)
                }

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
    parser.add_argument('--configuration', dest='configuration', required=True,
                        help='Name of the configuration module to use.')
    parser.add_argument('--outname', dest='outname', required=True,
                        help='Output file name.')
    args = utils.parse_args(parser)
        
    #### run main function ###
    if not args.plot:
        # for sample in args.masses:
        trigger_regions(args.indir, args.channel, args.year, args.outname)

    ###########################

