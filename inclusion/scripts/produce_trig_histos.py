# coding: utf-8

_all_ = [ 'build_histograms' ]

import re
import os
import functools
import argparse
import h5py
import itertools as it

from ROOT import (
    TFile,
    TH1D,
    TH2D,
)

import os
import sys
parent_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, parent_dir)

import inclusion
from inclusion import selection
from inclusion.utils import utils
from inclusion.utils.utils import join_name_trigger_intersection as joinNTC

from inclusion import config

def build_histograms(infile, outdir, dataset, sample, isdata,
                     channels, variables, triggers,
                     subtag, tprefix, binedges_fname):
    # -- Check if outdir exists, if not create it
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    if not os.path.exists( os.path.join(outdir, sample) ):
        os.makedirs( os.path.join(outdir, sample) )
    outdir = os.path.join(outdir, sample)

    if not os.path.exists(infile):
        mes = '[' + os.path.basename(__file__) + '] {} does not exist.'.format(infile)
        raise ValueError(mes)

    f_in = TFile(infile)
    t_in = f_in.Get('HTauTauTree')

    binedges, nbins = utils.load_binning(afile=binedges_fname, key=subtag,
                                         variables=variables, channels=channels)

    triggercomb = {}
    for chn in channels:
        triggercomb[chn] = utils.generate_trigger_combinations(chn, triggers)

    # Define 1D histograms
    #  hRef: pass the reference trigger
    #  hTrig: pass the reference trigger + trigger under study
    hRef, hTrig = ({} for _ in range(2))

    for chn in channels:
        hRef[chn], hTrig[chn] = ({} for _ in range(2))
        for j in variables:
            binning1D = (nbins[j][chn], binedges[j][chn])
            hTrig[chn][j], hRef[chn][j] = ({} for _ in range(2))
            for tcomb in triggercomb[chn]:
                hname = utils.get_hnames('Ref1D')(chn, j, tcomb)
                hRef[chn][j][joinNTC(tcomb)] = TH1D(hname, '', *binning1D)
                hTrig[chn][j][joinNTC(tcomb)] = {}

    # Define 2D histograms
    #  h2Ref: pass the reference trigger
    #  h2Trig: pass the reference trigger + trigger under study
    h2Ref, h2Trig = ({} for _ in range(2))
    
    for chn in channels:
        h2Ref[chn], h2Trig[chn] = ({} for _ in range(2))
        for onetrig in triggers:
            if onetrig in config.pairs2D.keys():
                combtrigs = {x for x in triggercomb[chn] if onetrig in x}
                for combtrig in combtrigs:
                    cstr = joinNTC(combtrig)
                    
                    for j in config.pairs2D[onetrig]:
                        bin2D = ( nbins[j[0]][chn], binedges[j[0]][chn],
                                  nbins[j[1]][chn], binedges[j[1]][chn] )
                        vname = utils.add_vnames(j[0], j[1])
                        if vname not in h2Ref[chn]:
                            h2Ref[chn][vname] = {}
                        hname = utils.get_hnames('Ref2D')(chn, vname, cstr)
                        h2Ref[chn][vname][cstr] = TH2D(hname, '', *bin2D)
                        if vname not in h2Trig[chn]:
                            h2Trig[chn][vname] = {}
                        h2Trig[chn][vname][cstr] = {}

    t_in.SetBranchStatus('*', 0)
    _entries = ('pass_triggerbit', 'RunNumber',
                'PUReweight', 'lumi', 'IdAndIsoSF_deep_pt',
                'HHKin_mass', 'pairType', 'dau1_eleMVAiso', 'dau1_iso', 'dau1_deepTauVsJet', 'dau2_deepTauVsJet',
                'nleps', 'nbjetscand', 'tauH_SVFIT_mass', 'bH_mass_raw',)
    _entries += tuple(variables)
    for ientry in _entries:
        t_in.SetBranchStatus(ientry, 1)

    nentries = t_in.GetEntriesFast()
    for ientry,entry in enumerate(t_in):
        if ientry%10000==0:
             print('Processed {} entries out of {}.'.format(ientry, nentries))

        # this is slow: do it once only
        entries = utils.dot_dict({x: getattr(entry, x) for x in _entries})
        
        sel = selection.EventSelection(entries, dataset, isdata)

        #mcweight   = entries.MC_weight
        pureweight = entries.PUReweight
        lumi       = entries.lumi
        idandiso   = entries.IdAndIsoSF_deep_pt
        
        #if utils.is_nan(mcweight): mcweight=1
        if utils.is_nan(pureweight): pureweight=1
        if utils.is_nan(lumi): lumi=1
        if utils.is_nan(idandiso): idandiso=1

        evt_weight = pureweight*lumi*idandiso
        if utils.is_nan(evt_weight) or isdata:
            evt_weight = 1

        fill_var = {}
        for v in variables:
            fill_var[v] = {}
            for chn in channels:
                fill_var[v].update({chn: entries[v]})
                if fill_var[v][chn]>binedges[v][chn][-1]:
                    fill_var[v][chn]=binedges[v][chn][-1] # include overflow

        # whether the event passes cuts (1 and 2-dimensional_ and single triggers
        pass_trigger, pcuts1D, pcuts2D = ({} for _ in range(3))
        for trig in triggers:
            pcuts2D[trig] = {}
        for trig in triggers:
            pass_trigger[trig] = sel.trigger_bits(trig)

            pcuts1D[trig] = {}
            for var in variables:
                pcuts1D[trig][var] = sel.var_cuts(trig, [var], args.nocut_dummy_str)

            if trig in config.pairs2D.keys():
                # combtrigs = tuple(x for x in triggercomb if trig in x)
                # for combtrig in combtrigs:
                # pcuts2D[joinNTC(combtrig)] = {}
                for j in config.pairs2D[trig]:
                    vname = utils.add_vnames(j[0],j[1])
                    for t in triggers:
                        pcuts2D[t][vname] = sel.var_cuts(t, [j[0], j[1]], args.nocut_dummy_str)

        #logic AND to intersect all triggers in this combination
        pass_trigger_intersection = {}
        for chn in channels:
            for tcomb in triggercomb[chn]:
                pass_trigger_intersection[joinNTC(tcomb)] = functools.reduce(
                    lambda x,y: x and y,
                    [ pass_trigger[x] for x in tcomb ] )

        for chn in channels:
            if utils.is_channel_consistent(chn, entries.pairType):
                
                # fill histograms for 1D efficiencies
                for j in variables:
                    binning1D = (nbins[j][chn], binedges[j][chn])

                    # The following is tricky, as we are considering, simultaneously:
                    # - all trigger intersection combinations
                    # - all cut combinations for each trigger combination (see '_cuts')
                        
                    # Logic AND to intersect all cuts for this trigger combination
                    # Each element will contain one possible cut combination
                    # for the trigger combination 'tcomb' being considered
                    for tcomb in triggercomb[chn]:
                        cstr = joinNTC(tcomb)

                        if not sel.dataset_cuts(tcomb, chn):
                            continue
                        if not sel.dataset_triggers(triggers, tcomb, chn)[0]:
                            continue
                        if not sel.match_inters_with_dataset(tcomb, chn):
                            continue

                        hRef[chn][j][cstr].Fill(fill_var[j][chn], evt_weight)
                        
                        cuts_combinations = list(it.product( *(pcuts1D[atrig][j].items()
                                                             for atrig in tcomb) ))

                        # One dict item per cut combination
                        # - key: all cut strings joined
                        # - value: logical and of all cuts
                        pcuts_inters = { (args.intersection_str).join(e[0] for e in elem): 
                                        functools.reduce(                 
                                            lambda x,y: x and y,
                                            [ e[1] for e in elem ]
                                            )
                                        for elem in cuts_combinations }

                        for key,val in pcuts_inters.items():
                            if key not in hTrig[chn][j][cstr]:
                                base_str = utils.get_hnames('Trig1D')(chn,j,cstr)
                                htrig_name = utils.rewrite_cut_string(base_str, key)
                                hTrig[chn][j][cstr][key] = TH1D(htrig_name, '', *binning1D)

                            if val and pass_trigger_intersection[cstr]:
                                hTrig[chn][j][cstr][key].Fill(fill_var[j][chn], evt_weight)

                # fill 2D efficiencies
                for onetrig in triggers:
                    if onetrig in config.pairs2D.keys():
                        combtrigs = tuple(x for x in triggercomb[chn] if onetrig in x)

                        for combtrig in combtrigs:
                            cstr = joinNTC(combtrig)
                            
                            if not sel.dataset_cuts(combtrig, chn):
                                continue
                            if not sel.dataset_triggers(triggers, combtrig, chn)[0]:
                                continue
                            if not sel.match_inters_with_dataset(combtrig, chn):
                                continue
                            
                            for j in config.pairs2D[onetrig]:
                                vname = utils.add_vnames(j[0],j[1])
                                fill_info = ( fill_var[j[0]][chn], fill_var[j[1]][chn],
                                              evt_weight )
                                try:
                                    cuts_combinations = list(it.product(
                                        *(pcuts2D[atrig][vname].items() for atrig in combtrig) ))
                                except KeyError:
                                    print('CHECK: ', combtrig)
                                    raise

                                if args.debug:
                                    print(combtrig, j)
                                    print(cuts_combinations)

                                # One dict item per cut combination
                                # - key: all cut strings joined
                                # - value: logical and of all cuts
                                pcuts_inters = { (args.intersection_str).join(e[0] for e in elem):
                                                 functools.reduce(
                                                     lambda x,y: x and y,
                                                     [ e[1] for e in elem ]
                                                     )
                                                 for elem in cuts_combinations }

                                h2Ref[chn][vname][cstr].Fill(*fill_info)
    
                                for key,val in pcuts_inters.items():
                                    if key not in h2Trig[chn][vname][cstr]:
                                        base_str = utils.get_hnames('Trig2D')(chn, vname, cstr)
                                        h2name = utils.rewrite_cut_string(base_str, key)
                                        bin2D = ( nbins[j[0]][chn], binedges[j[0]][chn],  
                                                  nbins[j[1]][chn], binedges[j[1]][chn] ) 
                                        h2Trig[chn][vname][cstr][key] = TH2D(h2name, '', *bin2D)
                                        
                                    if val and pass_trigger_intersection[cstr]:
                                        h2Trig[chn][vname][cstr][key].Fill(*fill_info)
    
    file_id = ''.join( c for c in infile[-10:] if c.isdigit() ) 
    outname = os.path.join(outdir, tprefix + sample + '_' + file_id + subtag + '.root')
    empty_files = True
    
    f_out = TFile(outname, 'RECREATE')
    f_out.cd()
    for chn in channels:
        for j in variables:
            for tcomb in triggercomb[chn]:
                cstr = joinNTC(tcomb)
                if hRef[chn][j][cstr].GetEntries() > 0:
                    empty_files = False
                
                hRef[chn][j][cstr].Write( utils.get_hnames('Ref1D')(chn,j,cstr) )

                for khist,vhist in hTrig[chn][j][cstr].items():
                    base_str = utils.get_hnames('Trig1D')(chn, j, cstr)
                    writename = utils.rewrite_cut_string(base_str, khist)
                    vhist.Write(writename)

        for vname in h2Ref[chn].keys():

            for tc in h2Trig[chn][vname].keys():
                h2Ref[chn][vname][tc].Write()
                
                for key in h2Trig[chn][vname][tc].keys():
                    base_str = h2Trig[chn][vname][tc][key].GetName()
                    writename = utils.rewrite_cut_string(base_str, key)
                    h2Trig[chn][vname][tc][key].Write(writename)

    if empty_files:
        mes = 'All 1D histograms are empty.'
        print('WARNING: ' + mes)

    f_in.Close()
    f_out.Close()
    print('Saving file {} at {} '.format(file_id, outname) )

# Parse input arguments
parser = argparse.ArgumentParser(description='Command line parser')

parser.add_argument('--binedges_fname', dest='binedges_fname', required=True, help='where the bin edges are stored')
parser.add_argument('--outdir', dest='outdir', required=True,
                    help='output directory')
parser.add_argument('--dataset', dest='dataset', required=True,
                    help='Dataset name as provided by the user: MET, EG, ...')
parser.add_argument('--sample', dest='sample', required=True,
                    help='Process name as in SKIM directory')
parser.add_argument('--isdata', dest='isdata', required=True, type=int,
                    help='Whether it is data or MC')
parser.add_argument('--file', dest='infile', required=True,
                    help='Full path of ROOT input file')
parser.add_argument('--subtag', dest='subtag', required=True,
                    help='Additional (sub)tag to differentiate similar runs within the same tag.')
parser.add_argument('--tprefix', dest='tprefix', required=True, help='Targets name prefix.')
parser.add_argument('--channels', dest='channels', required=True, nargs='+', type=str,  
                    help='Select the channels over which the workflow will be run.' )
parser.add_argument('--triggers', dest='triggers', required=True, nargs='+', type=str,
                    help='Select the triggers over which the workflow will be run.' )
parser.add_argument('--variables',   dest='variables', required=True, nargs='+', type=str,
                    help='Select the variables over which the workflow will be run.' )
parser.add_argument('--intersection_str', dest='intersection_str', required=False, default='_PLUS_',
                    help='String used to represent set intersection between triggers.')
parser.add_argument('--nocut_dummy_str', dest='nocut_dummy_str', required=True,
                    help='Dummy string associated to trigger histograms were no cuts are applied.')
parser.add_argument('--debug', action='store_true', help='debug verbosity')
args = utils.parse_args(parser)

build_histograms(args.infile, args.outdir, args.dataset, args.sample, args.isdata,
                 args.channels, args.variables, args.triggers,
                 args.subtag, args.tprefix, args.binedges_fname)
