"""
Script which calculates the trigger scale factors.
On production mode should run in the grid via scripts/submitTriggerEff.py. 
"""
import re
import os
import functools
import argparse
import fnmatch
import math
from array import array
import numpy as np
from ROOT import (
    TFile,
    TH1D,
    TH2D,
)
import h5py
from collections import defaultdict
import itertools as it

import sys
sys.path.append( os.environ['PWD'] ) 

from utils.utils import (
    add_vnames,
    generate_trigger_combinations,
    get_hnames,
    get_trigger_bit,
    is_channel_consistent,
    join_name_trigger_intersection as joinNTC,
    LeafManager,
    load_binning,
    pass_any_trigger,
    pass_selection_cuts,
    rewrite_cut_string,
    parse_args,
    pass_trigger_bits,
    print_configuration,
)

from luigi_conf import (
    _2Dpairs,
    _cuts,
    _cuts_ignored,
)

def passes_cuts(trig, variables, leavesmanager, debug):
    """
    Handles cuts on trigger variables that enter the histograms. 
    Variables being displayed are not cut (i.e., they pass the cut).
    Checks all combinations of cuts specified in '_cuts':
        example: _cuts = {'A': ('>', [10,20]), 'B': ('<', [50,40]))}
        `passes_cuts` will check 4 combinations and return a dict of length 4
        (unless some cuts are ignored according to '_cuts_ignored') 
    Works for both 1D and 2D efficiencies.
    """
    if debug:
        print('Trigger={}; Variables={}'.format(trig, variables))

    flagnameJoin = lambda var,sign,val: ('_'.join([str(x) for x in [var,sign,val]])).replace('.','p')
    
    dflags = defaultdict(lambda: [])
    
    try:
        trig_cuts = _cuts[trig]
    except KeyError: # the trigger has no cut associated
        if debug:
            print('KeyError')            
        return {args.nocut_dummy_str: True}

    for avar,acut in trig_cuts.items():
        # ignore cuts according to the user's definition in '_cuts_ignored'
        # example: do not cut on 'met_et' when displaying 'metnomu_et'
        ignore = functools.reduce( lambda x, y: x or y,
                                   [ avar in _cuts_ignored[k] for k in variables
                                     if k in _cuts_ignored ],
                                   False
                                  )

        # additionally, by default do not cut on the variable(s) being plotted
        if avar not in variables and not ignore:
            value = leavesmanager.get_leaf(avar) 

            for c in acut[1]:
                flagname = flagnameJoin(avar, acut[0], c)

                if debug:
                    print('Cut: {} {} {}'.format(avar, acut[0], c))

                if acut[0]=='>':
                    dflags[avar].append( (flagname, value > c) )
                elif acut[0]=='<':
                    dflags[avar].append( (flagname, value < c) )
                else:
                    raise ValueError("The operator for the cut is currently not supported: Use '>' or '<'.")

    def applyCutsCombinations(dflags):
        tmp = {}
        allNames = sorted(dflags)
        combinations = it.product(*(dflags[name] for name in allNames))
        combinations = list(combinations)

        if debug:
            print('--------- [applyCutsCombinations] ------------------')
            print('Variables being cut: {}'.format(allNames))
            print('All cut combinations: {}'.format(combinations))
            print('----------------------------------------------------')
            
        for comb in combinations:
            joinFlag = functools.reduce( lambda x,y: x and y, [k[1] for k in comb] )
            tmp[ '_AND_'.join([k[0] for k in comb]) ] = joinFlag

        return tmp

    if dflags:
        res = applyCutsCombinations(dflags)
    else:
        res = {args.nocut_dummy_str: True}
    return res
    
def get_trigger_eff_sig(indir, outdir, sample, fileName,
                        channels, variables, triggers,
                        subtag, tprefix, isdata, binedges_fname):
    # -- Check if outdir exists, if not create it
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    if not os.path.exists( os.path.join(outdir, sample) ):
        os.makedirs( os.path.join(outdir, sample) )
    outdir = os.path.join(outdir, sample)

    if not os.path.exists(fileName):
        raise ValueError('[' + os.path.basename(__file__) + '] {} does not exist.'.format(fileName))

    f_in = TFile( fileName )
    t_in = f_in.Get('HTauTauTree')

    binedges, nbins = load_binning( afile=binedges_fname, key=subtag,
                                    variables=variables, channels=channels )

    triggercomb = {}
    for chn in channels:
        triggercomb[chn] = generate_trigger_combinations(chn, triggers)
        
    # Define 1D histograms
    #  hRef: pass the reference trigger
    #  hTrig: pass the reference trigger + trigger under study
    hRef, hTrig = ({} for _ in range(2))

    for i in channels:
        hRef[i], hTrig[i] = ({} for _ in range(2))
        for j in variables:
            binning1D = (nbins[j][i], binedges[j][i])
            hTrig[i][j]={}
            hRef[i][j] = TH1D( get_hnames('Ref1D')(i, j), '', *binning1D)
            for tcomb in triggercomb[i]:
                hTrig[i][j][joinNTC(tcomb)]={}

    # Define 2D histograms
    #  h2Ref: pass the reference trigger
    #  h2Trig: pass the reference trigger + trigger under study
    h2Ref, h2Trig = ({} for _ in range(2))
    
    for i in channels:
        h2Ref[i], h2Trig[i] = ({} for _ in range(2))
        for onetrig in triggers:
            if onetrig in _2Dpairs.keys():
                combtrigs = {x for x in triggercomb[i] if onetrig in x}
                for combtrig in combtrigs:
                    for j in _2Dpairs[onetrig]:
                        binning2D = ( nbins[j[0]][i], binedges[j[0]][i],
                                      nbins[j[1]][i], binedges[j[1]][i] )
                        vname = add_vnames(j[0],j[1])
                        if vname not in h2Ref[i]:
                            h2Ref[i][vname] = TH2D(get_hnames('Ref2D')(i, vname), '', *binning2D)
                        if vname not in h2Trig[i]:
                            h2Trig[i][vname] = {}
                        h2Trig[i][vname][joinNTC(combtrig)] = {}

    lf = LeafManager( fileName, t_in )
    
    for entry in range(0,t_in.GetEntries()):
        t_in.GetEntry(entry)

        if not pass_selection_cuts(lf):
            continue

        trig_bit = lf.get_leaf('pass_triggerbit')
        run = lf.get_leaf('RunNumber')
        if not pass_any_trigger(triggers, trig_bit, run, isdata=isdata):
            continue

        #mcweight   = lf.get_leaf( "MC_weight" )
        pureweight = lf.get_leaf( "PUReweight" )
        trigsf     = lf.get_leaf( "trigSF" )
        lumi       = lf.get_leaf( "lumi" )
        idandiso   = lf.get_leaf( "IdAndIsoSF_deep_pt")
        
        #if np.isnan(mcweight): mcweight=1
        if np.isnan(pureweight): pureweight=1
        if np.isnan(trigsf): trigsf=1
        if np.isnan(lumi): lumi=1
        if np.isnan(idandiso): idandiso=1

        evt_weight = pureweight*trigsf*lumi*idandiso
        if np.isnan(evt_weight) or isdata:
            evt_weight = 1

        fill_var = {}
        for v in variables:
            fill_var[v] = {}
            for chn in channels:
                fill_var[v].update({chn: lf.get_leaf(v)})
                if fill_var[v][chn]>binedges[v][chn][-1]:
                    fill_var[v][chn]=binedges[v][chn][-1] # include overflow

        # whether the event passes cuts (1 and 2-dimensional_ and single triggers
        pass_trigger, pcuts1D, pcuts2D = ({} for _ in range(3))
        for trig in triggers:
            pcuts2D[trig] = {}    
        for trig in triggers:
            pass_trigger[trig] = pass_trigger_bits(trig, trig_bit, run, isdata)

            pcuts1D[trig] = {}
            for var in variables:
                pcuts1D[trig][var] = passes_cuts(trig, [var], lf, args.debug)

            if trig in _2Dpairs.keys():
                                # combtrigs = tuple(x for x in triggercomb if trig in x)
                # for combtrig in combtrigs:
                # pcuts2D[joinNTC(combtrig)] = {}
                for j in _2Dpairs[trig]:
                    vname = add_vnames(j[0],j[1])
                    for t in triggers:
                        pcuts2D[t][vname] = passes_cuts(t, [j[0], j[1]], lf, args.debug)

        #logic AND to intersect all triggers in this combination
        pass_trigger_intersection = {}
        for i in channels:
            for tcomb in triggercomb[i]:
                pass_trigger_intersection[joinNTC(tcomb)] = functools.reduce(
                    lambda x,y: x and y,
                    [ pass_trigger[x] for x in tcomb ] )

        for i in channels:
            if is_channel_consistent(i, lf.get_leaf('pairType')):

                # fill histograms for 1D efficiencies
                for j in variables:
                    binning1D = (nbins[j][i], binedges[j][i])
                    hRef[i][j].Fill(fill_var[j][i], evt_weight)

                    # The following is tricky, as we are considering, simultaneously:
                    # - all trigger intersection combinations
                    # - all cut combinations for each trigger combination (see '_cuts')
                        
                    # Logic AND to intersect all cuts for this trigger combination
                    # Each element will contain one possible cut combination
                    # for the trigger combination 'tcomb' being considered
                    for tcomb in triggercomb[i]:
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

                        if args.debug:
                            print(j)
                            print(cuts_combinations)
                            print(pcuts_inters)
                            print()

                        for key,val in pcuts_inters.items():
                            if key not in hTrig[i][j][joinNTC(tcomb)]:
                                base_str = get_hnames('Trig1D')(i,j,joinNTC(tcomb))
                                htrig_name = rewrite_cut_string(base_str, key)
                                hTrig[i][j][joinNTC(tcomb)][key] = TH1D(htrig_name, '', *binning1D)

                            if val and pass_trigger_intersection[joinNTC(tcomb)]:
                                hTrig[i][j][joinNTC(tcomb)][key].Fill(fill_var[j][i], evt_weight)

                # fill 2D efficiencies
                for onetrig in triggers:
                    if onetrig in _2Dpairs.keys():
                        combtrigs = tuple(x for x in triggercomb[i] if onetrig in x)

                        for combtrig in combtrigs:
                            for j in _2Dpairs[onetrig]:
                                vname = add_vnames(j[0],j[1])
                                fill_info = ( fill_var[j[0]][i], fill_var[j[1]][i], evt_weight )
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

                                if combtrig==combtrigs[0]: #avoid filling multiple times
                                    h2Ref[i][vname].Fill(*fill_info)
    
                                for key,val in pcuts_inters.items():
                                    if key not in h2Trig[i][vname][joinNTC(combtrig)]:
                                        base_str = get_hnames('Trig2D')(i, vname, joinNTC(combtrig))
                                        h2name = rewrite_cut_string(base_str, key)
                                        binning2D = ( nbins[j[0]][i], binedges[j[0]][i],  
                                                      nbins[j[1]][i], binedges[j[1]][i] ) 
                                        h2Trig[i][vname][joinNTC(combtrig)][key] = TH2D(h2name, '',
                                                                                        *binning2D)
                       
                                    if val and pass_trigger_intersection[joinNTC(combtrig)]:
                                        h2Trig[i][vname][joinNTC(combtrig)][key].Fill(*fill_info)

    file_id = ''.join( c for c in fileName[-10:] if c.isdigit() ) 
    outname = os.path.join(outdir, tprefix + sample + '_' + file_id + subtag + '.root')
    print('Saving file {} at {} '.format(file_id, outname) )

    f_out = TFile(outname, 'RECREATE')
    f_out.cd()
    for i in channels:
        for j in variables:
            hRef[i][j].Write( get_hnames('Ref1D')(i,j) )
            for tcomb in triggercomb[i]:
                for khist,vhist in hTrig[i][j][joinNTC(tcomb)].items():
                    base_str = get_hnames('Trig1D')(i,j,joinNTC(tcomb))
                    writename = rewrite_cut_string(base_str, khist)
                    vhist.Write(writename)

        for vname in h2Ref[i].keys():
            h2Ref[i][vname].Write()
            for tc in h2Trig[i][vname].keys():
                for key in h2Trig[i][vname][tc].keys():
                    #print('\t|\t'.join([i, vname, tc, key]))
                    #h2name = get_hnames('Trig2D')(i, vname, joinNTC(tcomb))
                    base_str = h2Trig[i][vname][tc][key].GetName()
                    writename = rewrite_cut_string(base_str, key)
                    h2Trig[i][vname][tc][key].Write(writename)
    f_out.Close()
    f_in.Close()

# Parse input arguments
parser = argparse.ArgumentParser(description='Command line parser')

parser.add_argument('--binedges_fname', dest='binedges_fname', required=True, help='where the bin edges are stored')
parser.add_argument('--indir',       dest='indir',       required=True, help='SKIM directory')
parser.add_argument('--outdir',      dest='outdir',      required=True, help='output directory')
parser.add_argument('--sample',      dest='sample',      required=True, help='Process name as in SKIM directory')
parser.add_argument('--isdata',      dest='isdata',      required=True, help='Whether it is data or MC', type=int)
parser.add_argument('--file',        dest='fileName',    required=True, help='Full path of ROOT input file')
parser.add_argument('--subtag',      dest='subtag',      required=True,
                    help='Additional (sub)tag to differentiate similar runs within the same tag.')
parser.add_argument('--tprefix',     dest='tprefix',     required=True, help='Targets name prefix.')
parser.add_argument('--channels',    dest='channels',    required=True, nargs='+', type=str,  
                    help='Select the channels over which the workflow will be run.' )
parser.add_argument('--triggers',    dest='triggers',    required=True, nargs='+', type=str,
                    help='Select the triggers over which the workflow will be run.' )
parser.add_argument('--variables',   dest='variables',   required=True, nargs='+', type=str,
                    help='Select the variables over which the workflow will be run.' )
parser.add_argument('--intersection_str', dest='intersection_str', required=False, default='_PLUS_',
                    help='String used to represent set intersection between triggers.')
parser.add_argument('--nocut_dummy_str', dest='nocut_dummy_str', required=True,
                    help='Dummy string associated to trigger histograms were no cuts are applied.')
parser.add_argument('--debug', action='store_true', help='debug verbosity')
args = parse_args(parser)

get_trigger_eff_sig(args.indir, args.outdir, args.sample, args.fileName,
                    args.channels, args.variables, args.triggers,
                    args.subtag, args.tprefix, args.isdata, args.binedges_fname)
