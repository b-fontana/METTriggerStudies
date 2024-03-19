# coding: utf-8

_all_ = [ 'build_histograms' ]

import os
import sys
parent_dir = os.path.abspath(__file__ + 3 * '/..')
sys.path.insert(0, parent_dir)

import inclusion
from inclusion import selection
from inclusion.config import main
from inclusion.utils import utils
from inclusion.utils.utils import join_name_trigger_intersection as joinNTC

import re
import functools
import argparse
import itertools as it
import importlib

import ROOT

def build_histograms(args):
    # -- Check if outdir exists, if not create it
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
    if not os.path.exists( os.path.join(args.outdir, args.sample) ):
        os.makedirs( os.path.join(args.outdir, args.sample) )
    outdir = os.path.join(args.outdir, args.sample)

    if not os.path.exists(args.infile):
        mes = '[' + os.path.basename(__file__) + '] {} does not exist.'.format(args.infile)
        raise ValueError(mes)

    config_module = importlib.import_module(args.configuration)

    f_in = ROOT.TFile.Open(args.infile)
    t_in = f_in.Get('HTauTauTree')

    binedges, nbins = utils.load_binning(afile=args.binedges_fname, key=args.subtag,
                                         variables=args.variables, channels=args.channels)
    triggercomb = {}
    for chn in args.channels:
        triggercomb[chn] = utils.generate_trigger_combinations(chn, config_module.triggers,
                                                               config_module.exclusive)

    # Define 1D histograms
    #  hRef: pass the reference trigger
    #  hTrg: pass the reference trigger + trigger under study
    hRef, hTrg = ({} for _ in range(2))

    for chn in args.channels:
        hRef[chn], hTrg[chn] = ({} for _ in range(2))
        for j in args.variables:
            binning1D = (nbins[j][chn], binedges[j][chn])
            hTrg[chn][j], hRef[chn][j] = ({} for _ in range(2))
            for tcomb in triggercomb[chn]:
                hname = utils.get_hnames('Ref1D')(chn, j, joinNTC(tcomb))

                hRef[chn][j][joinNTC(tcomb)] = ROOT.TH1D(hname, '', *binning1D)
                hTrg[chn][j][joinNTC(tcomb)] = {}

    # Define 2D histograms
    #  h2Ref: pass the reference trigger
    #  h2Trig: pass the reference trigger + trigger under study
    h2Ref, h2Trig = ({} for _ in range(2))
    
    for chn in args.channels:
        h2Ref[chn], h2Trig[chn] = ({} for _ in range(2))
        for onetrig in config_module.triggers:
            if onetrig in config_module.pairs2D.keys():
                combtrigs = {x for x in triggercomb[chn] if onetrig in x}
                for combtrig in combtrigs:
                    cstr = joinNTC(combtrig)
                    
                    for j in config_module.pairs2D[onetrig]:
                        bin2D = ( nbins[j[0]][chn], binedges[j[0]][chn],
                                  nbins[j[1]][chn], binedges[j[1]][chn] )
                        vname = utils.add_vnames(j[0], j[1])
                        if vname not in h2Ref[chn]:
                            h2Ref[chn][vname] = {}
                        hname = utils.get_hnames('Ref2D')(chn, vname, cstr)
                        h2Ref[chn][vname][cstr] = ROOT.TH2D(hname, '', *bin2D)
                        if vname not in h2Trig[chn]:
                            h2Trig[chn][vname] = {}
                        h2Trig[chn][vname][cstr] = {}

    t_in.SetBranchStatus('*', 0)
    _entries = utils.define_used_tree_variables(config_module.custom_cut)
    _entries += tuple(args.variables)
    for ientry in _entries:
        t_in.SetBranchStatus(ientry, 1)
    cmet = 0

    nentries = t_in.GetEntriesFast()
    for ientry,entry in enumerate(t_in):
        if ientry%10000==0:
             print('{} / {}'.format(ientry, nentries))

        # this is slow: do it once only
        entries = utils.dot_dict({x: getattr(entry, x) for x in _entries})

        w_mc     = entries.MC_weight
        w_pure   = entries.PUReweight
        w_l1pref = entries.L1pref_weight
        w_trig   = entries.trigSF
        w_idiso  = entries.IdSF_deep_2d
        w_jetpu  = entries.PUjetID_SF
        w_btag   = entries.bTagweightReshape
        
        if utils.is_nan(w_mc)     : w_mc=1
        if utils.is_nan(w_pure)   : w_pure=1
        if utils.is_nan(w_l1pref) : w_l1pref=1
        if utils.is_nan(w_trig)   : w_trig=1
        if utils.is_nan(w_idiso)  : w_idiso=1
        if utils.is_nan(w_jetpu)  : w_jetpu=1
        if utils.is_nan(w_btag)   : w_btag=1

        evt_weight = 1.#w_mc * w_pure * w_l1pref * w_trig * w_idiso * w_jetpu
        if args.isdata:
            if evt_weight != 1.:
                raise RuntimeError('The weight for data is {}!'.format(evt_weight))

        # if config_module.category != "baseline" and config_module.category != "boosted":
        #     # can be '-1' when no bjets are present
        #     evt_weight *= w_btag 

        if args.isdata:
            if abs(evt_weight) != 1.:
                raise RuntimeError('The weight for data is {}!'.format(evt_weight))

        sel = selection.EventSelection(entries, isdata=args.isdata, year=args.year,
                                       configuration=config_module)

        if not sel.sel_category(config_module.category):
            continue
        
        fill_var = {}
        for v in args.variables:
            fill_var[v] = {}
            for chn in args.channels:
                fill_var[v].update({chn: entries[v]})
                if fill_var[v][chn]>binedges[v][chn][-1]:
                    fill_var[v][chn]=binedges[v][chn][-1] # include overflow

        # whether the event passes cuts (1 and 2-dimensional) and single triggers
        pass_trigger, pcuts1D, pcuts2D = ({} for _ in range(3))
        for trig in config_module.triggers:
            pcuts2D[trig] = {}
        for trig in config_module.triggers:
            pass_trigger[trig] = sel.trigger_bits(trig)
            pcuts1D[trig] = {}
            for var in args.variables:
                pcuts1D[trig][var] = sel.var_cuts(trig, [var], args.nocut_dummy_str)

            if trig in config_module.pairs2D.keys():
                # combtrigs = tuple(x for x in triggercomb if trig in x)
                # for combtrig in combtrigs:
                # pcuts2D[joinNTC(combtrig)] = {}
                for j in config_module.pairs2D[trig]:
                    vname = utils.add_vnames(j[0],j[1])
                    for t in config_module.triggers:
                        pcuts2D[t][vname] = sel.var_cuts(t, [j[0], j[1]], args.nocut_dummy_str)

        #logic AND to intersect all triggers in this combination
        pass_trigger_intersection = {}
        for chn in args.channels:
            for tcomb in triggercomb[chn]:
                pass_trigger_intersection[joinNTC(tcomb)] = functools.reduce(
                    lambda x,y: x and y,
                    [ pass_trigger[x] for x in tcomb ] )

        for chn in args.channels:
            if utils.is_channel_consistent(chn, entries.pairType):

                # fill histograms for 1D efficiencies
                for j in args.variables:
                    binning1D = (nbins[j][chn], binedges[j][chn])

                    # The following is tricky, as we are considering, simultaneously:
                    # - all trigger intersection combinations
                    # - all cut combinations for each trigger combination (see 'main.cuts')
                        
                    # Logic AND to intersect all cuts for this trigger combination
                    # Each element will contain one possible cut combination
                    # for the trigger combination 'tcomb' being considered
                    for tcomb in triggercomb[chn]:
                        cstr = joinNTC(tcomb)

                        if not sel.check_inters_with_dataset(tcomb, chn, args.dataset):
                            continue
                        if not sel.dataset_cuts(tcomb, chn):
                            continue
                        if not sel.dataset_triggers(tcomb, chn, config_module.triggers,
                                                    args.dataset)[0]:
                            continue

                        # avoid underflow bin with negative weights crashing efficiency calculation
                        underflow = fill_var[j][chn] < binedges[j][chn][0]
                        fill_info = fill_var[j][chn], 1. if underflow else evt_weight
                        hRef[chn][j][cstr].Fill(*fill_info)
                        
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
                            if key not in hTrg[chn][j][cstr]:
                                base_str = utils.get_hnames('Trig1D')(chn,j,cstr)
                                htrig_name = utils.rewrite_cut_string(base_str, key)
                                hTrg[chn][j][cstr][key] = ROOT.TH1D(htrig_name, '', *binning1D)

                            if val and pass_trigger_intersection[cstr]:
                                hTrg[chn][j][cstr][key].Fill(*fill_info)

                # fill 2D efficiencies
                for onetrig in config_module.triggers:
                    if onetrig in config_module.pairs2D.keys():
                        combtrigs = tuple(x for x in triggercomb[chn] if onetrig in x)

                        for combtrig in combtrigs:
                            cstr = joinNTC(combtrig)
                            
                            if not sel.check_inters_with_dataset(combtrig, chn, args.dataset):
                                continue
                            if not sel.dataset_cuts(combtrig, chn):
                                continue
                            if not sel.dataset_triggers(combtrig, chn, config_module.triggers, args.dataset)[0]:
                                continue
                            
                            for j in config_module.pairs2D[onetrig]:
                                vname = utils.add_vnames(j[0],j[1])
                                # avoid underflow bin with negative weights crashing efficiency calculation
                                underflow = (fill_var[j[0]][chn] < binedges[j[0]][chn][0] or
                                             fill_var[j[1]][chn] < binedges[j[1]][chn][0])
                                fill_info = (fill_var[j[0]][chn], fill_var[j[1]][chn], 1. if underflow else evt_weight)

                                try:
                                    cuts_combinations = list(it.product(
                                        *(pcuts2D[atrig][vname].items() for atrig in combtrig) ))
                                except KeyError:
                                    print('CHECK: ', combtrig)
                                    raise

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
                                        h2Trig[chn][vname][cstr][key] = ROOT.TH2D(h2name, '', *bin2D)
                                        
                                    if val and pass_trigger_intersection[cstr]:
                                        h2Trig[chn][vname][cstr][key].Fill(*fill_info)
    
    file_id = ''.join(c for c in args.infile[-10:] if c.isdigit())
    outname = os.path.join(outdir, args.tprefix + args.sample + '_' + file_id + args.subtag + '.root')
    empty_files = True

    # normalize all histograms with luminosity and sum of weights
    if not args.isdata:
        norm = utils.get_lumi(args.year) / utils.total_sum_weights(args.infile, isdata=False)
        for chn in args.channels:
            for vname in hRef[chn].keys():
                for cstr in hRef[chn][vname].keys():
                    hRef[chn][vname][cstr].Scale(norm)
                    for key in hTrg[chn][vname][cstr].keys():
                        hTrg[chn][vname][cstr][key].Scale(norm)
            for vname in h2Ref[chn].keys():
                for cstr in h2Ref[chn][vname].keys():
                    h2Ref[chn][vname][cstr].Scale(norm)
                    for key in h2Trig[chn][vname][cstr].keys():
                        h2Trig[chn][vname][cstr][key].Scale(norm)

    # sanity check: reference counts must be always equal or larger than after applying some cut
    # remove negative weights under special conditions
    for chn in args.channels:
        for vname in hRef[chn].keys():
            for cstr in hRef[chn][vname].keys():
                old_ref = hRef[chn][vname][cstr].Integral()
                for key in hTrg[chn][vname][cstr].keys():
                    old_trg = hTrg[chn][vname][cstr][key].Integral()
                    neg_ref, neg_trg = False, False
                    for ibin in range(hRef[chn][vname][cstr].GetNbinsX()+1):
                        val_ref = hRef[chn][vname][cstr].GetBinContent(ibin)
                        val_trg = hTrg[chn][vname][cstr][key].GetBinContent(ibin)
                        if val_ref < 0.:
                            neg_ref = True
                            hRef[chn][vname][cstr].SetBinContent(ibin, 0.)
                        if val_trg < 0.:
                            neg_trg = True
                            hTrg[chn][vname][cstr][key].SetBinContent(ibin, 0.)

                        val_ref = hRef[chn][vname][cstr].GetBinContent(ibin)
                        val_trg = hTrg[chn][vname][cstr][key].GetBinContent(ibin)
                        if val_ref < val_trg:
                            print("Denominator: {} | Numerator: {}".format(val_ref, val_trg))
                            print("NegRef? {} | NegTrg? {}".format(neg_ref, neg_trg))
                            raise AssertionError()
                    if neg_ref:
                        new = hRef[chn][vname][cstr].Integral()
                        hRef[chn][vname][cstr].Scale(old_ref/new)
                    if neg_trg:
                        new = hTrg[chn][vname][cstr][key].Integral()
                        hTrg[chn][vname][cstr][key].Scale(old_trg/new)

    f_out = ROOT.TFile(outname, 'RECREATE')
    f_out.cd()
    for chn in args.channels:
        for j in args.variables:
            for tcomb in triggercomb[chn]:
                cstr = joinNTC(tcomb)
                if hRef[chn][j][cstr].GetEntries() > 0:
                    empty_files = False
                
                hRef[chn][j][cstr].Write( utils.get_hnames('Ref1D')(chn,j,cstr) )

                for khist,vhist in hTrg[chn][j][cstr].items():
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
parser = argparse.ArgumentParser(description='Producer trigger histograms.')

parser.add_argument('--binedges_fname', required=True, help='where the bin edges are stored')
parser.add_argument('--outdir', required=True, help='output directory')
parser.add_argument('--dataset', required=True,
                    help='Dataset name as provided by the user: MET, EG, ...')
parser.add_argument('--sample', required=True,
                    help='Process name as in SKIM directory')
parser.add_argument('--isdata', required=True, type=int, help='Whether it is data or MC')
parser.add_argument('--file', dest='infile', required=True, help='Full path of ROOT input file')
parser.add_argument('--year', required=True, type=str, choices=('2016', '2016APV', '2017', '2018'),
                    help='Data year: impact thresholds and selections.')
parser.add_argument('--subtag', required=True,
                    help='Additional (sub)tag to differentiate similar runs within the same tag.')
parser.add_argument('--tprefix', required=True, help='Targets name prefix.')
parser.add_argument('--channels', required=True, nargs='+', type=str,  
                    help='Select the channels over which the workflow will be run.' )
parser.add_argument('--variables', required=True, nargs='+', type=str,
                    help='Select the variables over which the workflow will be run.' )
parser.add_argument('--intersection_str', required=False, default=main.inters_str,
                    help='String used to represent set intersection between triggers.')
parser.add_argument('--nocut_dummy_str', required=True,
                    help='Dummy string associated to trigger histograms were no cuts are applied.')
parser.add_argument('--configuration', required=True,
                    help='Name of the configuration module to use.')
args = utils.parse_args(parser)

build_histograms(args)
