# coding: utf-8

_all_ = [ 'get_trig_counts' ]

import os
import sys
parent_dir = os.path.abspath(__file__ + 3 * '/..')
sys.path.insert(0, parent_dir)
import glob
    
import inclusion
from inclusion import selection
from inclusion.utils import utils
from inclusion.utils.utils import join_name_trigger_intersection as joinNTC

import functools
import argparse
import importlib

import ROOT

def get_trig_counts(args):
    # -- Check if outdir exists, if not create it
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
    if not os.path.exists( os.path.join(args.outdir, args.sample) ):
        os.makedirs( os.path.join(args.outdir, args.sample) )
    outdir = os.path.join(args.outdir, args.sample)
    
    if not os.path.exists(args.filename):
        mes = '[' + os.path.basename(__file__) + '] {} does not exist.'.format(args.filename)
        raise ValueError(mes)

    config_module = importlib.import_module(args.configuration)
    
    xsec_norm = utils.total_cross_section(args.filename, args.isdata)
    
    f_in = ROOT.TFile(args.filename)
    t_in = f_in.Get('HTauTauTree')

    triggercomb = {}
    for chn in args.channels:
        triggercomb[chn] = utils.generate_trigger_combinations(chn, config_module.triggers,
                                                               config_module.exclusive)

    c_ref , c_inters  = ({} for _ in range(2))
    w_ref , w_inters  = ({} for _ in range(2))
    w2_ref, w2_inters = ({} for _ in range(2))
    for chn in args.channels:
        c_ref[chn], c_inters[chn] = ({} for _ in range(2))
        w_ref[chn], w_inters[chn] = ({} for _ in range(2))
        w2_ref[chn], w2_inters[chn] = ({} for _ in range(2))
        for tcomb in triggercomb[chn]:
            tstr = joinNTC(tcomb)
            c_ref[chn][tstr] = 0
            c_inters[chn][tstr] = 0
            w_ref[chn][tstr] = 0.
            w_inters[chn][tstr] = 0.
            w2_ref[chn][tstr] = 0.
            w2_inters[chn][tstr] = 0.

    t_in.SetBranchStatus('*', 0)
    _entries = utils.define_used_tree_variables(config_module.custom_cut)
    for ientry in _entries:
        t_in.SetBranchStatus(ientry, 1)

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

        evt_weight = w_mc * w_pure * w_l1pref * w_trig * w_idiso * w_jetpu * w_btag
        if args.isdata:
            assert evt_weight == 1.

        sel = selection.EventSelection(entries, args.isdata, configuration=config_module)

        if not sel.sel_category(config_module.category):
            continue
        
        pass_trigger = {}
        for trig in config_module.triggers:
            pass_trigger[trig] = sel.trigger_bits(trig)

        for chn in args.channels:
            if utils.is_channel_consistent(chn, entries.pairType):
                    
                for tcomb in triggercomb[chn]:                
                    pass_trigger_intersection = functools.reduce(
                        lambda x,y: x and y, #logic AND to join all triggers in this option
                        [ pass_trigger[x] for x in tcomb ]
                    )

                    if not sel.check_inters_with_dataset(tcomb, chn, args.dataset):
                        continue
                    if not sel.dataset_cuts(tcomb, chn):
                        continue
                    if not sel.dataset_triggers(tcomb, chn, config_module.triggers, args.dataset)[0]:
                        continue

                    tstr = joinNTC(tcomb)
                    c_ref[chn][tstr] += 1
                    w_ref[chn][tstr] += evt_weight
                    w2_ref[chn][tstr] += evt_weight*evt_weight

                    if pass_trigger_intersection:
                        c_inters[chn][tstr] += 1
                        w_inters[chn][tstr] += evt_weight
                        w2_inters[chn][tstr] += evt_weight*evt_weight
                                            

    file_id = ''.join( c for c in args.filename[-10:] if c.isdigit() )

    proc_folder = os.path.dirname(args.filename).split('/')[-1]
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    outName = os.path.join(outdir, args.tprefix + proc_folder + '_' + file_id + args.subtag + '.csv')
    print('Saving file {} at {} '.format(file_id, outName) )

    sep = ','
    with open(outName, 'w') as f:
        for chn in args.channels:
            for tcomb in triggercomb[chn]:
                try:
                    reftrig = sel.dataset_triggers(tcomb, chn, config_module.triggers, args.dataset)[1]
                except OverflowError:
                    continue
                
                reftrig = joinNTC(reftrig)
                
                tstr = joinNTC(tcomb)
                basestr = sep.join((tstr, chn, reftrig))

                if not args.isdata:
                    norm_factor = utils.get_lumi(args.year) / utils.total_sum_weights(args.file, isdata=False)
                    w_ref[chn][tstr] *= norm_factor
                    w_inters[chn][tstr] *= norm_factor
                    w2_ref[chn][tstr] *= norm_factor**2
                    w2_inters[chn][tstr] *= norm_factor**2

                counts_ref = str(int(c_ref[chn][tstr]))
                counts_int = str(int(c_inters[chn][tstr]))
                f.write( sep.join(('Reference', basestr, counts_ref)) + '\n' )
                f.write( sep.join(('Intersection', basestr, counts_int)) + '\n' )

                weights_ref = str(float(w_ref[chn][tstr]))
                weights_int = str(float(w_inters[chn][tstr]))
                f.write( sep.join(('Reference_weighted', basestr, weights_ref)) + '\n' )
                f.write( sep.join(('Intersection_weighted', basestr, weights_int)) + '\n' )

                w2_ref_str = str(float(w2_ref[chn][tstr]))
                w2_int_str = str(float(w2_inters[chn][tstr]))
                f.write( sep.join(('Reference_w2', basestr, w2_ref_str)) + '\n' )
                f.write( sep.join(('Intersection_w2', basestr, w2_int_str)) + '\n' )

# -- Parse input arguments
parser = argparse.ArgumentParser(description='Produce trigger counts.')

parser.add_argument('--outdir',      dest='outdir',      required=True, help='output directory')
parser.add_argument('--dataset', dest='dataset', required=True,
                    help='Dataset name as provided by the user: MET, EG, ...')
parser.add_argument('--sample',      dest='sample',      required=True, help='Process name as in SKIM directory')
parser.add_argument('--isdata',      dest='isdata',      required=True, help='Whether it is data or MC', type=int)
parser.add_argument('--year', required=True, type=str, choices=('2016', '2016APV', '2017', '2018'),
                    help='Data year: impact thresholds and selections.')
parser.add_argument('--file',        dest='filename',    required=True, help='ID of input root file')
parser.add_argument('--subtag',      dest='subtag',      required=True,
                    help='Additional (sub)tag to differ  entiate similar runs within the same tag.')
parser.add_argument('--tprefix',     dest='tprefix',     required=True, help='Targets name prefix.')
parser.add_argument('--channels',    dest='channels',    required=True, nargs='+', type=str,  
                    help='Select the channels over which the workflow will be run.' )
parser.add_argument('--configuration', dest='configuration', required=True,
                    help='Name of the configuration module to use.')
args = utils.parse_args(parser)

get_trig_counts(args)
