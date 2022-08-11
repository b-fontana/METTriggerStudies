# coding: utf-8

_all_ = [ "get_trig_counts" ]

import re
import os
import sys
import functools
import argparse
import fnmatch
from array import array
from ROOT import TFile

import sys
sys.path.append( os.environ['PWD'] ) 

import inclusion
from inclusion import selection
from inclusion.utils import utils
from inclusion.utils.utils import join_name_trigger_intersection as joinNTC

def get_trig_counts(outdir, dataset, sample, filename,
                    channels, triggers,
                    subtag, tprefix, isdata):
    # -- Check if outdir exists, if not create it
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    if not os.path.exists( os.path.join(outdir, sample) ):
        os.makedirs( os.path.join(outdir, sample) )
    outdir = os.path.join(outdir, sample)

    if not os.path.exists(filename):
        mes = '[' + os.path.basename(__file__) + '] {} does not exist.'.format(filename)
        raise ValueError(mes)

    f_in = TFile(filename)
    t_in = f_in.Get('HTauTauTree')

    triggercomb = {}
    for chn in channels:
        triggercomb[chn] = utils.generate_trigger_combinations(chn, triggers)

    c_ref, c_inters = ({} for _ in range(2))
    for chn in channels:
        c_ref[chn], c_inters[chn] = ({} for _ in range(2))
        for tcomb in triggercomb[chn]:
            tcomb_str = joinNTC(tcomb)
            c_ref[chn][tcomb_str] = 0
            c_inters[chn][tcomb_str] = 0

    t_in.SetBranchStatus('*', 0)
    _entries = ('pass_triggerbit', 'RunNumber',
                'HHKin_mass', 'pairType', 'dau1_eleMVAiso', 'dau1_iso', 'dau1_deepTauVsJet', 'dau2_deepTauVsJet',
                'nleps', 'nbjetscand', 'tauH_SVFIT_mass', 'bH_mass_raw',)
    for ientry in _entries:
        t_in.SetBranchStatus(ientry, 1)

    nentries = t_in.GetEntriesFast()
    for ientry,entry in enumerate(t_in):
        if ientry%10000==0:
             print('Processed {} entries out of {}.'.format(ientry, nentries))

        # this is slow: do it once only
        entries = utils.dot_dict({x: getattr(entry, x) for x in _entries})
        
        sel = selection.EventSelection(entries, dataset, isdata)
        
        pass_trigger = {}
        for trig in triggers:
            pass_trigger[trig] = sel.trigger_bits(trig)

        for chn in channels:
            if utils.is_channel_consistent(chn, entries.pairType):

                for tcomb in triggercomb[chn]:                
                    pass_trigger_intersection = functools.reduce(
                        lambda x,y: x and y, #logic AND to join all triggers in this option
                        [ pass_trigger[x] for x in tcomb ]
                    )

                    if not sel.dataset_cuts(tcomb, chn):
                        continue
                    if not sel.dataset_triggers(triggers, tcomb, chn)[0]:
                        continue
                    if not sel.match_inters_with_dataset(tcomb, chn):
                        continue
                    
                    tcomb_str = joinNTC(tcomb)
                    c_ref[chn][tcomb_str] += 1

                    if pass_trigger_intersection:
                        c_inters[chn][tcomb_str] += 1
                                            
    file_id = ''.join( c for c in filename[-10:] if c.isdigit() )

    proc_folder = os.path.dirname(filename).split('/')[-1]
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    outName = os.path.join(outdir, tprefix + proc_folder + '_' + file_id + subtag + '.csv')
    print('Saving file {} at {} '.format(file_id, outName) )

    sep = ','
    with open(outName, 'w') as f:
        for chn in channels:
            for tcomb in triggercomb[chn]:
                try:
                    reftrig = sel.dataset_triggers(triggers, tcomb, chn)[1]
                except OverflowError:
                    continue
                
                reftrig = joinNTC(reftrig)
                
                tcomb_str = joinNTC(tcomb)
                basestr = tcomb_str + sep + chn + sep + reftrig + sep
                f.write( 'Reference' + sep + basestr + str(int(c_ref[chn][tcomb_str])) + '\n' )
                f.write( 'Intersection' + sep + basestr + str(int(c_inters[chn][tcomb_str])) + '\n' )
    
# -- Parse input arguments
parser = argparse.ArgumentParser(description='Command line parser')

parser.add_argument('--outdir',      dest='outdir',      required=True, help='output directory')
parser.add_argument('--dataset', dest='dataset', required=True,
                    help='Dataset name as provided by the user: MET, EG, ...')
parser.add_argument('--sample',      dest='sample',      required=True, help='Process name as in SKIM directory')
parser.add_argument('--isdata',      dest='isdata',      required=True, help='Whether it is data or MC', type=int)
parser.add_argument('--file',        dest='filename',    required=True, help='ID of input root file')
parser.add_argument('--subtag',      dest='subtag',      required=True,
                    help='Additional (sub)tag to differ  entiate similar runs within the same tag.')
parser.add_argument('--tprefix',     dest='tprefix',     required=True, help='Targets name prefix.')
parser.add_argument('--channels',    dest='channels',    required=True, nargs='+', type=str,  
                    help='Select the channels over which the workflow will be run.' )
parser.add_argument('--triggers',    dest='triggers',    required=True, nargs='+', type=str,
                    help='Select the triggers over which the workflow will be run.' )
parser.add_argument('--debug', action='store_true', help='debug verbosity')
args = utils.parse_args(parser)

get_trig_counts( outdir=args.outdir, dataset=args.dataset, sample=args.sample,
                 filename=args.filename, channels=args.channels, triggers=args.triggers,
                 subtag=args.subtag, tprefix=args.tprefix, isdata=args.isdata )
