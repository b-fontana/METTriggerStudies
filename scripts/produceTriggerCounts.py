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

from writeHTCondorProcessingFiles import runTrigger_outputs_sample

from utils.utils import (
    generate_trigger_combinations,
    is_channel_consistent,
    join_name_trigger_intersection as joinNTC,
    LeafManager,
    pass_any_trigger,
    pass_selection_cuts,
    parse_args,
    pass_trigger_bits,
)

from luigi_conf import _triggers_custom

def getTriggerCounts(indir, outdir, sample, fileName,
                     channels, triggers,
                     subtag, tprefix, isdata ):
    # -- Check if outdir exists, if not create it
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    if not os.path.exists( os.path.join(outdir, sample) ):
        os.makedirs( os.path.join(outdir, sample) )
    outdir = os.path.join(outdir, sample)

    if not os.path.exists(fileName):
        raise ValueError('[' + os.path.basename(__file__) + '] {} does not exist.'.format(fileName))

    f_in = TFile(fileName)
    t_in = f_in.Get('HTauTauTree')

    triggercomb = {}
    for chn in channels:
        triggercomb[chn] = generate_trigger_combinations(chn, triggers)

    counter, counterRef = ({} for _ in range(2))
    for chn in channels:
        counter[chn] = {}
        counterRef.setdefault(chn, 0.)
        for tcomb in triggercomb[chn]:
            tcomb_str = joinNTC(tcomb)
            counter[chn][tcomb_str] = 0

    lf = LeafManager(fileName, t_in)
    
    for entry in range(0,t_in.GetEntries()):
        t_in.GetEntry(entry)

        if not pass_selection_cuts(lf):
            continue

        trig_bit = lf.get_leaf('pass_triggerbit')
        run = lf.get_leaf('RunNumber')
        if not pass_any_trigger(triggers, trig_bit, run, isdata=isdata):
            continue

        pass_trigger = {}
        for trig in triggers:
            pass_trigger[trig] = pass_trigger_bits(trig, trig_bit, run, isdata)

        for chn in channels:
            if is_channel_consistent(chn, lf.get_leaf('pairType')):
                counterRef[chn] += 1
                
                for tcomb in triggercomb[chn]:
                    pass_trigger_intersection = functools.reduce(
                        lambda x,y: x and y, #logic AND to join all triggers in this option
                        [ pass_trigger[x] for x in tcomb ]
                    )

                    if pass_trigger_intersection:
                        tcomb_str = joinNTC(tcomb)
                        counter[chn][tcomb_str] += 1
                                            
    file_id = ''.join( c for c in fileName[-10:] if c.isdigit() )

    proc_folder = os.path.dirname(fileName).split('/')[-1]
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    outName = os.path.join(outdir, tprefix + proc_folder + '_' + file_id + subtag + '.csv')
    print('Saving file {} at {} '.format(file_id, outName) )

    sep = ','
    with open(outName, 'w') as f:
        for chn in channels:
            f.write( 'Total' + sep + chn + sep + str(int(counterRef[chn])) + '\n' )
        for chn in channels:
            for tcomb in triggercomb[chn]:
                tcomb_str = joinNTC(tcomb)
                basestr = tcomb_str + sep + chn + sep
                f.write(basestr + str(int(counter[chn][tcomb_str])) + '\n' )
    
# -- Parse input arguments
parser = argparse.ArgumentParser(description='Command line parser')

parser.add_argument('--indir',       dest='indir',       required=True, help='SKIM directory')
parser.add_argument('--outdir',      dest='outdir',      required=True, help='output directory')
parser.add_argument('--sample',      dest='sample',      required=True, help='Process name as in SKIM directory')
parser.add_argument('--isdata',      dest='isdata',      required=True, help='Whether it is data or MC', type=int)
parser.add_argument('--file',        dest='fileName',    required=True, help='ID of input root file')
parser.add_argument('--subtag',      dest='subtag',      required=True,
                    help='Additional (sub)tag to differ  entiate similar runs within the same tag.')
parser.add_argument('--tprefix',     dest='tprefix',     required=True, help='Targets name prefix.')
parser.add_argument('--channels',    dest='channels',    required=True, nargs='+', type=str,  
                    help='Select the channels over which the workflow will be run.' )
parser.add_argument('--triggers',    dest='triggers',    required=True, nargs='+', type=str,
                    help='Select the triggers over which the workflow will be run.' )
parser.add_argument('--debug', action='store_true', help='debug verbosity')
args = parse_args(parser)

getTriggerCounts( args.indir, args.outdir, args.sample, args.fileName,
                  args.channels, args.triggers,
                  args.subtag, args.tprefix, args.isdata )
