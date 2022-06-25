import os
import re
import numpy as np
import argparse
import sys
sys.path.append(os.path.join(os.environ['CMSSW_BASE'], 'src', 'METTriggerStudies'))
from pathlib import Path

from utils import utils

@utils.set_pure_input_namespace
def addTriggerCounts(args):
    def are_there_files(files, regex):
        if len(files)==0:
            m =  '\nThe walk was performed in {} .\n'.format(os.path.join(args.indir, smpl))
            m += 'The regular expression was {}.\n'.format(regex.pattern)
            m += 'The regular expression could not retrieve any files.'
            raise ValueError(m)

    inputs_join = []
    if args.aggr:
        regex = re.compile(args.tprefix + '.+_Sum.*' + args.subtag + '.csv')
        walk_path = args.indir

        for root, d, files in os.walk(walk_path, topdown=True):
            if root[len(walk_path):].count(os.sep) < 1:

                for afile in files:
                    if regex.match( os.path.basename(afile) ):
                        afile_full = os.path.join(root, afile)
                        if afile_full in args.infile_counts:
                            print('FULL: ', afile_full)
                            inputs_join.append( afile_full )

        are_there_files(files, regex)
        print(args.infile_counts)
        assert set(inputs_join) == set(args.infile_counts)

    else:
        regex = re.compile( args.tprefix + '.+_[0-9]{1,5}' + args.subtag + '.csv' )
        walk_path = os.path.join(args.indir, args.sample)
        for root, _, files in os.walk( walk_path ):
            for afile in files:
                if regex.match( os.path.basename(afile) ):
                    inputs_join.append( os.path.join(root, afile) )
            are_there_files(files, regex)

    sep = ','
    counter, counter_ref = ({} for _ in range(2))
    for afile in inputs_join:
        with open(afile, 'r') as f:
            for line in f.readlines():
                if line.strip(): #ignore empty lines
                    trig, chn, count = [x.replace('\n', '') for x in line.split(sep)]
                    print(trig, chn, count)
                    
                    if trig != 'Total':
                        counter.setdefault(chn, {})
                        counter[chn].setdefault(trig, 0)
                        counter[chn][trig] += int(count)
                    else:
                        counter_ref.setdefault(chn, 0)
                        counter_ref[chn] += int(count)

    outs = dict()
    outputs_csv = args.outfile_counts
    channels = list(counter.keys())
    for chn in channels:
        if args.aggr:
            suboutdir = os.path.join(args.outdir, chn, 'Counts_' + args.dataset_name)
            if not os.path.exists(suboutdir):
                os.makedirs(suboutdir)
            outs[chn] = os.path.join(suboutdir, 'table.csv')
        else:
            pref, suf = outputs_csv.split('.')
            outs[chn] = pref + '_' + chn + '.' + suf
        
    for ic,chn in enumerate(counter):

        with open(outs[chn], 'w') as fcsv:
            trigs, vals = ([] for _ in range(2))
            trigs.append('Total')
            vals.append(counter_ref[chn])
            for trig,val in counter[chn].items():
                trigs.append(trig)
                vals.append(val)

            trigs = np.array(trigs)
            vals = np.array(vals)
            
            #sort
            vals, trigs = (np.array(t[::-1]) for t in zip(*sorted(zip(vals, trigs))))

            #remove zeros
            zeromask = vals == 0
            vals = vals[~zeromask]
            trigs = trigs[~zeromask]
            refval = vals[trigs.tolist().index('Total')]

            fcsv.write('Trigger Intersection' + sep + 'Counts' + sep + 'Efficiency' + '\n')
            for i,j in zip(vals,trigs):
                if i != 0: #do not print the padding
                    #remove the extra info after the line break
                    eff = float(i) / refval
                    newline = ( str(j).replace('_PLUS_', '  AND  ') + sep + str(i) +
                                sep + str(round(eff,4)) + '\n' )
                    fcsv.write(newline)

            fcsv.write('\n')
                    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Command line parser')

    parser.add_argument('--indir',       dest='indir',       required=True, help='SKIM directory')
    parser.add_argument('--outdir',      dest='outdir',      required=True, help='output directory')
    parser.add_argument('--subtag',      dest='subtag',      required=True,
                        help='Additional (sub)tag to differ  entiate similar runs within the same tag.')
    parser.add_argument('--tprefix',     dest='tprefix',     required=True, help='Targets name prefix.')
    parser.add_argument('--aggregation_step', dest='aggr',   required=True, type=int,
                        help='Whether to run the sample aggregation step or the "per sample step"')
    parser.add_argument('--dataset_name', dest='dataset_name', required=True,
                        help='Name of the dataset being used.')
    parser.add_argument('--debug', action='store_true', help='debug verbosity')

    parser.add_argument('--sample',     dest='sample',       required=False,
                        help='Process name as in SKIM directory. Used for the first step only.')
    
    parser.add_argument('--infile_counts',  dest='infile_counts', required=False, nargs='+', type=str,
                        help='Name of input csv files with counts. Used for the aggrgeation step only.')
    parser.add_argument('--outfile_counts', dest='outfile_counts', required=True,
                        help='Name of output csv files with counts.')

    args = parser.parse_args()
    utils.print_configuration(args)
    
    addTriggerCounts(args)
