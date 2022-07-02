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
    reftrigs = {}
    c_ref, c_inters = ({} for _ in range(2))
    for afile in inputs_join:
        with open(afile, 'r') as f:
            for line in f.readlines():
                if line.strip(): #ignore empty lines
                    title, comb, chn, reftrig, count = [x.replace('\n', '') for x in line.split(sep)]
                    print(title, comb, chn, reftrig, count)

                    if chn not in reftrigs:
                        reftrigs[chn] = {}
                    if comb not in reftrigs[chn]:
                        reftrigs[chn][comb] = reftrig
                        
                    if title == 'Reference':
                        c_ref.setdefault(chn, {})
                        c_ref[chn].setdefault(comb, 0)
                        c_ref[chn][comb] += int(count)
                    elif title == 'Intersection':
                        c_inters.setdefault(chn, {})
                        c_inters[chn].setdefault(comb, 0)
                        c_inters[chn][comb] += int(count)
                    else:
                        mes = 'Column {} is not supported.'
                        raise ValueError(mes.format(mes))

    outs = dict()
    outputs_csv = args.outfile_counts
    channels = list(c_inters.keys())
    for chn in channels:
        if args.aggr:
            suboutdir = os.path.join(args.outdir, chn, 'Counts_' + args.dataset_name)
            if not os.path.exists(suboutdir):
                os.makedirs(suboutdir)
            outs[chn] = os.path.join(suboutdir, 'table.csv')
        else:
            pref, suf = outputs_csv.split('.')
            outs[chn] = pref + '_' + chn + '.' + suf
        
    for ic,chn in enumerate(channels):
        with open(outs[chn], 'w') as fcsv:
            ref_combs, ref_vals = ([] for _ in range(2))
            int_combs, int_vals = ([] for _ in range(2))
            for comb, val in c_inters[chn].items():
                ref_combs.append('Reference_' + comb)
                ref_vals.append(c_ref[chn][comb])
                int_combs.append(comb)
                int_vals.append(val)

            ref_combs = np.array(ref_combs)
            ref_vals  = np.array(ref_vals)
            int_combs = np.array(int_combs)
            int_vals  = np.array(int_vals)
            
            #sort
            ref_vals, ref_combs = (np.array(t[::-1])
                                   for t in zip(*sorted(zip(ref_vals, ref_combs))))
            int_vals, int_combs = (np.array(t[::-1])
                                   for t in zip(*sorted(zip(int_vals, int_combs))))

            #remove zeros
            zeromask = ref_vals == 0
            ref_vals     = ref_vals[~zeromask]
            ref_combs    = ref_combs[~zeromask]
            inters_vals  = inters_vals[~zeromask]
            inters_combs = inters_combs[~zeromask]
            if not all(inters_vals):
                mes = 'At least one of the trigger intersections has zero entries.'
                mes += ' This should not happen after filtering using the reference '
                mes += '(denominator of the efficiency).\n'
                mess += 'Intersections: {}\n'
                mess += 'Counts: {}\n'
                raise RuntimeError(mes.format(inters_combs, inters_vals))
            #refval = vals[trigs.tolist().index('Total')]

            line = sep.join('Trigger Intersection', 'Reference', 'Counts', 'Efficiency\n')
            fcsv.write(line)
            for refv, refc, intv, intc in zip(ref_vals, ref_combs, inters_vals, inters_combs):
                eff = float(intv) / float(refv)
                newline = ( str(refc).replace('_PLUS_', '  AND  ') + sep + str(i) +
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
