# coding: utf-8

_all_ = [ 'add_trigger_counts' ]

import os
import sys
parent_dir = os.path.abspath(__file__ + 3 * '/..')
sys.path.insert(0, parent_dir)

import re
import numpy as np
import argparse

import inclusion
from inclusion.utils import utils

@utils.set_pure_input_namespace
def add_trigger_counts(args):
    def are_there_files(files, regex):
        if len(files)==0:
            m =  '\nThe walk was performed in {} .\n'.format(os.path.join(args.indir))
            m += 'The regular expression was {}.\n'.format(regex.pattern)
            m += 'The regular expression could not retrieve any files.'
            raise ValueError(m)

    inputs_join = []
    if args.aggr:
        regex = re.compile(args.tprefix + '.+_Sum.*' + args.subtag + '_' + args.channel + '.csv')
        walk_path = args.indir

        for root, d, files in os.walk(walk_path, topdown=True):
            if root[len(walk_path):].count(os.sep) < 1:

                for afile in files:
                    if regex.match( os.path.basename(afile) ):
                        afile_full = os.path.join(root, afile)
                        if afile_full in args.infile_counts:
                            inputs_join.append( afile_full )

        are_there_files(files, regex)
        assert set(inputs_join) == set(args.infile_counts)

    else:
        regex = re.compile( args.tprefix + '.+_[0-9]{1,5}' + args.subtag + '.csv' )
        walk_path = os.path.join(args.indir, args.sample)
        for root, _, files in os.walk( walk_path ):
            for afile in files:
                if regex.match( os.path.basename(afile) ):
                    inputs_join.append( os.path.join(root, afile) )
            are_there_files(files, regex)
            
    regex2 = re.compile(args.tprefix + '(.+)_Sum.*' + args.subtag + '_' + args.channel + '.csv')
    aggr_outs = []
    sep = ','
    reftrigs = {}
    c_ref, c_inters = ({} for _ in range(2))
    
    for afile in inputs_join:
        filetype = regex2.findall(afile)

        with open(afile, 'r') as f:
            for line in f.readlines():
                if line.strip(): #ignore empty lines
                    if args.aggr:
                        aggr_outs.append( filetype[0] + sep + line )

                    else:
                        title, comb, chn, reftrig, count = [x.replace('\n', '') for x in line.split(sep)]

                        if chn != args.channel:
                            continue
                        if ( title == '' and comb == '' and chn == '' and
                            reftrig == '' and count == '' ):
                            continue

                        if comb not in reftrigs:
                            reftrigs[comb] = reftrig                        

                        if title == 'Reference':
                            c_ref.setdefault(comb, 0)
                            c_ref[comb] += int(count)
                        elif title == 'Intersection':
                            c_inters.setdefault(comb, 0)
                            c_inters[comb] += int(count)
                        else:
                            mes = 'Column {} is not supported.'
                            raise ValueError(mes.format(mes))

    outputs_csv = args.outfile_counts
    if args.aggr:
        suboutdir = os.path.join(args.outdir, args.channel, 'Counts_' + args.dataset_name)
        if not os.path.exists(suboutdir):
            os.makedirs(suboutdir)
        outs = os.path.join(suboutdir, 'table.csv')
    else:
        pref, suf = outputs_csv.split('.')
        outs = outputs_csv

    if args.aggr:
        with open(outs, 'w') as fcsv:
            for il,l in enumerate(aggr_outs):
                ftype, atype, comb, ref, c, eff = [x.replace('\n', '') for x in l.split(sep)]
                
                if ( atype == '' and comb == '' and ref == '' and
                     c == '' and eff == '' ):
                    fcsv.write(',,,,,\n')
                    continue

                if atype != 'Type':
                    fcsv.write(l)

                if atype == 'Type' and il==0:
                    newline = ( 'File Type' + sep + atype + sep + comb + sep +
                                ref + sep + c + sep + eff + '\n')
                    fcsv.write(newline)

    else:
        with open(outs, 'w') as fcsv:
            ref_combs, ref_vals = ([] for _ in range(2))
            int_combs, int_vals = ([] for _ in range(2))
            references, ftypes = ([] for _ in range(2))
            for comb, val in c_inters.items():
                ref_combs.append(comb)
                ref_vals.append(c_ref[comb])
                int_combs.append(comb)
                int_vals.append(val)
                references.append(reftrigs[comb])
     
            ref_combs = np.array(ref_combs)
            ref_vals  = np.array(ref_vals)
            int_combs = np.array(int_combs)
            int_vals  = np.array(int_vals)
                
            #sort
            gzip = zip(ref_vals, ref_combs, int_vals, int_combs)
            ref_vals, ref_combs, int_vals, int_combs = (np.array(t[::-1])
                                                        for t in zip(*sorted(gzip)))
     
            #remove zeros
            zeromask = int_vals == 0
            ref_vals  = ref_vals[~zeromask]
            ref_combs = ref_combs[~zeromask]
            int_vals  = int_vals[~zeromask]
            int_combs = int_combs[~zeromask]
            if not all(ref_vals):
                mes = 'At least one of the trigger references has zero entries.'
                mes += ' This should not happen after filtering using the intersections '
                mes += '(denominator of the efficiency).\n'
                mes += 'Intersections:\n  {}\n'
                mes += 'Intersection counts:\n {}\n'
                mes += 'References:\n {}\n'
                mes += 'References counts:\n {}\n'
                raise RuntimeError(mes.format(int_combs, int_vals, ref_combs, ref_vals))
            #refval = vals[trigs.tolist().index('Total')]
     
            line = sep.join(('Type', 'Trigger Intersection', 'Reference', 'Counts', 'Efficiency\n'))
            fcsv.write(line)
            gzip = zip(ref_vals, ref_combs, int_vals, int_combs, references)
            for refv, refc, intv, intc, refref in gzip:
                eff = float(intv) / float(refv)

                newline = ( 'Numerator' + sep +
                            str(intc).replace('_PLUS_', '  AND  ') + sep +
                            refref + sep + str(intv) +
                            sep + str(round(eff,4)) + '\n' )
                fcsv.write(newline)
     
                newline = ( 'Denominator' + sep + str(refc).replace('_PLUS_', '  AND  ')
                            + sep + refref + sep + str(refv) +
                            sep + '1\n' )
                fcsv.write(newline)
     
                fcsv.write(',,,,\n')
                fcsv.write(',,,,\n')        
                fcsv.write('\n')

    print('Save file {}.'.format(outs))
        
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Command line parser')

    parser.add_argument('--indir', dest='indir', required=True, help='SKIM directory')
    parser.add_argument('--outdir', dest='outdir', required=True, help='output directory')
    parser.add_argument('--subtag', dest='subtag',      required=True,
                        help='Additional (sub)tag to differ  entiate similar runs within the same tag.')
    parser.add_argument('--tprefix',     dest='tprefix',     required=True, help='Targets name prefix.')
    parser.add_argument('--aggregation_step', dest='aggr',   required=True, type=int,
                        help='Whether to run the sample aggregation step or the "per sample step"')
    parser.add_argument('--dataset_name', dest='dataset_name', required=True,
                        help='Name of the dataset being used.')
    parser.add_argument('--debug', action='store_true', help='debug verbosity')

    parser.add_argument('--sample',     dest='sample',       required=False,
                        help='Process name as in SKIM directory. Used for the first step only.')
    
    parser.add_argument('--infile_counts', dest='infile_counts', required=False, nargs='+', type=str,
                        help='Name of input csv files with counts. Used for the aggregation step only.')
    parser.add_argument('--outfile_counts', dest='outfile_counts', required=True,
                        help='Name of output csv files with counts.')
    parser.add_argument('--channel', dest='channel', required=False, help='Channel to be used for the aggregation.')
    args = utils.parse_args(parser)
    
    add_trigger_counts(args)
