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
from inclusion.config import main

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
        regex = re.compile( args.tprefix + args.sample + '.+_[0-9]{1,5}' + args.subtag + '.csv' )
        walk_path = os.path.join(args.indir, args.sample)
        for root, _, files in os.walk( walk_path ):
            for afile in files:
                if regex.match( os.path.basename(afile) ):
                    inputs_join.append( os.path.join(root, afile) )
            are_there_files(files, regex)

    aggr_outs = []
    sep = ','
    reftrigs = {}
    c_ref, c_inters = ({} for _ in range(2))

    if args.aggr:
        regex2 = re.compile(args.tprefix + '(.+)_Sum.*' + args.subtag + '_' + args.channel + '.csv')    
        for afile in inputs_join:
            filetype = regex2.findall(afile)

            with open(afile, 'r') as f:
                for line in f.readlines():
                    if line.strip(): #ignore empty lines
                        aggr_outs.append(filetype[0] + sep + line)


    else:
        for afile in inputs_join:
            filetype = regex2.findall(afile)

            with open(afile, 'r') as f:
                for line in f.readlines():
                    if not line.strip(): #ignore empty lines
                        continue
                        
                    split_line = [x.replace('\n', '') for x in line.split(sep)]
                    if all( not x for x in split_line ):
                        continue
                            
                    title, comb, chn, reftrig, count = split_line
                        
                    if chn != args.channel:
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


    if args.aggr:
        suboutdir = os.path.join(args.outdir, args.channel, 'Counts_' + args.dataset_name)
        if not os.path.exists(suboutdir):
            os.makedirs(suboutdir)
        outs = os.path.join(suboutdir, 'table.csv')
    else:
        outputs_csv = args.outfile_counts
        pref, suf = outputs_csv.split('.')
        outs = outputs_csv

    if args.aggr:
        with open(outs, 'w') as fcsv:
            for il,l in enumerate(aggr_outs):
                ftype, atype, comb, ref, c, eff = [x.replace('\n', '') for x in l.split(sep)]
                
                if ( atype == '' and comb == '' and ref == '' and
                     c == '' and eff == '' ):
                    continue

                if atype != 'Type':
                    fcsv.write(l)

                if atype == 'Type' and il==0:
                    newline = sep.join(('File Type', atype, comb, ref, c, eff)) + '\n'
                    fcsv.write(newline)

    else:
        with open(outs, 'w') as fcsv:
            ref_combs, ref_vals = ([] for _ in range(2))
            int_combs, int_vals = ([] for _ in range(2))
            refs = []
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
            refs      = np.array(refs)g
     
            #remove zeros
            zeromask = int_vals == 0
            ref_vals  = ref_vals[~zeromask]
            ref_combs = ref_combs[~zeromask]
            int_vals  = int_vals[~zeromask]
            int_combs = int_combs[~zeromask]
            refs      = refs[~zeromask]
     
            line = sep.join(('Type', 'Trigger Intersection', 'Reference', 'Counts', 'Efficiency\n'))
            fcsv.write(line)

            # calculate efficiencies
            effs = []
            gzip = zip(ref_vals, ref_combs, int_vals, int_combs, refs)
            for refv, refc, intv, intc, refref in gzip:
                effs.append(float(intv) / float(refv))

            effs = np.sort(np.array(effs))[::-1]
            
            # sort other columns following efficiencies descending order
            idxs_sort = effs.argsort(axis=None)[::-1]
            ref_vals  = ref_vals[idxs_sort]
            ref_combs = ref_combs[idxs_sort]
            int_vals  = int_vals[idxs_sort]
            int_combs = int_combs[idxs_sort]

            gzip = zip(effs, ref_vals, ref_combs, int_vals, int_combs, refs)
            for eff, refv, refc, intv, intc, refref in gzip:
                newline = sep.join(('Numerator',
                                    str(intc).replace(main.inters_str, '  AND  '),
                                    refref, str(intv), str(round(eff,4)))) + '\n'
                fcsv.write(newline)
     
                newline = sep.join(('Denominator',
                                    str(refc).replace(main.inters_str, '  AND  '),
                                    refref, str(refv), '1\n' )
                fcsv.write(newline)

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

    parser.add_argument('--sample',     dest='sample',       required=False,
                        help='Process name as in SKIM directory. Used for the first step only.')
    
    parser.add_argument('--infile_counts', dest='infile_counts', required=False, nargs='+', type=str,
                        help='Name of input csv files with counts. Used for the aggregation step only.')
    parser.add_argument('--outfile_counts', dest='outfile_counts', help='Name of output csv files with counts.')
    parser.add_argument('--channel', dest='channel', required=False, help='Channel to be used for the aggregation.')
    args = utils.parse_args(parser)
    
    add_trigger_counts(args)
