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

import ROOT

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
        regex = re.compile( args.tprefix + args.sample + '.*_[0-9]{1,5}' + args.subtag + '.csv' )
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
    w_ref, w_inters = ({} for _ in range(2))
    w2_ref, w2_inters = ({} for _ in range(2))

    if args.aggr:
        regex2 = re.compile(args.tprefix + '(.+)_Sum.*' + args.subtag + '_' + args.channel + '.csv')
        for afile in inputs_join:
            filetype = regex2.findall(afile)

            with open(afile, 'r') as f:
                for line in f.readlines():
                    if not line.strip(): #ignore empty lines
                        continue

                    split_line = [x.replace('\n', '') for x in line.split(sep)]
                    if all( not x for x in split_line ):
                        continue

                    new_line = line.replace('\n', '')
                    aggr_outs.append(filetype[0] + sep + new_line)

    else: # paired with 'if args.aggr:'
        for afile in inputs_join:
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
                    elif title == 'Reference_weighted':
                        w_ref.setdefault(comb, 0.)
                        w_ref[comb] += float(count)
                    elif title == 'Intersection_weighted':
                        w_inters.setdefault(comb, 0.)
                        w_inters[comb] += float(count)
                    elif title == 'Reference_w2':
                        w2_ref.setdefault(comb, 0.)
                        w2_ref[comb] += float(count)
                    elif title == 'Intersection_w2':
                        w2_inters.setdefault(comb, 0.)
                        w2_inters[comb] += float(count)
                    else:
                        mes = 'Column {} is not supported.'
                        raise ValueError(mes.format(mes))


    if args.aggr:
        table_name = 'table.csv'
        sub = os.path.join(args.outdir, args.channel, 'Tables')
        sub_c = os.path.join(sub, 'Counts_' + args.dataset_name)
        sub_w = os.path.join(sub, 'Weights_' + args.dataset_name)
        sub_c_squash = os.path.join(sub, 'CountsSquash_' + args.dataset_name)
        sub_w_squash = os.path.join(sub, 'WeightsSquash_' + args.dataset_name)
        utils.create_single_dir(sub_c)
        utils.create_single_dir(sub_w)
        utils.create_single_dir(sub_c_squash)
        utils.create_single_dir(sub_w_squash)
        outs_c = os.path.join(sub_c, table_name)
        outs_w = os.path.join(sub_w, table_name)
        outs_c_squash = os.path.join(sub_c_squash, table_name)
        outs_w_squash = os.path.join(sub_w_squash, table_name)
    else:
        outputs_csv = args.outfile_counts
        pref, suf = outputs_csv.split('.')
        outs = outputs_csv

    aggr_c_squash, aggr_w_squash = ([] for _ in range(2))

    passed, w2_pass, w2_total = 0, 0., 0.
    
    if args.aggr:
        # Unweighted counts
        with open(outs_c, 'w') as fcsv1:
            for il,l in enumerate(aggr_outs):
                split_line = [x.replace('\n', '') for x in l.split(sep)]
                if all( not x for x in split_line ):
                        continue

                # write a new line only once per numerator/denominator pair of lines
                dataset, atype, comb, ref, c, eff = split_line
                if atype != 'Type':
                    if atype=='Numerator':
                        passed = ROOT.TH1I('h_pass'+str(il), 'h_pass'+str(il), 1, 0., 1.)
                        for _ in range(int(c)):
                            passed.Fill(0.5)
                        continue
                    elif atype=='Denominator':
                        total = ROOT.TH1I('h_pass'+str(il), 'h_pass'+str(il), 1, 0., 1.)
                        for _ in range(int(c)):
                            total.Fill(0.5)
                    else:
                        if (atype != 'Numerator_weighted' and atype != 'Denominator_weighted' and
                            atype != 'Numerator_w2' and atype != 'Denominator_w2'):
                            mes = 'Type {} does not exist.'
                            raise ValueError(mes.format(atype))

                    # only lines with "Denominator" reach the following
                    if not ROOT.TEfficiency.CheckConsistency(passed,total):
                        raise ValueError('Bad histogram for TEfficiency')
                    eff = ROOT.TEfficiency(passed, total)
                    efflow = str(round(eff.GetEfficiencyErrorLow(1),3))
                    effup  = str(round(eff.GetEfficiencyErrorUp(1),3))
                    effval = (str(round(eff.GetEfficiency(1),3)) +
                              ' +' + effup + ' -' + efflow)
                    newline = sep.join((dataset, ref, comb,
                                        str(passed.GetBinContent(1)), c, effval))
                    fcsv1.write(newline + '\n')
                    aggr_c_squash.append(newline)

                if atype=='Type' and il==0:
                    newline = sep.join(('File Type', 'Reference', 'Intersection',
                                        'Pass', 'Total', 'Efficiency'))
                    fcsv1.write(newline + '\n')

        pass_vals, tot_vals = ({} for _ in range(2))
        eff_up_vals, eff_down_vals = ({} for _ in range(2))
            
        with open(outs_c_squash, 'w') as fcsv1_squash:
            # sum values from different samples for the same trigger intersection ("squash")
            for il,l in enumerate(aggr_c_squash):
                split_line = [x.replace('\n', '') for x in l.split(sep)]
                if all( not x for x in split_line ):
                    continue

                _, ref, comb, npass, ntot, eff = split_line
                pass_vals.setdefault((comb,ref), 0.)
                tot_vals.setdefault((comb,ref), 0.)
                eff_up_vals.setdefault((comb,ref), 0.)
                eff_down_vals.setdefault((comb,ref), 0.)

                pass_vals[(comb,ref)] += float(npass)
                tot_vals[(comb,ref)] += float(ntot)

                eff_up = float(re.findall('.+\+(.+)\s.+', eff)[0])
                eff_down = float(re.findall('.+-(.+)$', eff)[0])
                eff_up_vals[(comb,ref)] += eff_up*eff_up
                eff_down_vals[(comb,ref)] += eff_down*eff_down

            for k in eff_up_vals:
                eff_up_vals[k] = np.sqrt(eff_up_vals[k])
            for k in eff_down_vals:
                eff_down_vals[k] = np.sqrt(eff_down_vals[k])
                
            newline = sep.join(('File Type', 'Reference', 'Intersection',
                                'Pass', 'Total', 'Efficiency'))
            fcsv1_squash.write(newline + '\n')
    
            # calculate the new efficiency and write
            gzip = zip(pass_vals.items(), tot_vals.items(), eff_up_vals.items(), eff_down_vals.items())
            for (pk,pv),(tk,tv),(euk,euv),(edk,edv) in gzip:
                assert pk == tk
                assert pk == euk
                assert euk == edk
                eff = str(float(pv)/float(tv)) + '+' + str(euv) + ' -' + str(edv)
                newline = sep.join((pk[1], pk[0], str(pv), str(tv), eff))
                fcsv1_squash.write(newline + '\n')


        # Weighted counts
        with open(outs_w, 'w') as fcsv2:
            for il,l in enumerate(aggr_outs):
                split_line = [x.replace('\n', '') for x in l.split(sep)]
                if all( not x for x in split_line ):
                        continue

                # write a new line only once per numerator/denominator pair of lines
                dataset, atype, comb, ref, c, eff = split_line
                if atype != 'Type':
                    if atype=='Numerator_weighted':
                        passed = ROOT.TH1F('hw_pass'+str(il), 'hw_pass'+str(il), 1, 0., 1.)
                        for _ in range(int(c)):
                            passed.Fill(0.5)
                        continue                            
                    elif atype=='Denominator_weighted':
                        total = ROOT.TH1F('hw_pass'+str(il), 'hw_pass'+str(il), 1, 0., 1.)
                        for _ in range(int(c)):
                            total.Fill(0.5)
                        continue
                    elif atype=='Numerator_w2':
                        w2_pass = np.sqrt(c)
                        continue
                    elif atype=='Denominator_w2':
                        w2_total = np.sqrt(c)
                    else:
                        if atype != 'Numerator' and atype != 'Denominator':
                            mes = 'Type {} does not exist.'
                            raise ValueError(mes.format(atype))

                    # only lines with "Denominator_w2" reach the following
                    if not ROOT.TEfficiency.CheckConsistency(passed,total):
                        raise ValueError('Bad histogram for TEfficiency')
                    eff = ROOT.TEfficiency(passed, total)
                    efflow = str(round(eff.GetEfficiencyErrorLow(1),3))
                    effup  = str(round(eff.GetEfficiencyErrorUp(1),3))
                    effval = (str(round(eff.GetEfficiency(1),3)) +
                              ' +' + effup + ' -' + efflow)
                    pass_str = (str(round(passed.GetBinContent(1),3)) + ' +-' + str(w2_pass))
                    total_str = (str(round(total.GetBinContent(1),3)) + ' +-' + str(w2_total))
                    
                    newline = sep.join((dataset, ref, comb, pass_str, total_str, effval))
                    fcsv2.write(newline + '\n')
                    aggr_w_squash.append(newline)

                if atype=='Type' and il==0:
                    newline = sep.join(('File Type', 'Reference', 'Intersection',
                                        'Weighted Pass', 'Weighted Total', 'Efficiency'))
                    fcsv2.write(newline + '\n')

        pass_vals, tot_vals = ({} for _ in range(2))
        eff_up_vals, eff_down_vals = ({} for _ in range(2))
            
        with open(outs_w_squash, 'w') as fcsv2_squash:
            # sum values from different samples for the same trigger intersection ("squash")
            for il,l in enumerate(aggr_c_squash):
                split_line = [x.replace('\n', '') for x in l.split(sep)]
                if all( not x for x in split_line ):
                    continue

                _, ref, comb, npass, ntot, eff = split_line
                pass_vals.setdefault((comb,ref), 0.)
                tot_vals.setdefault((comb,ref), 0.)
                eff_up_vals.setdefault((comb,ref), 0.)
                eff_down_vals.setdefault((comb,ref), 0.)

                npass, npass_err = re.findall('(.+)\+-(.+)$', npass)[0]
                ntot, ntot_err = re.findall('(.+)\+-(.+)$', ntot)[0]
                pass_vals[(comb,ref)] += float(npass)
                tot_vals[(comb,ref)] += float(ntot)
                pass_err_vals[(comb,ref)] += float(npass_err)*float(npass_err)
                tot_err_vals[(comb,ref)] += float(ntot_err)*float(ntot_err)
                
                eff_up = float(re.findall('.+\+(.+)\s.+', eff)[0])
                eff_down = float(re.findall('.+-(.+)$', eff)[0])
                eff_up_vals[(comb,ref)] += eff_up*eff_up
                eff_down_vals[(comb,ref)] += eff_down*eff_down

            for i in range(len(pass_err_vals)):
                pass_err_vals[i] = np.sqrt(pass_err_vals[i])                
            for i in range(len(tot_err_vals)):
                tot_err_vals[i] = np.sqrt(tot_err_vals[i])                
            for i in range(len(eff_up_vals)):
                eff_up_vals[i] = np.sqrt(eff_up_vals[i])
            for i in range(len(eff_down_vals)):
                eff_down_vals[i] = np.sqrt(eff_down_vals[i])
                
            newline = sep.join(('File Type', 'Reference', 'Intersection',
                                'Pass', 'Total', 'Efficiency'))
            fcsv2_squash.write(newline + '\n')
    
            # calculate the new efficiency and write
            gzip = zip(pass_vals.items(), tot_vals.items(),
                       pass_err_vals.items(), tot_err_vals.items(),
                       eff_up_vals.items(), eff_down_vals.items())
            for (pk,pv),(pek,pev),(tek,tev),(tk,tv),(euk,euv),(edk,edv) in gzip:
                assert pk == tk
                assert pk == euk
                assert euk == edk
                pass_v = str(pv) + '+-' + str(pev)
                tot_v = str(tv) + '+-' + str(tev) 
                eff = str(float(pv)/float(tv)) + '+' + float(euv) + ' -' + float(edv)
                newline = sep.join((pk[1], pk[0], pass_v, tot_v, eff))
                fcsv2_squash.write(newline + '\n')


    else: # paired with 'if args.aggr'
        with open(outs, 'w') as fcsv:
            c_ref_combs, c_ref_vals   = ([] for _ in range(2))
            c_int_combs, c_int_vals   = ([] for _ in range(2))
            w_ref_combs, w_ref_vals   = ([] for _ in range(2))
            w_int_combs, w_int_vals   = ([] for _ in range(2))
            w2_ref_combs, w2_ref_vals = ([] for _ in range(2))
            w2_int_combs, w2_int_vals = ([] for _ in range(2))
            c_refs, w_refs, w2_refs = ([] for _ in range(3))
            for comb, val in c_inters.items():
                c_ref_combs.append(comb)
                c_ref_vals.append(c_ref[comb])
                c_int_combs.append(comb)
                c_int_vals.append(val)
                c_refs.append(reftrigs[comb])
            for comb, val in w_inters.items():
                w_ref_combs.append(comb)
                w_ref_vals.append(w_ref[comb])
                w_int_combs.append(comb)
                w_int_vals.append(val)
                w_refs.append(reftrigs[comb])
            for comb, val in w2_inters.items():
                w2_ref_combs.append(comb)
                w2_ref_vals.append(w2_ref[comb])
                w2_int_combs.append(comb)
                w2_int_vals.append(val)
                w2_refs.append(reftrigs[comb])
     
            c_ref_combs = np.array(c_ref_combs)
            c_ref_vals  = np.array(c_ref_vals)
            c_int_combs = np.array(c_int_combs)
            c_int_vals  = np.array(c_int_vals)
            c_refs      = np.array(c_refs)
            w_ref_combs = np.array(w_ref_combs)
            w_ref_vals  = np.array(w_ref_vals)
            w_int_combs = np.array(w_int_combs)
            w_int_vals  = np.array(w_int_vals)
            w_refs      = np.array(w_refs)
            w2_ref_combs = np.array(w2_ref_combs)
            w2_ref_vals  = np.array(w2_ref_vals)
            w2_int_combs = np.array(w2_int_combs)
            w2_int_vals  = np.array(w2_int_vals)
            w2_refs      = np.array(w2_refs)
     
            #remove zeros
            c_zeromask   = c_int_vals == 0
            c_ref_vals   = c_ref_vals[~c_zeromask]
            c_ref_combs  = c_ref_combs[~c_zeromask]
            c_int_vals   = c_int_vals[~c_zeromask]
            c_int_combs  = c_int_combs[~c_zeromask]
            c_refs       = c_refs[~c_zeromask]

            w_ref_vals   = w_ref_vals[~c_zeromask]
            w_ref_combs  = w_ref_combs[~c_zeromask]
            w_int_vals   = w_int_vals[~c_zeromask]
            w_int_combs  = w_int_combs[~c_zeromask]
            w_refs       = w_refs[~c_zeromask]

            w2_ref_vals  = w2_ref_vals[~c_zeromask]
            w2_ref_combs = w2_ref_combs[~c_zeromask]
            w2_int_vals  = w2_int_vals[~c_zeromask]
            w2_int_combs = w2_int_combs[~c_zeromask]
            w2_refs      = w2_refs[~c_zeromask]

            line = sep.join(('Type', 'Trigger Intersection', 'Reference', 'Counts', 'Efficiency\n'))
            fcsv.write(line)

            # calculate efficiencies
            c_effs, w_effs, w2_effs = ([] for _ in range(3))

            c_gzip = zip(c_ref_vals, c_ref_combs, c_int_vals, c_int_combs, c_refs)
            for refv, refc, intv, intc, refref in c_gzip:
                c_effs.append(float(intv) / float(refv))
            c_effs = np.sort(np.array(c_effs))[::-1]

            w_gzip = zip(w_ref_vals, w_ref_combs, w_int_vals, w_int_combs, w_refs)
            for refv, refc, intv, intc, refref in w_gzip:
                w_effs.append(float(intv) / float(refv))
            w_effs = np.sort(np.array(w_effs))[::-1]

            # sort other columns following unweighted efficiencies descending order
            c_idxs_sort  = c_effs.argsort(axis=None)[::-1]
            c_ref_vals   = c_ref_vals[c_idxs_sort]
            c_ref_combs  = c_ref_combs[c_idxs_sort]
            c_int_vals   = c_int_vals[c_idxs_sort]
            c_int_combs  = c_int_combs[c_idxs_sort]
            
            w_ref_vals   = w_ref_vals[c_idxs_sort]
            w_ref_combs  = w_ref_combs[c_idxs_sort]
            w_int_vals   = w_int_vals[c_idxs_sort]
            w_int_combs  = w_int_combs[c_idxs_sort]
            
            w2_ref_vals  = w2_ref_vals[c_idxs_sort]
            w2_ref_combs = w2_ref_combs[c_idxs_sort]
            w2_int_vals  = w2_int_vals[c_idxs_sort]
            w2_int_combs = w2_int_combs[c_idxs_sort]

            c_gzip = zip(c_effs, c_ref_vals, c_ref_combs, c_int_vals, c_int_combs, c_refs)
            for eff, refv, refc, intv, intc, refref in c_gzip:
                newline = sep.join(('Numerator',
                                    str(intc).replace(main.inters_str, '  AND  '),
                                    refref, str(intv), str(round(eff,4)))) + '\n'
                fcsv.write(newline)
     
                newline = sep.join(('Denominator',
                                    str(refc).replace(main.inters_str, '  AND  '),
                                    refref, str(refv), '1\n' ))
                fcsv.write(newline)
                
            w_gzip = zip(w_effs, w_ref_vals, w_ref_combs, w_int_vals, w_int_combs, w_refs)
            for eff, refv, refc, intv, intc, refref in w_gzip:
                newline = sep.join(('Numerator_weighted',
                                    str(intc).replace(main.inters_str, '  AND  '),
                                    refref, str(intv), str(round(eff,4)))) + '\n'
                fcsv.write(newline)
     
                newline = sep.join(('Denominator_weighted',
                                    str(refc).replace(main.inters_str, '  AND  '),
                                    refref, str(refv), '1\n' ))
                fcsv.write(newline)

    print('Save file {}.'.format(outs))
        
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Command line parser')

    parser.add_argument('--indir', dest='indir', required=True, help='SKIM directory')
    parser.add_argument('--outdir', dest='outdir', required=True, help='output directory')
    parser.add_argument('--subtag', dest='subtag', required=True,
                        help='Additional (sub)tag to differ  entiate similar runs within the same tag.')
    parser.add_argument('--tprefix', dest='tprefix', required=True, help='Targets name prefix.')
    parser.add_argument('--aggregation_step', dest='aggr', required=True, type=int,
                        help='Whether to run the sample aggregation step or the "per sample step"')
    parser.add_argument('--dataset_name', dest='dataset_name', required=True,
                        help='Name of the dataset being used.')

    parser.add_argument('--sample', dest='sample', required=False,
                        help='Process name as in SKIM directory. Used for the first step only.')
    
    parser.add_argument('--infile_counts', dest='infile_counts', required=False, nargs='+', type=str,
                        help='Name of input csv files with counts. Used for the aggregation step only.')
    parser.add_argument('--outfile_counts', dest='outfile_counts',
                        help='Name of output csv files with counts.')
    parser.add_argument('--channel', dest='channel', required=False,
                        help='Channel to be used for the aggregation.')
    args = utils.parse_args(parser)
    
    add_trigger_counts(args)
