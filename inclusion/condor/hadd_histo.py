# coding: utf-8

_all_ = [ 'hadd_histo', 'hadd_histo_outputs' ]

import os
import sys
from utils import utils
from condor.job_writer import JobWriter

@utils.set_pure_input_namespace
def run_hadd_histo_outputs(args):
    targets = []

    # add the merge of all the samples first
    _tbase1, _tbase2 = utils.hadd_subpaths(args)
    tbase = _tbase1 + _tbase2
    t = os.path.join( args.indir, tbase + '.root' )
    targets.append( t )

    # add individual sample merges
    for smpl in args.samples:
        tbase = _tbase1 + '_' + smpl + _tbase2
        t = os.path.join( args.indir, tbase + '.root' )
        targets.append( t )
        
    return targets

@utils.set_pure_input_namespace
def hadd_histo_outputs(args):
    """
    Outputs are guaranteed to have the same length.
    Returns all separate paths to avoid code duplication.
    """
    ret = JobWriter.define_output( localdir=args.localdir,
                                   data_folders=['HaddHisto' + args.dataset_name,
                                                 'HaddHistoAgg' + args.dataset_name],
                                   tag=args.tag )
    return ret


@utils.set_pure_input_namespace
def hadd_histo(args):
    """Adds ROOT histograms"""
    script = os.path.basename(__file__)
    targets = run_hadd_histo_outputs(args)
    outs_job, outs_submit, outs_check, outs_log = hadd_histo_outputs(args)
    jw = JobWriter()
    comm = 'hadd -f ${1} ${@:2}'

    for out in outs_job:
        jw.write_shell(filename=out, command=comm, localdir=args.localdir)
        if out == outs_job[0]:
            jw.add_string('echo "{} without aggregation (dataset {}) done."'.format(script, args.dataset_name))
        elif out == outs_job[1]:
            jw.add_string('echo "{} with aggregation (dataset {}) done."'.format(script, args.dataset_name))

    #### Write submission file
    inputs_join = []
    for out1,out2,out3,out4 in zip(outs_job,outs_submit,outs_check,outs_log):
        jw.write_condor( filename=out2,
                         executable=out1,
                         outfile=out3,
                         logfile=out4,
                         queue='long',
                         machine='llrt3condor' )

        qlines = []
        if out1 == outs_job[0]:
            for t,smpl in zip(targets[1:], args.samples):
                inputs = os.path.join(args.indir, smpl, args.tprefix + '*' + args.subtag + '.root')
                inputs_join.append(t)
                # join subdatasets (different MC or Data subfolders, ex: TT_fullyHad, TT_semiLep, ...)
                qlines.append('  {}, {}'.format(t, inputs))
        elif out1 == outs_job[1]:
            # join MC or Data subdatasets into a single one (ex: TT)
            qlines.append('  {}, {}'.format(targets[0], ' '.join(inputs_join)))
        
        jw.write_queue( qvars=('myoutput', 'myinputs'),
                        qlines=qlines )
