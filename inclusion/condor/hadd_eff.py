# coding: utf-8

_all_ = [ 'hadd_eff', 'hadd_eff_outputs' ]

import os
import sys
parent_dir = os.path.abspath(__file__ + 3 * '/..')
sys.path.insert(0, parent_dir)

import inclusion
from inclusion import config
from inclusion.utils import utils
from inclusion.condor.job_writer import JobWriter

def _subpaths(args):
    _tbase1 = args.outprefix
    _tbase2 = '_Sum' + args.subtag
    return _tbase1, _tbase2

@utils.set_pure_input_namespace
def run_hadd_eff_outputs(args):
    targets = []

    # add the merge of all the samples first
    _tbase1, _tbase2 = _subpaths(args)
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
def hadd_eff_outputs(args):
    """
    Outputs are guaranteed to have the same length.
    Returns all separate paths to avoid code duplication.
    """
    return JobWriter.define_output( localdir=args.localdir,
                                    data_folders='HaddEff',
                                    tag=args.tag )

@utils.set_pure_input_namespace
def hadd_eff(args):
    """Adds ROOT histograms"""
    targets = run_hadd_eff_outputs(args)
    outs_job, outs_submit, outs_check, outs_log = hadd_eff_outputs(args)
    jw = JobWriter()
    comm = 'hadd -f ${1} ${@:2}'

    for out in outs_job:
        jw.write_shell(filename=out, command=comm, localdir=args.localdir)
        if out == outs_job[0]:
            jw.add_string('echo "HaddEff done."')
        elif out == outs_job[1]:
            jw.add_string('echo "HaddEff Agg done."')

    #### Write submission file
    inputs_join = []
    for out1,out2,out3,out4 in zip(outs_job,outs_submit,outs_check,outs_log):
        jw.write_condor(filename=out2,
                        real_exec='/dev/null',
                        shell_exec=out1,
                        outfile=out3,
                        logfile=out4,
                        queue=config.queue,
                        machine='llrt3condor')

        qlines = []
        if out1 == outs_job[0]:
            for t,smpl in zip(targets[1:], args.samples):
                inputs = os.path.join(args.indir, smpl, args.outprefix + '*' + args.subtag + '.root ')
                inputs_join.append(t)
                # join subdatasets (different MC or Data subfolders, ex: TT_fullyHad, TT_semiLep, ...)
                jw.add_string('  {}, {}'.format(t, inputs))
        elif out1 == outs_job[1]:
            # join MC or Data subdatasets into a single one (ex: TT)
            jw.add_string('  {}, {}'.format(targets[0], ' '.join(inputs_join)))

        jw.write_queue( qvars=('myoutput', 'myinputs'),
                        qlines=qlines )
