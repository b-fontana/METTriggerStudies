# coding: utf-8

_all_ = [ 'hadd_counts', 'hadd_counts_outputs' ]

import os
import sys
parent_dir = os.path.abspath(__file__ + 3 * '/..')
sys.path.insert(0, parent_dir)

import inclusion
from inclusion.config import main
from inclusion.utils import utils
from inclusion.condor.job_writer import JobWriter
from inclusion.utils.utils import build_script_command as bsc

@utils.set_pure_input_namespace
def run_hadd_counts_outputs(args):
    targets = []

    # add the merge of all the samples first
    _tbase1, _tbase2 = utils.hadd_subpaths(args)
    tbase = _tbase1 + _tbase2
    for chn in args.channels:
        t = os.path.join( args.indir, tbase + '_' + chn + '.csv' )
        targets.append( t )

    # add individual sample merges
    for smpl in args.samples:
        tbase = _tbase1 + '_' + smpl + _tbase2
        for chn in args.channels:
            t = os.path.join( args.indir, tbase + '_' + chn + '.csv' )
            targets.append( t )
    return targets

@utils.set_pure_input_namespace
def hadd_counts_outputs(args):
    return JobWriter.define_output( localdir=args.localdir,
                                    data_folders=['HaddCounts' + args.dataset_name,
                                                  'HaddCountsAgg' + args.dataset_name],
                                    tag=args.tag )

@utils.set_pure_input_namespace
def hadd_counts(args):
    targets = run_hadd_counts_outputs(args)
    outs_job, outs_submit, outs_check, outs_log = hadd_counts_outputs(args)

    script = 'add_trig_counts.py'
    pars = {'indir'          : args.indir,
            'outdir'         : args.outdir,
            'subtag'         : args.subtag,
            'tprefix'        : args.tprefix,
            'dataset_name'   : args.dataset_name}
    comm_base = bsc(name=script, sep=' ', **pars)

    pars1 = {'outfile_counts'   : '${1}',
             'sample'           : '${2}',
             'channel'          : '${3}',
             'aggregation_step' : '0'}
    comm1 = comm_base + bsc(name=None, sep=' ', **pars1)

    pars2 = {'channel'          : '${1}',
             'infile_counts'    : '${2}',
             'aggregation_step' : '1'}
    comm2 = comm_base + bsc(name=None, sep=' ', **pars2)
    
    #### Write shell executable (python scripts must be wrapped in shell files to run on HTCondor)
    jw = JobWriter()
    for out in outs_job:
        if out == outs_job[0]:
            jw.write_shell(filename=out, command=comm1,
                           localdir=args.localdir)
            jw.add_string('echo "{} without aggregation (dataset {}) done."'.format(script, args.dataset_name))
        elif out == outs_job[1]:
            jw.write_shell(filename=out, command=comm2,
                           localdir=args.localdir)
            jw.add_string('echo "{} with aggregation {} done."'.format(script, args.dataset_name))

    #### Write submission file
    inputs_join = {}
    nchannels = len(args.channels)
    for out1,out2,out3,out4 in zip(outs_job,outs_submit,outs_check,outs_log):
        jw.write_condor(filename=out2,
                        real_exec=utils.build_script_path(script),
                        shell_exec=out1,
                        outfile=out3,
                        logfile=out4,
                        queue=main.queue,
                        machine=main.machine)

        qvars = None
        qlines = []
        if out1 == outs_job[0]:
            qvars = ('myoutput', 'channel', 'sample')
            for it,t in enumerate(targets[nchannels:]):
                smpl = args.samples[ int(it/nchannels) ]
                chn = args.channels[ int(it%nchannels) ]
                if chn not in inputs_join:
                    inputs_join[chn] = []
                inputs_join[chn].append(t)
                qlines.append('  {}, {}, {}'.format(t,smpl,chn))

        elif out1 == outs_job[1]:
            qvars = ('channel', 'myinputs')
            for ichn,chn in enumerate(args.channels):
                qlines.append(" {}, '{}'".format(chn, ' '.join(inputs_join[chn])))

        jw.write_queue( qvars=qvars, qlines=qlines )
