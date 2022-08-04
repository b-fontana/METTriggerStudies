# coding: utf-8

_all_ = [ 'closure', 'closure_outputs' ]

import sys
sys.path.append('..')

import os
import argparse

from utils import utils
from condor.job_writer import JobWriter

@utils.set_pure_input_namespace
def closure_outputs(args):
    job_f, subm_f, check_f, log_f = JobWriter.define_output( localdir=args.localdir,
                                                             data_folders='Closure',
                                                             tag=args.tag )
    return job_f[0], subm_f[0], check_f[0], log_f[0]

@utils.set_pure_input_namespace
def closure(args):
    outs_job, outs_submit, outs_check, outs_log = closure_outputs(args)
    jw = JobWriter()

    script = 'run_closure.py'
    prog = utils.build_prog_path(args.localdir, script)
    command = utils.join_strings( '{}'.format(prog),
                                 '--indir_eff {}'.format(args.indir_eff),
                                 '--indir_union {}'.format(args.indir_union),
                                 '--indir_json {}'.format(args.indir_json),
                                 '--mc_processes {}'.format(' '.join(args.mc_processes)),
                                 '--outdir {}'.format(args.outdir),
                                 '--in_prefix {}'.format(args.inprefix),
                                 '--channel ${1}',
                                 '--closure_single_trigger ${2}',
                                 '--variables {}'.format(' '.join(args.variables)),
                                 '--subtag {}'.format(args.subtag),
                                 '--binedges_fname {}'.format(args.binedges_filename),
                                 '--data_name {}'.format(args.data_name),
                                 '--mc_name {}'.format(args.mc_name),
                                 '--eff_prefix {}'.format(args.eff_prefix),
                                 sep=' ')
    
    if args.debug:
        command += '--debug '

    jw.write_shell(filename=outs_job, command=command, localdir=args.localdir)
    jw.add_string('echo "{} for channel ${1} and single trigger ${2} done."'.format(script))

    #### Write submission file
    jw.write_condor( filename=outs_submit,
                     executable=outs_job,
                     outfile=outs_check,
                     logfile=outs_log,
                     queue='long',
                     machine='llrt3condor' )

    qlines = []
    for chn in args.channels:
        for trig in args.closure_single_triggers:
            qlines.append('  {},{}'.format(chn,trig))

    jw.write_queue( qvars=('channel', 'closure_single_trigger'),
                    qlines=qlines )
