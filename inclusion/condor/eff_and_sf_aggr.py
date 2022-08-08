# coding: utf-8

_all_ = [ 'eff_and_sf_aggr', 'eff_and_sf_aggr_outputs' ]

import sys
sys.path.append("..")

import os
import argparse

from utils import utils
from condor.job_writer import JobWriter

@utils.set_pure_input_namespace
def eff_and_sf_aggr_outputs(args):
    job_f, subm_f, check_f, log_f = JobWriter.define_output( localdir=args.localdir,
                                                             data_folders='EffAndSFAgg',
                                                             tag=args.tag )
    return job_f[0], subm_f[0], check_f[0], log_f[0]

@utils.set_pure_input_namespace
def eff_and_sf_aggr(args):
    script = 'aggr_eff_and_sf.py'
    prog = utils.build_prog_path(args.localdir, script)
    outs_job, outs_submit, outs_check, outs_log = eff_and_sf_aggr_outputs(args)
    jw = JobWriter()

    #### Write shell executable (python scripts must be wrapped in shell files to run on HTCondor)
    command = utils.join_strings('{}'.format(prog),
                                 '--indir {}'.format(args.indir),
                                 '--outdir {}'.format(args.outdir),
                                 '--channel ${1}',
                                 '--file_prefix {}'.format(args.file_prefix),
                                 '--variables {}'.format(' '.join(args.variables,)),
                                 '--subtag {}'.format(args.subtag),
                                 sep=' ')

    if args.debug:
        command += '--debug '

    jw.write_shell(filename=outs_job, command=command, localdir=args.localdir)
    jw.add_string('echo "{} done."'.format(script))

    #### Write submission file
    jw.write_condor( filename=outs_submit,
                     executable=outs_job,
                     outfile=outs_check,
                     logfile=outs_log,
                     queue='long',
                     machine='llrt3condor' )

    qlines = []
    for chn in args.channels:
        qlines.append('  {}'.format(chn))

    jw.write_queue( qvars=('channel',),
                    qlines=qlines )
