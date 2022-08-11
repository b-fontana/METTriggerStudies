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
    outs_job, outs_submit, outs_check, outs_log = eff_and_sf_aggr_outputs(args)

    #### Write shell executable (python scripts must be wrapped in shell files to run on HTCondor)
    pars = {'--indir'       : args.indir,
            '--outdir'      : args.outdir,
            '--channel'     : '${1}',
            '--file_prefix' : args.file_prefix,
            '--variables'   : ' '.join(args.variables,),
            '--subtag'      : args.subtag}

    script = 'aggr_eff_and_sf.py'
    comm = utils.build_script_command(name=script, sep=' ', **pars)
    if args.debug:
        comm += '--debug '

    jw = JobWriter()
    jw.write_shell(filename=outs_job, command=comm, localdir=args.localdir)
    jw.add_string('echo "{} done."'.format(script))

    #### Write submission file
    jw.write_condor(filename=outs_submit,
                    real_exec=utils.build_script_path(script),
                    shell_exec=outs_job,
                    outfile=outs_check,
                    logfile=outs_log,
                    queue='long',
                    machine='llrt3condor')

    qlines = []
    for chn in args.channels:
        qlines.append('  {}'.format(chn))

    jw.write_queue( qvars=('channel',),
                    qlines=qlines )
