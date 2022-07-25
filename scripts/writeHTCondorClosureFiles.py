import sys
sys.path.append("..")

import os
import argparse

from utils.utils import (
    build_prog_path,
    join_strings,
    set_pure_input_namespace,
)

from scripts.jobWriter import JobWriter

@set_pure_input_namespace
def writeHTCondorClosureFiles_outputs(args):
    job_f, subm_f, check_f = JobWriter.define_output( localdir=args.localdir,
                                                      data_folders='Closure',
                                                      tag=args.tag )
    return job_f[0], subm_f[0], check_f[0]

@set_pure_input_namespace
def writeHTCondorClosureFiles(args):
    outs_job, outs_submit, outs_check = writeHTCondorClosureFiles_outputs(args)
    jw = JobWriter()

    #### Write shell executable (python scripts must be wrapped in shell files to run on HTCondor)
    prog = build_prog_path(args.localdir, 'runClosure.py')
    command = join_strings( '{} '.format(prog),
                            '--indir_eff {} '.format(args.indir_eff),
                            '--indir_union {} '.format(args.indir_union),
                            '--indir_json {} '.format(args.indir_json),
                            '--mc_processes {} '.format(' '.join(args.mc_processes)),
                            '--outdir {} '.format(args.outdir),
                            '--in_prefix {} '.format(args.inprefix),
                            '--channel ${1} ',
                            '--closure_single_trigger ${2} ',
                            '--variables {} '.format(' '.join(args.variables)),
                            '--subtag {} '.format(args.subtag),
                            '--binedges_fname {} '.format(args.binedges_filename),
                            '--data_name {} '.format(args.data_name),
                            '--mc_name {} '.format(args.mc_name),
                            '--eff_prefix {} '.format(args.eff_prefix) )
    
    if args.debug:
        command += '--debug '

    jw.write_shell(filename=outs_job, command=command, localdir=args.localdir)
    jw.add_string('echo "runClosure for channel ${1} and single trigger ${2} done."')

    #### Write submission file
    jw.write_condor( filename=outs_submit,
                     executable=outs_job,
                     outfile=outs_check,
                     queue='short',
                     machine='llrt3condor7' )

    qlines = []
    for chn in args.channels:
        for trig in args.closure_single_triggers:
            qlines.append('  {},{}'.format(chn,trig))

    jw.write_queue( qvars=('channel', 'closure_single_trigger'),
                    qlines=qlines )
