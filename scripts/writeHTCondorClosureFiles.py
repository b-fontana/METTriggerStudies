import sys
sys.path.append("..")

import os
import argparse

from utils.utils import (
    build_prog_path,
    setPureInputNamespace,
)

from scripts.jobWriter import JobWriter

@setPureInputNamespace
def writeHTCondorClosureFiles_outputs(args):
    job_f, subm_f, check_f = JobWriter.define_output( localdir=args.localdir,
                                                      data_folders='Closure',
                                                      tag=args.tag )
    return job_f[0], subm_f[0], check_f[0]

@setPureInputNamespace
def writeHTCondorClosureFiles(args):
    outs_job, outs_submit, outs_check = writeHTCondorClosureFiles_outputs(args)
    jw = JobWriter()

    #### Write shell executable (python scripts must be wrapped in shell files to run on HTCondor)
    prog = build_prog_path(args.localdir, 'runClosure.py')
    command =  ( '{prog} --indir_eff {inref} '.format(prog=prog, inref=args.indir_eff)
                 + '--indir_union {inunion} '.format(inunion=args.indir_union)
                 + '--indir_json {injson} '.format(injson=args.indir_json)
                 + '--mc_processes {procs} '.format(procs=' '.join(args.mc_processes))
                 + '--outdir {outdir} '.format(outdir=args.outdir)
                 + '--in_prefix {inprefix} '.format(inprefix=args.inprefix)
                 + '--channel ${1} '
                 + '--closure_single_trigger ${2} '
                 + '--variables {variables} '.format(variables=' '.join(args.variables))
                 + '--subtag {subtag} '.format(subtag=args.subtag)
                 + '--binedges_fname {be} '.format(be=args.binedges_filename)
                 + '--data_name {dn} '.format(dn=args.data_name)
                 + '--mc_name {mn} '.format(mn=args.mc_name)
                 + '--eff_prefix {effprefix} '.format(effprefix=args.eff_prefix)
                )
    
    if args.debug:
        command += '--debug '

    jw.write_init(outs_job, command, args.localdir)
    jw.add_string('echo "runClosure for channel ${1} and single trigger ${2} done."')

    #### Write submission file
    jw.write_init( filename=outs_submit,
                   executable=outs_job,
                   outfile=outs_check,
                   queue='short' )

    qlines = []
    for chn in args.channels:
        for trig in args.closure_single_triggers:
            qlines.append('  {},{}'.format(chn,trig))

    jw.write_queue( qvars=('channel', 'closure_single_trigger'),
                    qlines=qlines )
