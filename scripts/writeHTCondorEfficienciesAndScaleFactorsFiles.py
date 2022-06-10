
import sys
sys.path.append("..")

import os
import argparse

from utils.utils import (
    build_prog_path,
    generate_trigger_combinations,
    join_name_trigger_intersection as joinNTC,
    set_pure_input_namespace,
)
from scripts.jobWriter import JobWriter

@set_pure_input_namespace
def writeHTCondorEfficienciesAndScaleFactorsFiles_outputs(args):
    """
    Outputs are guaranteed to have the same length.
    Returns all separate paths to avoid code duplication.
    """
    job_f, subm_f, check_f = JobWriter.define_output( localdir=args.localdir,
                                                      data_folders='EffAndScaleFactors',
                                                      tag=args.tag,
                                                      names='EfficienciesAndSF' )
    return job_f[0], subm_f[0], check_f[0]

@set_pure_input_namespace
def writeHTCondorEfficienciesAndScaleFactorsFiles(args):
    prog = build_prog_path(args.localdir, 'runEfficienciesAndScaleFactors.py')
    outs_job, outs_submit, outs_check = writeHTCondorEfficienciesAndScaleFactorsFiles_outputs(args)
    jw = JobWriter()

    #### Write shell executable (python scripts must be wrapped in shell files to run on HTCondor)
    command =  ( ( '{prog} --indir {indir} --outdir {outdir} '
                   '--mc_processes {mc_processes} '
                   '--mc_name {mc_name} --data_name {data_name} '
                   '--triggercomb ${{1}} '
                   '--channels {channels} --variables {variables} '
                   '--subtag {subtag} '
                   '--tprefix {tprefix} '
                   '--canvas_prefix {cprefix} '
                  ).format( prog=prog, indir=args.indir, outdir=args.outdir,
                            mc_processes=' '.join(args.mc_processes,),
                            mc_name=args.mc_name, data_name=args.data_name,
                            channels=' '.join(args.channels,), variables=' '.join(args.variables,),
                            subtag=args.subtag,
                            draw_independent_MCs=1 if args.draw_independent_MCs else 0,
                            tprefix=args.tprefix,
                            cprefix=args.canvas_prefix,
                           )
                )

    if args.draw_independent_MCs:
        command += '--draw_independent_MCs '
    if args.debug:
        command += '--debug '

    jw.write_init(outs_job, command, args.localdir)
    jw.add_string('echo "runEfficienciesAndScaleFactors done."')

    #### Write submission file
    jw.write_init( filename=outs_submit,
                   executable=outs_job,
                   outfile=outs_check,
                   queue='short' )

    qlines = []
    triggercomb = generate_trigger_combinations(args.triggers)
    for tcomb in triggercomb:
        qlines.append('  {}'.format(joinNTC(tcomb)))

    jw.write_queue( qvars=('triggercomb',),
                    qlines=qlines )
