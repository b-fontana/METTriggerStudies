# coding: utf-8

_all_ = [ "eff_and_sf", "eff_and_sf_outputs" ]

import sys
sys.path.append("..")

import os
import argparse

from utils.utils import (
    build_prog_path,
    generate_trigger_combinations,
    join_name_trigger_intersection as joinNTC,
    join_strings,
    set_pure_input_namespace,
)
from condor.job_writer import JobWriter

@set_pure_input_namespace
def eff_and_sf_outputs(args):
    """
    Outputs are guaranteed to have the same length.
    Returns all separate paths to avoid code duplication.
    """
    job_f, subm_f, check_f = JobWriter.define_output( localdir=args.localdir,
                                                      data_folders='EffAndScaleFactors',
                                                      tag=args.tag )
    return job_f[0], subm_f[0], check_f[0]

@set_pure_input_namespace
def eff_and_sf(args):
    prog = build_prog_path(args.localdir, 'runEfficienciesAndScaleFactors.py')
    outs_job, outs_submit, outs_check = eff_and_sf_outputs(args)
    jw = JobWriter()

    #### Write shell executable (python scripts must be wrapped in shell files to run on HTCondor)
    command = join_strings( '{} '                .format(prog),
                            '--indir {} '        .format(args.indir),
                            '--outdir {} '       .format(args.outdir),
                            '--mc_keys {} '      .format(' '.join(args.mc_keys)),
                            '--mc_vals {} '      .format(' '.join(args.mc_vals)),
                            '--data_keys {} '    .format(' '.join(args.data_keys)),
                            '--data_vals {} '    .format(' '.join(args.data_vals)),
                            '--triggercomb ${1} ',
                            '--channels {} '     .format(' '.join(args.channels)),
                            '--variables {} '    .format(' '.join(args.variables)),
                            '--subtag {} '       .format(args.subtag),
                            '--tprefix {} '      .format(args.tprefix),
                            '--canvas_prefix {} '.format(args.canvas_prefix),
                        )

    if args.draw_independent_MCs:
        command += '--draw_independent_MCs '
    if args.debug:
        command += '--debug '

    jw.write_shell(outs_job, command=command, localdir=args.localdir)
    jw.add_string('echo "runEfficienciesAndScaleFactors done."')

    #### Write submission file
    jw.write_condor( filename=outs_submit,
                     executable=outs_job,
                     outfile=outs_check,
                     queue='short',
                     machine='llrt3condor' )

    qlines = []
    for chn in args.channels:
        if chn == args.channels[0]:
            triggercomb = generate_trigger_combinations(chn, args.triggers)
        else:
            triggercomb += generate_trigger_combinations(chn, args.triggers)
            
    for tcomb in triggercomb:
        qlines.append('  {}'.format(joinNTC(tcomb)))

    jw.write_queue( qvars=('triggercomb',),
                    qlines=qlines )
