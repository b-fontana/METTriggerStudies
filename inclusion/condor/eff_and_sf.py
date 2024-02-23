# coding: utf-8

_all_ = [ 'eff_and_sf', 'eff_and_sf_outputs' ]

import os
import sys
parent_dir = os.path.abspath(__file__ + 3 * '/..')
sys.path.insert(0, parent_dir)

import importlib

import inclusion
from inclusion.config import main
from inclusion.utils import utils
from inclusion.condor.job_writer import JobWriter

@utils.set_pure_input_namespace
def eff_and_sf_outputs(args):
    """
    Outputs are guaranteed to have the same length.
    Returns all separate paths to avoid code duplication.
    """
    tmp = JobWriter.define_output( localdir=args.localdir,
                                   data_folders='EffAndScaleFactors',
                                   tag=args.tag )
    job_f, subm_f, check_f, log_f = tmp
    return job_f[0], subm_f[0], check_f[0], log_f[0]

@utils.set_pure_input_namespace
def eff_and_sf(args):
    outs_job, outs_submit, outs_check, outs_log = eff_and_sf_outputs(args)

    #### Write shell executable (python scripts must be wrapped in shell files to run on HTCondor)
    pars = {'indir'         : args.indir,
            'outdir'        : args.outdir,
            'mc_name'       : args.mc_name,
            'data_name'     : args.data_name,
            'triggercomb'   : '${1}',
            'channels'      : ' '.join(args.channels),
            'variables'     : ' '.join(args.variables),
            'configuration' : args.configuration,
            'subtag'        : args.subtag,
            'tprefix'       : args.tprefix,
            'year'          : args.year,
            'canvas_prefix' : args.canvas_prefix}

    script = 'run_eff_and_sf.py'
    comm = utils.build_script_command(name=script, sep=' ', **pars)
    if args.draw_independent_MCs:
        comm += '--draw_independent_MCs '
    if args.debug:
        comm += '--debug '

    jw = JobWriter()
    jw.write_shell(outs_job, command=comm, localdir=args.localdir)
    jw.add_string('echo "{} done."'.format(script))

    #### Write submission file
    jw.write_condor(filename=outs_submit,
                    real_exec=utils.build_script_path(script),
                    shell_exec=outs_job,
                    outfile=outs_check,
                    logfile=outs_log,
                    queue=main.queue,
                    machine=main.machine)

    cfg = importlib.import_module(args.configuration)
    qlines = []
    for chn in args.channels:
        if chn == args.channels[0]:
            triggercomb = utils.generate_trigger_combinations(chn, cfg.triggers,
                                                              cfg.exclusive)
        else:
            triggercomb += utils.generate_trigger_combinations(chn, cfg.triggers,
                                                               cfg.exclusive)
            
    for tcomb in set(triggercomb):
        qlines.append('  {}'.format( utils.join_name_trigger_intersection(tcomb)) )

    jw.write_queue( qvars=('triggercomb',),
                    qlines=qlines )
