# coding: utf-8

_all_ = [ 'closure', 'closure_outputs' ]

import os
import sys
parent_dir = os.path.abspath(__file__ + 3 * '/..')
sys.path.insert(0, parent_dir)

import inclusion
from inclusion.utils import utils
from inclusion.condor.job_writer import JobWriter

@utils.set_pure_input_namespace
def closure_outputs(args):
    job_f, subm_f, check_f, log_f = JobWriter.define_output( localdir=args.localdir,
                                                             data_folders='Closure',
                                                             tag=args.tag )
    return job_f[0], subm_f[0], check_f[0], log_f[0]

@utils.set_pure_input_namespace
def closure(args)                    :
    outs_job, outs_submit, outs_check, outs_log = closure_outputs(args)

    pars = {'indir_eff'              : args.indir_eff,
            'indir_union'            : args.indir_union,
            'indir_json'             : args.indir_json,
            'mc_processes'           : ' '.join(args.mc_processes),
            'outdir'                 : args.outdir,
            'in_prefix'              : args.inprefix,
            'channel'                : '${1}',
            'closure_single_trigger' : '${2}',
            'variables'              : ' '.join(args.variables),
            'subtag'                 : args.subtag,
            'binedges_fname'         : args.binedges_filename,
            'data_name'              : args.data_name,
            'mc_name'                : args.mc_name,
            'eff_prefix'             : args.eff_prefix }
    script = 'run_closure.py'
    comm = utils.build_script_command(name=script, sep=' ', **pars)
    if args.debug:
        comm += '--debug '

    jw = JobWriter()
    jw.write_shell(filename=outs_job, command=comm, localdir=args.localdir)
    jw.add_string('echo "{} for channel ${{1}} and single trigger ${{2}} done."'.format(script))

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
        for trig in args.closure_single_triggers:
            qlines.append('  {},{}'.format(chn,trig))

    jw.write_queue( qvars=('channel', 'closure_single_trigger'),
                    qlines=qlines )
