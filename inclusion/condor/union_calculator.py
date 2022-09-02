# coding: utf-8

_all_ = [ 'union_calculator', 'union_calculator_outputs' ]

import os
import sys
parent_dir = os.path.abspath(__file__ + 3 * '/..')
sys.path.insert(0, parent_dir)

import inclusion
from inclusion.config import main
from inclusion.utils import utils
from inclusion.condor.job_writer import JobWriter

@utils.set_pure_input_namespace
def union_calculator_outputs(args):
    """
    Outputs are guaranteed to have the same length.
    Returns all separate paths to avoid code duplication.
    """
    base_name = 'UnionWeightsCalculator'
    data_folders = [ os.path.join( base_name, proc) for proc in args.mc_processes ]
    return JobWriter.define_output( localdir=args.localdir,
                                    data_folders=data_folders,
                                    tag=args.tag )

@utils.set_pure_input_namespace
def union_calculator(args):
    jobs, subs, checks, logs = union_calculator_outputs(args)
    jw = JobWriter()

    for i,proc in enumerate(args.mc_processes):
        filelist, inputdir = utils.get_root_inputs(proc, args.indir_root)

        #### Write shell executable (python scripts must be wrapped in shell files to run on HTCondor)
        pars = {'indir_root'             : inputdir,
                'indir_json'             : args.indir_json,
                'indir_eff'              : args.indir_eff,
                'outdir'                 : args.outdir,
                'outprefix'              : args.outprefix,
                'sample'                 : proc,
                'channels'               : ' '.join(args.channels,),
                'triggers'               : ' '.join(args.triggers,),
                'file_name'              : '${1}',
                'closure_single_trigger' : '${2}',
                'variables'              : ' '.join(args.variables,),
                'tag'                    : args.tag,
                'subtag'                 : args.subtag,
                'data_name'              : args.data_name,
                'mc_name'                : args.mc_name,
                'binedges_fname'         : args.binedges_filename,
                'configuration'          : args.configuration}

        script = 'run_union_calculator.py'
        comm = utils.build_script_command(name=script, sep=' ', **pars)
        if args.debug:
            comm += '--debug '

        jw.write_shell(filename=jobs[i], command=comm, localdir=args.localdir)
        jw.add_string('echo "Process {} done."'.format(proc))

        #### Write submission file
        jw.write_condor(filename=subs[i],
                        real_exec=utils.build_script_path(script),
                        shell_exec=jobs[i],
                        outfile=checks[i],
                        logfile=logs[i],
                        queue=main.queue,
                        machine='llrt3condor')

        qlines = []
        for listname in filelist:
            for trig in args.closure_single_triggers:
                qlines.append('  {},{}'.format( os.path.basename(listname).replace('\n',''), trig ))
                        
        jw.write_queue( qvars=('filename', 'closure_single_trigger'),
                        qlines=qlines )
