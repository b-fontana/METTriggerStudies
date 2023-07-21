# coding: utf-8

_all_ = [ 'discriminator', 'discriminator_outputs' ]

import os
import sys
parent_dir = os.path.abspath(__file__ + 3 * '/..')
sys.path.insert(0, parent_dir)

import argparse

import inclusion
from inclusion.config import main
from inclusion.utils import utils
from inclusion.condor.job_writer import JobWriter

@utils.set_pure_input_namespace
def discriminator_outputs(args):
    """
    One output per channel. Allows channel parallellization with DAGMAN.
    """
    return JobWriter.define_output( localdir=args.localdir,
                                    data_folders=['Discriminator_' + x for x in args.channels],
                                    tag=args.tag )

@utils.set_pure_input_namespace
def discriminator(args):
    script = 'run_var_discriminator.py'

    outs_job, outs_submit, outs_check, outs_log = discriminator_outputs(args)
    jw = JobWriter()
    #### Write shell executable (python scripts must be wrapped in shell files to run on HTCondor)
    for i,chn in enumerate(args.channels):
        pars = {'indir'     : args.indir,
                'outdir'    : args.outdir,
                'channel'   : chn,
                'variables' : ' '.join(args.variables,),
                'tag'       : args.tag,
                'subtag'    : args.subtag,
                'configuration' : args.configuration}
        comm = utils.build_script_command(name=script, sep=' ', **pars)
        if args.debug:
            comm += '--debug '
        
        jw.write_shell(filename=outs_job[i], command=comm, localdir=args.localdir)
        jw.add_string('echo "Script {} with channel {} done."'.format(script, args.channels[i]))

        #### Write submission file
        jw.write_condor(filename=outs_submit[i],
                        real_exec=utils.build_script_path(script),
                        shell_exec=outs_job[i],
                        outfile=outs_check[i],
                        logfile=outs_log[i],
                        queue=main.queue,
                        machine='llrt3condor')
        jw.write_queue()

# -- Parse options
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Command line parser')

    parser.add_argument('--localdir', dest='localdir', default=os.getcwd(),
                        help='out directory')
    parser.add_argument('--indir', dest='indir', required=True, help='in directory')
    parser.add_argument('--outdir', dest='outdir', required=True, help='out directory')
    parser.add_argument('--tag', dest='tag', required=True, help='tag')
    parser.add_argument('--subtag', dest='subtag', required=True, help='subtag')
    parser.add_argument('--data_name', dest='data_name', required=True, help='Data sample name')
    parser.add_argument('--mc_name', dest='mc_name', required=True, help='MC sample name')
    parser.add_argument('--channels',   dest='channels', required=True, nargs='+', type=str,
                        help='Select the channels over which the workflow will be run.' )
    parser.add_argument('--triggers', dest='triggers', required=True, nargs='+', type=str,
                        help='Select the triggers over which the workflow will be run.' )
    parser.add_argument('--variables', dest='variables', required=True, nargs='+', type=str,
                        help='Select the variables over which the workflow will be run.' )
    parser.add_argument('--configuration', dest='configuration', required=True,
                        help='Name of the configuration module to use.')
    parser.add_argument('--debug', action='store_true', help='debug verbosity')
    args = parser.parse_args()

    submitTriggerEff( args )
