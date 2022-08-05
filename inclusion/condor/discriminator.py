# coding: utf-8

_all_ = [ "discriminator", "discriminator_outputs" ]

import sys
sys.path.append("..")

import os
import argparse
import ROOT

from utils import utils
from condor.job_writer import JobWriter

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
    prog = utils.build_prog_path(args.localdir, script)
    outs_job, outs_submit, outs_check, outs_log = discriminator_outputs(args)
    jw = JobWriter()

    #### Write shell executable (python scripts must be wrapped in shell files to run on HTCondor)
    for i,chn in enumerate(args.channels):
        command = utils.join_strings( '{}'.format(prog),
                                      '--indir {}'.format(args.indir),
                                      '--outdir {}'.format(args.outdir),
                                      '--channel {}'.format(chn),
                                      '--triggers {}'.format(' '.join(args.triggers,)),
                                      '--variables {}'.format(' '.join(args.variables,)),
                                      '--tag {}'.format(args.tag),
                                      '--subtag {}'.format(args.subtag),
                                      sep=' ')

        if args.debug:
            command += '--debug '

        jw.write_shell(filename=outs_job[i], command=command, localdir=args.localdir)
        jw.add_string('echo "Script {} with channel {} done."'.format(script, args.channels[i]))

        #### Write submission file
        jw.write_condor( filename=outs_submit[i],
                         executable=outs_job[i],
                         outfile=outs_check[i],
                         logfile=outs_log[i],
                         queue='long',
                         machine='llrt3condor' )
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
    parser.add_argument('--debug', action='store_true', help='debug verbosity')
    args = parser.parse_args()

    submitTriggerEff( args )
