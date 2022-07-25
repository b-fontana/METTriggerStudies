###### DOCSTRING ####################################################
# Submits all the jobs required to obtain the trigger scale factors
# Run example:
####################################################################

import sys
sys.path.append("..")

import os
import argparse
import ROOT

from utils import utils
from scripts.jobWriter import JobWriter

@utils.set_pure_input_namespace
def writeHTCondorUnionWeightsCalculatorFiles_outputs(args):
    """
    Outputs are guaranteed to have the same length.
    Returns all separate paths to avoid code duplication.
    """
    base_name = 'UnionWeightsCalculator'
    data_folders = [ os.path.join( base_name, proc) for proc in args.mc_processes ]
    names = [ base_name + '_' + x for x in args.mc_processes ]
    return JobWriter.define_output( localdir=args.localdir,
                                    data_folders=data_folders,
                                    tag=args.tag,
                                    names=names )

@utils.set_pure_input_namespace
def writeHTCondorUnionWeightsCalculatorFiles(args):
    prog = utils.build_prog_path(args.localdir, 'runUnionWeightsCalculator.py')
    jobs, subs, checks = writeHTCondorUnionWeightsCalculatorFiles_outputs(args)
    jw = JobWriter()

    for i,proc in enumerate(args.mc_processes):
        filelist, inputdir = utils.get_root_input_files(proc, args.indir_root)

        #### Write shell executable (python scripts must be wrapped in shell files to run on HTCondor)
        command =  ( '{prog} --indir_root {indir_root} '.format(prog=prog, indir_root=inputdir)
                     + '--indir_json {indir_json} '.format(indir_json=args.indir_json)
                     + '--indir_eff {indir_eff} '.format(indir_eff=args.indir_eff)
                     + '--outdir {outr} '.format(outr=args.outdir)
                     + '--outprefix {outprefix} '.format(outprefix=args.outprefix)
                     + '--sample {sample} '.format(sample=proc)
                     + '--channels {channels} '.format(channels=' '.join(args.channels,))
                     + '--triggers {triggers} '.format(triggers=' '.join(args.triggers,))
                     + '--file_name ${1} '
                     + '--closure_single_trigger ${2} '
                     + '--variables {variables} '.format(variables=' '.join(args.variables,))
                     + '--tag {tag} '.format(tag=args.tag)
                     + '--subtag {subtag} '.format(subtag=args.subtag)
                     + '--data_name {dataname} '.format(dataname=args.data_name)
                     + '--mc_name {mcname} '.format(mcname=args.mc_name)
                     + '--binedges_fname {be}'.format(be=args.binedges_filename)
                    )

        if args.debug:
            command += '--debug '

        jw.write_shell(filename=jobs[i], command=command, localdir=args.localdir)
        jw.add_string('echo "Process {} done."'.format(proc))

        #### Write submission file
        jw.write_condor( filename=subs[i],
                         executable=jobs[i],
                         outfile=checks[i],
                         queue='short',
                         machine='llrt3condor7' )

        qlines = []
        for listname in filelist:
            for trig in args.closure_single_triggers:
                qlines.append('  {},{}'.format( os.path.basename(listname).replace('\n',''), trig ))
                        
        jw.write_queue( qvars=('filename', 'closure_single_trigger'),
                        qlines=qlines )
