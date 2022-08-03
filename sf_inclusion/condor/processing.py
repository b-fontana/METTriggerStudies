# coding: utf-8

_all_ = [ "processing", "processing_outputs" ]

import sys
sys.path.append("..")

import os
import re
import argparse
import ROOT

from utils import utils
from condor.job_writer import JobWriter

def produce_trigger_outputs_sample(args, sample, ext):
    """
    Produces all outputs of the submitTriggerEff task.
    Limitation: As soon as one file is not produced, luigi
    reruns everything.
    """
    assert(ext in ('root', 'txt'))
    extension = '.' + ext
    t = []
    exp = re.compile('.+output(_[0-9]{1,5}).root')

    inputs, _ = utils.get_root_input_files(sample, args.indir)

    folder = os.path.join( args.outdir, proc )
    for inp in inputs:
        number = exp.search(inp)
        proc_folder = os.path.dirname(inp).split('/')[-1]
        basename = args.tprefix + '_' + proc_folder + number.group(1)
        basename += args.subtag + extension
        t.append( os.path.join(folder, basename) )
    return t
    
@utils.set_pure_input_namespace
def produce_trigger_outputs(args, ext='root'):
    """
    Produces all outputs of the submitTriggerEff task.
    Limitation: As soon as one file is not produced, luigi
    reruns everything.
    """
    tdata, tmc = ([] for _ in range(2))
    for proc in args.data_vals:
        tdata.extend( produce_trigger_outputs_sample(args, proc, ext) )
    for proc in args.mc_vals:
        tmc.extend( produce_trigger_outputs_sample(args, proc, ext) )
    return tdata, tmc

@utils.set_pure_input_namespace
def processing_outputs(args):
    if args.mode == 'histos':
        name = 'Histos'
    elif args.mode == 'counts':
        name = 'Counts'
    else:
        raise ValueError('Mode {} is not supported.'.format(args.mode))

    _data_tup = tuple((k,v) for k,v in zip(args.data_keys,args.data_vals))
    _mc_tup = tuple((k,v) for k,v in zip(args.mc_keys,args.mc_vals))

    data_folders = [ name + '_' + v for v in args.data_vals ]
    mc_folders   = [ name + '_' + v for v in args.mc_vals ]
    job_opt = dict(localdir=args.localdir, tag=args.tag)
    return ( JobWriter.define_output( data_folders=data_folders, **job_opt),
             JobWriter.define_output( data_folders=mc_folders,   **job_opt),
             _data_tup, _mc_tup )

@utils.set_pure_input_namespace
def processing(args):
    prog = utils.build_prog_path(args.localdir,
                                 ('produceTriggerHistograms.py' if args.mode == 'histos'
                                  else 'produceTriggerCounts.py'))
    jw = JobWriter()

    outs_data, outs_mc, _data_procs, _mc_procs = processing_outputs(args)

    # unite Data and MC lists
    outs_job    = outs_data[0] + outs_mc[0]
    outs_submit = outs_data[1] + outs_mc[1]
    outs_check  = outs_data[2] + outs_mc[2]
    outs_log    = outs_data[3] + outs_mc[3]
    _all_processes = _data_procs + _mc_procs

    for i, (kproc, vproc) in enumerate(_all_processes):
        filelist, _ = utils.get_root_input_files(vproc, args.indir)
        
        #### Write shell executable (python scripts must be wrapped in shell files to run on HTCondor)
        command =  ( '{} '            .format(prog) +
                     '--outdir {} '   .format(args.outdir) +
                     '--dataset {} '  .format(kproc) +
                     '--sample {} '   .format(vproc) +
                     '--isdata {} '   .format(int(vproc in args.data_vals)) +
                     '--file ${1} ' +
                     '--subtag {} '   .format(args.subtag) +
                     '--channels {} ' .format(' '.join(args.channels,)) +
                     '--triggers {} ' .format(' '.join(args.triggers,)) +
                     '--tprefix {} '  .format(args.tprefix)
                    )
        
        if args.debug:
            command += '--debug '

        if args.mode == 'histos':
            command += ( '--binedges_fname {} '.format(args.binedges_filename) +
                         '--intersection_str {} '.format(args.intersection_str) +
                         '--variables {} '.format(' '.join(args.variables,)) +
                         '--nocut_dummy_str {}'.format(args.nocut_dummy_str)
                        )

        jw.write_shell(filename=outs_job[i], command=command, localdir=args.localdir)
        jw.add_string('echo "Process {} done in mode {}."'.format(vproc,args.mode))

        #### Write submission file
        jw.write_condor( filename=outs_submit[i],
                         executable=outs_job[i],
                         outfile=outs_check[i],
                         logfile=outs_log[i],
                         queue='long',
                         machine='llrt3condor' )
        
        qlines = []
        for listname in filelist:
            qlines.append(' {}'.format( listname.replace('\n','') ))
        
        jw.write_queue( qvars=('filename',),
                        qlines=qlines )

# -- Parse options
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Command line parser')

    parser.add_argument('--binedges_dataset', dest='binedges_dataset', required=True, help='in directory')
    parser.add_argument('--localdir', dest='localdir', default=os.getcwd(), help='out directory')
    parser.add_argument('--indir', dest='indir', required=True, help='in directory')
    parser.add_argument('--outdir', dest='outdir', required=True, help='out directory')
    parser.add_argument('--tag', dest='tag', required=True, help='tag')
    parser.add_argument('--subtag', dest='subtag', required=True, help='subtag')
    parser.add_argument('--tprefix', dest='tprefix', required=True, help='target prefix')
    parser.add_argument('--mc_processes', dest='mc_processes', required=True, nargs='+', type=str,
                        help='list of MC process names')                
    parser.add_argument('--data_keys', dest='data_keys', required=True, nargs='+', type=str,
                        help='list of datasets')
    parser.add_argument('--data_vals', dest='data_vals', required=True, nargs='+', type=str,
                        help='list of datasets')
    parser.add_argument('--channels',   dest='channels', required=True, nargs='+', type=str,
                        help='Select the channels over which the workflow will be run.' )
    parser.add_argument('--triggers', dest='triggers', required=True, nargs='+', type=str,
                        help='Select the triggers over which the workflow will be run.' )
    parser.add_argument('--variables', dest='variables', required=True, nargs='+', type=str,
                        help='Select the variables over which the workflow will be run.' )
    parser.add_argument('--intersection_str', dest='intersection_str', required=False, default='_PLUS_',
                        help='String used to represent set intersection between triggers.')
    parser.add_argument('--nocut_dummy_str', dest='nocut_dummy_str', required=True,
                        help='Dummy string associated to trigger histograms were no cuts are applied.')
    parser.add_argument('--debug', action='store_true', help='debug verbosity')
    args = parser.parse_args()

    submitTriggerEff( args )
