import os
import sys
import sys
sys.path.append( os.path.join(os.environ['CMSSW_BASE'], 'src', 'METTriggerStudies'))
from utils.utils import (
    build_prog_path,
    hadd_subpaths,
    join_strings,
    set_pure_input_namespace,
    )
from condor.jobWriter import JobWriter

@set_pure_input_namespace
def runHaddCounts_outputs(args):
    targets = []

    # add the merge of all the samples first
    _tbase1, _tbase2 = hadd_subpaths(args)
    tbase = _tbase1 + _tbase2
    for chn in args.channels:
        t = os.path.join( args.indir, tbase + '_' + chn + '.csv' )
        targets.append( t )

    # add individual sample merges
    for smpl in args.samples:
        tbase = _tbase1 + '_' + smpl + _tbase2
        for chn in args.channels:
            t = os.path.join( args.indir, tbase + '_' + chn + '.csv' )
            targets.append( t )
    return targets

@set_pure_input_namespace
def writeHTCondorHaddCountsFiles_outputs(args):
    return JobWriter.define_output( localdir=args.localdir,
                                    data_folders=['HaddCounts' + args.dataset_name,
                                                  'HaddCountsAgg' + args.dataset_name],
                                    tag=args.tag )

@set_pure_input_namespace
def writeHTCondorHaddCountsFiles(args):
    """Adds TXT count files"""
    targets = runHaddCounts_outputs(args)
    outs_job, outs_submit, outs_check = writeHTCondorHaddCountsFiles_outputs(args)
    jw = JobWriter()
    
    prog = build_prog_path(args.localdir, 'addTriggerCounts.py')
    command_base = join_strings( '{} '.format(prog),
                                 '--indir {} '.format(args.indir),
                                 '--outdir {} '.format(args.outdir),
                                 '--subtag {} '.format(args.subtag),
                                 '--tprefix {} '.format(args.tprefix),
                                 '--dataset_name {} '.format(args.dataset_name),
                                 '--outfile_counts ${1} ' )

    command_first_step = ( command_base +
                           '--sample ${2} ' +
                           '--channel ${3} ' +
                           '--aggregation_step 0' )
    command_aggregation_step = ( command_base + '--infile_counts ${2} --channel ${3} --aggregation_step 1')
    
    #### Write shell executable (python scripts must be wrapped in shell files to run on HTCondor)
    for out in outs_job:
        if out == outs_job[0]:
            jw.write_shell(filename=out, command=command_first_step, localdir=args.localdir)
            jw.add_string('echo "HaddCounts {} done."'.format(args.dataset_name))
        elif out == outs_job[1]:
            jw.write_shell(filename=out, command=command_aggregation_step, localdir=args.localdir)
            jw.add_string('echo "HaddCounts Agg {} done."'.format(args.dataset_name))

    #### Write submission file
    inputs_join = {}
    nchannels = len(args.channels)
    for out1,out2,out3 in zip(outs_job,outs_submit,outs_check):
        jw.write_condor( filename=out2,
                         executable=out1,
                         outfile=out3,
                         queue='short',
                         machine='llrt3condor' )

        qvars = None
        qlines = []
        if out1 == outs_job[0]:
            qvars = ('myoutput', 'channel', 'sample')
            for it,t in enumerate(targets[nchannels:]):
                smpl = args.samples[ int(it/nchannels) ]
                chn = args.channels[ int(it%nchannels) ]
                if chn not in inputs_join:
                    inputs_join[chn] = []
                inputs_join[chn].append(t)
                qlines.append('  {}, {}, {}'.format(t,smpl,chn))

        elif out1 == outs_job[1]:
            qvars = ('myoutput', 'myinputs', 'channel')
            for ichn,chn in enumerate(args.channels):
                qlines.append('  {}, {}, {}'.format(targets[ichn], chn, ' '.join(inputs_join[chn])))

        jw.write_queue( qvars=qvars, qlines=qlines )
