import os
import sys
from utils import utils
from scripts.jobWriter import JobWriter

@utils.setPureInputNamespace
def runHaddCounts_outputs(args):
    targets = []

    # add the merge of all the samples first
    _tbase1, _tbase2 = utils.hadd_subpaths(args)
    tbase = _tbase1 + _tbase2
    t = os.path.join( args.indir, tbase + '.txt' )
    targets.append( t )

    # add individual sample merges
    for smpl in args.samples:
        tbase = _tbase1 + '_' + smpl + _tbase2
        t = os.path.join( args.indir, tbase + '.txt' )
        targets.append( t )
        
    return targets

@utils.setPureInputNamespace
def writeHTCondorHaddCountsFiles_outputs(args):
    return JobWriter.define_output( localdir=args.localdir,
                                    data_folders=['HaddCounts' + args.dataset_name,
                                                  'HaddCountsAgg' + args.dataset_name],
                                    tag=args.tag )

@utils.setPureInputNamespace
def writeHTCondorHaddCountsFiles(args):
    """Adds TXT count files"""
    targets = runHaddCounts_outputs(args)
    outs_job, outs_submit, outs_check = writeHTCondorHaddCountsFiles_outputs(args)
    jw = JobWriter()
    
    prog = utils.build_prog_path(args.localdir, 'addTriggerCounts.py')
    command_base =  ( '{prog} --indir {indir} '.format( prog=prog, indir=args.indir) +
                      '--outdir {outdir} '.format(outdir=args.outdir) +
                      '--subtag {subtag} '.format(subtag=args.subtag) +
                      '--tprefix {tprefix} '.format(tprefix=args.tprefix) +
                      '--dataset_name {dn} '.format(dn=args.dataset_name) +
                      '--outfile_counts ${1} '
                     )

    command_first_step = ( command_base +
                           '--sample ${2} ' +
                           ' --aggregation_step 0' )
    command_aggregation_step = ( command_base + '--infile_counts ${2} --aggregation_step 1')
    
    #### Write shell executable (python scripts must be wrapped in shell files to run on HTCondor)
    for out in outs_job:
        if out == outs_job[0]:
            jw.write_init(out, command_first_step, args.localdir)
            jw.add_string('echo "HaddCounts {} done."'.format(args.dataset_name))
        elif out == outs_job[1]:
            jw.write_init(out, command_aggregation_step, args.localdir)
            jw.add_string('echo "HaddCounts Agg {} done."'.format(args.dataset_name))

    #### Write submission file
    inputs_join = []

    for out1,out2,out3 in zip(outs_job,outs_submit,outs_check):
        jw.write_init( filename=out2,
                       executable=out1,
                       outfile=out3,
                       queue='short' )

        qvars = None
        qlines = []
        if out1 == outs_job[0]:
            qvars = ('myoutput', 'sample')
            for t,smpl in zip(targets[1:], args.samples):
                inputs = os.path.join(args.indir, smpl, args.tprefix + '*' + args.subtag + '.txt')
                inputs_join.append(t)
                qlines.append('  {}, {}'.format(t,smpl))

        elif out1 == outs_job[1]:
            qvars = ('myoutput', 'myinputs')
            qlines.append('  {}, {}'.format(targets[0], ' '.join(inputs_join)))
            
        jw.write_queue( qvars=qvars, qlines=qlines )
