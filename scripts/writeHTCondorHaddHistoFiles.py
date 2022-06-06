import os
import sys
from utils import utils
from scripts.jobWriter import JobWriter

@utils.set_pure_input_namespace
def runHaddHisto_outputs(args):
    targets = []

    # add the merge of all the samples first
    _tbase1, _tbase2 = utils.hadd_subpaths(args)
    tbase = _tbase1 + _tbase2
    t = os.path.join( args.indir, tbase + '.root' )
    targets.append( t )

    # add individual sample merges
    for smpl in args.samples:
        tbase = _tbase1 + '_' + smpl + _tbase2
        t = os.path.join( args.indir, tbase + '.root' )
        targets.append( t )
        
    return targets

@utils.set_pure_input_namespace
def writeHTCondorHaddHistoFiles_outputs(args):
    """
    Outputs are guaranteed to have the same length.
    Returns all separate paths to avoid code duplication.
    """
    return JobWriter.define_output( localdir=args.localdir,
                                    data_folders=[ 'HaddHisto' + args.dataset_name,
                                                   'HaddHistoAgg' + args.dataset_name],
                                    tag=args.tag )

@utils.set_pure_input_namespace
def writeHTCondorHaddHistoFiles(args):
    """Adds ROOT histograms"""
    targets = runHaddHisto_outputs(args)
    outs_job, outs_submit, outs_check = writeHTCondorHaddHistoFiles_outputs(args)
    jw = JobWriter()
    
    #### Write shell executable (python scripts must be wrapped in shell files to run on HTCondor)
    command = 'hadd -f ${1} ${@:2}' #bash: ${@:POS} captures all arguments starting from POS

    for out in outs_job:
        jw.write_init(out, command, args.localdir)
        if out == outs_job[0]:
            jw.add_string('echo "HaddHisto {} done."'.format(args.dataset_name))
        elif out == outs_job[1]:
            jw.add_string('echo "HaddHisto Agg {} done."'.format(args.dataset_name))

    #### Write submission file
    inputs_join = []
    for out1,out2,out3 in zip(outs_job,outs_submit,outs_check):
        jw.write_init( filename=out2,
                       executable=out1,
                       outfile=out3,
                       queue='short' )

        qlines = []
        if out1 == outs_job[0]:
            for t,smpl in zip(targets[1:], args.samples):
                inputs = os.path.join(args.indir, smpl, args.tprefix + '*' + args.subtag + '.root')
                inputs_join.append(t)
                # join subdatasets (different MC or Data subfolders, ex: TT_fullyHad, TT_semiLep, ...)
                qlines.append('  {}, {}'.format(t, inputs))
        elif out1 == outs_job[1]:
            # join MC or Data subdatasets into a single one (ex: TT)
            qlines.append('  {}, {}'.format(targets[0], ' '.join(inputs_join)))
        
        jw.write_queue( qvars=('myoutput', 'myinputs'),
                        qlines=qlines )
