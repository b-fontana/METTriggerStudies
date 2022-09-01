# coding: utf-8

_all_ = [ ]
    
import os
import sys
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, parent_dir)

import re
import luigi
import time
import argparse
from argparse import RawTextHelpFormatter

import config
from config import main
from utils import utils, luigi_utils as lutils
from scripts import def_bins
from condor import (
    closure,
    dag,
    discriminator,   
    eff_and_sf,
    eff_and_sf_aggr,
    hadd_counts,
    hadd_histo,
    job_writer,
    processing,
    union_calculator,
    )

descr = 'Run example: `copython inclusion/run.py --tag abc --data MET EG --mc_process TT`'
parser = argparse.ArgumentParser(description=descr, formatter_class=RawTextHelpFormatter)
parser.add_argument(
    '--nbins',
    type=int,
    default=6,
    help="Number of histogram bins. If fine-grained control is required modify the variable `binedges` in the configuration file."
    )
parser.add_argument(
    '--workers',
    type=int,
    default=1,
    help="Maximum number of worker which can be used to run the pipeline."
    )
parser.add_argument(
    '--scheduler',
    type=str,
    choices=['local', 'central'],
    default='local',
    help='Select the scheduler for luigi.'
    )
parser.add_argument(
    '--data',
    type=str,
    nargs='+',
    required=True,
    choices=main.data.keys(),
    help='Select the data over which the workflow will be run.'
    )
parser.add_argument(
    '--mc_processes',
    type=str,
    nargs='+',
    required=True,
    choices=main.mc_processes.keys(),
    help='Select the MC processes over which the workflow will be run.'
    )
parser.add_argument(
    '--triggers_closure',
    nargs='+', #1 or more arguments
    type=str,
    required=False,
    default=main.triggers,
    choices=main.triggers,
    help=( 'Select the triggers considered for the closure.' +
          'The default is to consider all of them.\n' +
          'To do a closure for one single trigger efficiency (which ' +
          'should provide a perfect match between original and ' +
          'weighted MC), one must specify here the trigger to consider.' )
          )
parser.add_argument(
    '--channels',
    nargs='+', #1 or more arguments
    type=str,
    default=main.channels,
    help='Select the channels over which the workflow will be run.'
    )
parser.add_argument(
    '--variables_for_efficiencies',
    nargs='+', #1 or more arguments
    type=str,
    default=main.var_eff,
    help='Select the variables to be used for the calculation of efficiencies.'
    )
parser.add_argument(
    '--variables_for_distributions',
    nargs='+', #1 or more arguments
    type=str,
    default=main.var_dist,
    help='Select the variables to be used for the display of distributions.'
    )
parser.add_argument(
    '-t',
    '--tag',
    type=str,
    required=True,
    help='Specifies a tag to differentiate workflow runs.'
    )
parser.add_argument(
    '--subtag',
    type=str,
    default='default',
    help='Specifies a subtag, for instance an additional cut within the same tag. We force its first character to be an underscore.'
    )
parser.add_argument(
    '--distributions',
    type=int,
    choices=[0,1,2],
    default=0,
    help="0: Does not draw the distributions (default).\n1: Also draws the distributions.\n2: Only draws the distributions."
    )
parser.add_argument(
    '--counts',
    action='store_true',
    help="Only runs the 'counting' workflow: check how many events pass each intersection of triggers. The default is to run the full workflow."
    )
parser.add_argument(
    '--configuration',
    type=str,
    default='sel_default',
    help='Specifies a subtag, for instance an additional cut within the same tag. We force its first character to be an underscore.'
    )
parser.add_argument(
    '--debug_workflow',
    action='store_true',
    help="Explicitly print the functions being run for each task, for workflow debugging purposes."
    )
FLAGS, _ = parser.parse_known_args()
assert set(FLAGS.triggers_closure).issubset(set(main.triggers))
         
user_data = {k:v for k,v in main.data.items()
             if k in FLAGS.data}
user_mc   = {k:v for k,v in main.mc_processes.items()
             if k in FLAGS.mc_processes}

data_keys, data_vals = utils.flatten_nested_dict(user_data)
mc_keys, mc_vals = utils.flatten_nested_dict(user_mc)
assert len(data_keys)==len(data_vals)
assert len(mc_keys)==len(mc_vals)

data_name = 'Data_' + '_'.join(FLAGS.data)
mc_name   = 'MC_'   + '_'.join(FLAGS.mc_processes)

data_storage = os.path.join(main.storage, FLAGS.tag, 'Data' )
out_storage = os.path.join(main.storage, FLAGS.tag, 'Outputs')
targets_folder = os.path.join(data_storage, 'targets')    
binedges_filename = os.path.join(data_storage, 'binedges.hdf5')

subtag = ( FLAGS.subtag if FLAGS.subtag==''
           else ( '_' + FLAGS.subtag if FLAGS.subtag[0] != '_' else FLAGS.subtag ) )


variables_join = tuple(set(FLAGS.variables_for_efficiencies + FLAGS.variables_for_distributions))

#### scripts/def_bins
bins_params = {'nbins'             : FLAGS.nbins,
               'binedges_filename' : binedges_filename,
               'indir'             : main.inputs,
               'outdir'            : data_storage,
               'data_vals'         : data_vals,
               'variables'         : variables_join,
               'channels'          : FLAGS.channels,
               'tag'               : FLAGS.tag,
               'subtag'            : subtag,
               'debug'             : FLAGS.debug_workflow}


#### condor/dag
write_params = {'data_name' : data_name,
                'localdir'  : main.base_folder,
                'tag'       : FLAGS.tag}

#### scripts/produce_trig_histos
histos_params = {'binedges_filename' : binedges_filename,
                 'indir'             : main.inputs,
                 'outdir'            : data_storage,
                 'localdir'          : main.base_folder,
                 'data_keys'         : data_keys,
                 'data_vals'         : data_vals,
                 'mc_keys'           : mc_keys,
                 'mc_vals'           : mc_vals,
                 'triggers'          : main.triggers,
                 'channels'          : FLAGS.channels,
                 'variables'         : variables_join,
                 'tag'               : FLAGS.tag,
                 'subtag'            : subtag,
                 'intersection_str'  : main.inters_str,
                 'nocut_dummy_str'   : main.nocut_dummy,
                 'configuration'     : 'inclusion.main.' + FLAGS.configuration,
                 'debug'             : FLAGS.debug_workflow}

#### scripts/hadd_histo
haddhisto_params = {'indir'    : data_storage,
                    'localdir' : main.base_folder,
                    'tag'      : FLAGS.tag,
                    'subtag'   : subtag, }

#### scripts/add_counts
haddcounts_params = {'indir'    : data_storage,
                     'outdir'   : out_storage,
                     'localdir'  : main.base_folder,
                     'tag'      : FLAGS.tag,
                     'subtag'   : subtag,
                     'channels' : FLAGS.channels, }

#### drawTriggerScaleFactors
sf_params = {'data_name'            : data_name,
             'mc_name'              : mc_name,
             'draw_independent_MCs' : False,
             'indir'                : data_storage,
             'outdir'               : out_storage,
             'localdir'             : main.base_folder,
             'triggers'             : main.triggers,
             'channels'             : FLAGS.channels,
             'variables'            : FLAGS.variables_for_efficiencies,
             'binedges_filename'    : binedges_filename,
             'tag'                  : FLAGS.tag,
             'subtag'               : subtag,
             'canvas_prefix'        : main.pref['canvas'],
             'intersection_str'     : main.inters_str,
             'nocut_dummy_str'      : main.nocut_dummy,
             'debug'                : FLAGS.debug_workflow,}

sfagg_params = {'indir'       : out_storage,
                'outdir'      : out_storage,
                'localdir'    : main.base_folder,
                'channels'    : FLAGS.channels,
                'variables'   : FLAGS.variables_for_efficiencies,
                'tag'         : FLAGS.tag,
                'subtag'      : subtag,
                'file_prefix' : main.pref['sf'],
                'debug'       : FLAGS.debug_workflow,}

#### scripts/discriminator
discriminator_params = {'indir'            : data_storage,
                        'outdir'           : data_storage,
                        'localdir'         : main.base_folder,
                        'triggers'         : main.triggers,
                        'channels'         : FLAGS.channels,
                        'variables'        : FLAGS.variables_for_efficiencies,
                        'tag'              : FLAGS.tag,
                        'subtag'           : subtag,
                        'intersection_str' : main.inters_str,
                        'debug'            : FLAGS.debug_workflow,}

#### scripts/calculator
calculator_params = {'binedges_filename'       : binedges_filename,
                     'indir_root'              : main.inputs,
                     'indir_json'              : data_storage,
                     'indir_eff'               : out_storage,
                     'outdir'                  : data_storage,
                     'outprefix'               : main.pref['clos'],
                     'data_name'               : data_name,
                     'mc_name'                 : mc_name,
                     'mc_processes'            : mc_vals,
                     'localdir'                : main.base_folder,
                     'triggers'                : main.triggers,
                     'closure_single_triggers' : FLAGS.triggers_closure,
                     'channels'                : FLAGS.channels,
                     'variables'               : FLAGS.variables_for_efficiencies,
                     'tag'                     : FLAGS.tag,
                     'subtag'                  : subtag,
                     'configuration'           : 'inclusion.main.' + FLAGS.configuration,
                     'debug'                   : FLAGS.debug_workflow,}

#### scripts/closure
closure_params = {'data_name'               : data_name,
                  'mc_name'                 : mc_name,
                  'binedges_filename'       : binedges_filename,
                  'indir_eff'               : out_storage,
                  'indir_union'             : data_storage,
                  'indir_json'              : data_storage,
                  'outdir'                  : out_storage,
                  'inprefix'                : main.pref['clos'],
                  'eff_prefix'              : main.pref['sf'],
                  'mc_processes'            : mc_vals,
                  'out_weighted_prefix'     : main.pref['clos'],
                  'out_original_prefix'     : main.pref['histos'],
                  'localdir'                : main.base_folder,
                  'closure_single_triggers' : FLAGS.triggers_closure,
                  'channels'                : FLAGS.channels,
                  'variables'               : FLAGS.variables_for_efficiencies,
                  'tag'                     : FLAGS.tag,
                  'subtag'                  : subtag,
                  'debug'                   : FLAGS.debug_workflow }

#### Helper functions

def get_target_path(taskname, targets_dir):
    re_txt = re.compile('\.txt')
    target_path = os.path.join(targets_dir,
                               re_txt.sub('_'+taskname+'.txt', main.targ_def)) 
    return target_path


#### Tasks
 
class DefineBinning(luigi.Task):
    """Calculate the most adequate binning based on data."""
    args = utils.dot_dict(bins_params)
    
    @lutils.WorkflowDebugger(flag=FLAGS.debug_workflow)
    def output(self):
        target = def_bins.define_binning_outputs(self.args)
 
        #write the target files for debugging
        target_path = get_target_path(self.__class__.__name__, targets_folder)
        utils.remove(target_path)
        with open(target_path, 'w') as f:
            f.write( target + '\n' )
 
        return luigi.LocalTarget(target)
 
    @lutils.WorkflowDebugger(flag=FLAGS.debug_workflow)
    def run(self):
        def_bins.define_binning(self.args)
 
 
class Processing(lutils.ForceRun):
    """Write htcondor files for total and passed trigger histograms."""
    params = utils.dot_dict(histos_params)
 
    mode = luigi.ChoiceParameter(choices=('histos', 'counts'),
                                 var_type=str)
    
    @lutils.WorkflowDebugger(flag=FLAGS.debug_workflow)
    def output(self):
        self.params['mode'] = self.mode
        self.params['tprefix'] = main.pref[self.mode]
        obj_data, obj_mc, _, _ = processing.processing_outputs(self.params)
        o1_1, o1_2, _, _ = obj_data
        o2_1, o2_2, _, _ = obj_mc
 
        target_path = get_target_path(self.__class__.__name__,
                                             targets_folder)
        utils.remove( target_path )
 
        objs = (o1_1, o1_2, o2_1, o2_2)
        with open( target_path, 'w' ) as f:
            for x in objs:
                for t in x:
                    f.write(t + '\n')
 
        return sum([lutils.convert_to_luigi_local_targets(x) for x in objs], [])
    
    @lutils.WorkflowDebugger(flag=FLAGS.debug_workflow)
    def run(self):
        self.params['mode'] = self.mode
        self.params['tprefix'] = main.pref[self.mode]
        processing.processing(self.params)
 
    @lutils.WorkflowDebugger(flag=FLAGS.debug_workflow)
    def requires(self):
        return DefineBinning()
 
 
class HaddHisto(lutils.ForceRun):
    """Write htcondor files for hadd'ing histograms in root files."""
    samples = luigi.ListParameter()
    dataset_name = luigi.Parameter()
    args = utils.dot_dict(haddhisto_params)
    args['tprefix'] = main.pref['histos']
    
    @lutils.WorkflowDebugger(flag=FLAGS.debug_workflow)
    def output(self):
        self.args['samples'] = lutils.luigi_to_raw( self.samples )
        self.args['dataset_name'] = self.dataset_name
        o1, o2, _, _ = hadd_histo.hadd_histo_outputs(self.args)
        
        target_path = get_target_path(self.__class__.__name__, targets_folder)
        utils.remove(target_path)
        with open(target_path, 'w') as f:
            for t in o1: f.write( t + '\n' )
            for t in o2: f.write( t + '\n' )
 
        _c1 = lutils.convert_to_luigi_local_targets(o1)
        _c2 = lutils.convert_to_luigi_local_targets(o2)
        return _c1 + _c2
 
    @lutils.WorkflowDebugger(flag=FLAGS.debug_workflow)
    def run(self):
        self.args['samples'] = lutils.luigi_to_raw( self.samples )
        self.args['dataset_name'] = self.dataset_name
        hadd_histo.hadd_histo( self.args )
 
 
class HaddCounts(lutils.ForceRun):
    """Write htcondor files for hadding txt count files."""
    samples = luigi.ListParameter()
    dataset_name = luigi.Parameter()
    args = utils.dot_dict(haddcounts_params)
    args['tprefix'] = main.pref['counts']
    
    @lutils.WorkflowDebugger(flag=FLAGS.debug_workflow)
    def output(self):
        self.args['samples'] = lutils.luigi_to_raw( self.samples )
        self.args['dataset_name'] = self.dataset_name
        o1, o2, _, _ = hadd_counts.hadd_counts_outputs( self.args )
        
        target_path = get_target_path(self.__class__.__name__,
                                             targets_folder)
        utils.remove( target_path )
        with open( target_path, 'w' ) as f:
            for t in o1: f.write( t + '\n' )
            for t in o2: f.write( t + '\n' )
 
        _c1 = lutils.convert_to_luigi_local_targets(o1)
        _c2 = lutils.convert_to_luigi_local_targets(o2)
        return _c1 + _c2
 
    @lutils.WorkflowDebugger(flag=FLAGS.debug_workflow)
    def run(self):
        self.args['samples'] = lutils.luigi_to_raw(self.samples)
        self.args['dataset_name'] = self.dataset_name
        hadd_counts.hadd_counts( self.args )
 
 
class EffAndSF(lutils.ForceRun):
    """Write htcondor files for efficiencies and scale factors."""
    params = utils.dot_dict(sf_params)
    params['tprefix'] = main.pref['histos']
    
    @lutils.WorkflowDebugger(flag=FLAGS.debug_workflow)
    def output(self):
        o1, o2, _, _ = eff_and_sf.eff_and_sf_outputs(self.params)
 
        target_path = get_target_path(self.__class__.__name__,
                                             targets_folder)
        utils.remove( target_path )
        with open( target_path, 'w' ) as f:
            f.write( o1 + '\n' )
            f.write( o2 + '\n' )
 
        _c1 = lutils.convert_to_luigi_local_targets(o1)
        _c2 = lutils.convert_to_luigi_local_targets(o2)
        return _c1 + _c2
 
    @lutils.WorkflowDebugger(flag=FLAGS.debug_workflow)
    def run(self):
        eff_and_sf.eff_and_sf(self.params)
 
class EffAndSFAggr(lutils.ForceRun):
    """
    Aggregate htcondor files for efficiencies and scale factors.
    Useful for transfering the intersection efficiencies to the KLUB framework.
    Not needed for the following steps.
    """
    params = utils.dot_dict(sfagg_params)
    
    @lutils.WorkflowDebugger(flag=FLAGS.debug_workflow)
    def output(self):
        o1, o2, _, _ = eff_and_sf_aggr.eff_and_sf_aggr_outputs(self.params)
 
        target_path = get_target_path(self.__class__.__name__,
                                             targets_folder)
        utils.remove( target_path )
        with open( target_path, 'w' ) as f:
            f.write( o1 + '\n' )
            f.write( o2 + '\n' )
 
        _c1 = lutils.convert_to_luigi_local_targets(o1)
        _c2 = lutils.convert_to_luigi_local_targets(o2)
        return _c1 + _c2
 
    @lutils.WorkflowDebugger(flag=FLAGS.debug_workflow)
    def run(self):
        eff_and_sf_aggr.eff_and_sf_aggr(self.params)
 
 
class Discriminator(lutils.ForceRun):
    """Write htcondor files for the variable discriminator."""
    params = utils.dot_dict(discriminator_params)
    
    @lutils.WorkflowDebugger(flag=FLAGS.debug_workflow)
    def output(self):
        o1, o2, _, _ = discriminator.discriminator_outputs(self.params)
 
        target_path = get_target_path(self.__class__.__name__,
                                             targets_folder)
        utils.remove( target_path )
        with open( target_path, 'w' ) as f:
            for t in o1: f.write( t + '\n' )
            for t in o2: f.write( t + '\n' )
 
        _c1 = lutils.convert_to_luigi_local_targets(o1)
        _c2 = lutils.convert_to_luigi_local_targets(o2)
        return _c1 + _c2
 
    @lutils.WorkflowDebugger(flag=FLAGS.debug_workflow)
    def run(self):
        discriminator.discriminator(self.params)
 
 
class UnionCalculator(lutils.ForceRun):
    """Write htcondor files for scale factor calculator."""
    params = utils.dot_dict(calculator_params)
    
    @lutils.WorkflowDebugger(flag=FLAGS.debug_workflow)
    def output(self):
        o1, o2, _, _ = union_calculator.union_calculator_outputs(self.params)
 
        #write the target files for debugging
        target_path = get_target_path(self.__class__.__name__,targets_folder)
        utils.remove( target_path )
        with open( target_path, 'w' ) as f:
            for t in o1: f.write( t + '\n' )
            for t in o2: f.write( t + '\n' )
 
        _c1 = lutils.convert_to_luigi_local_targets(o1)
        _c2 = lutils.convert_to_luigi_local_targets(o2)
        return _c1 + _c2
 
    @lutils.WorkflowDebugger(flag=FLAGS.debug_workflow)
    def run(self):
        union_calculator.union_calculator(self.params)
 
 
class Closure(lutils.ForceRun):
    """Write htcondor files for displaying closure plots."""
    params = utils.dot_dict(closure_params)
    
    @lutils.WorkflowDebugger(flag=FLAGS.debug_workflow)
    def output(self):
        o1, o2, _, _ = closure.closure_outputs(self.params)
 
        target_path = get_target_path(self.__class__.__name__,targets_folder)
        utils.remove( target_path )
        with open( target_path, 'w' ) as f:
            f.write( o1 + '\n' )
            f.write( o2 + '\n' )
 
        _c1 = lutils.convert_to_luigi_local_targets(o1)
        _c2 = lutils.convert_to_luigi_local_targets(o2)
        return _c1 + _c2
 
    @lutils.WorkflowDebugger(flag=FLAGS.debug_workflow)
    def run(self):
        closure.closure(self.params)
 
 
class Dag(lutils.ForceRun):
    """Triggering all htcondor writing classes."""
    params        = utils.dot_dict(write_params)
    p_histos      = utils.dot_dict(histos_params)
    p_hadd_histo  = utils.dot_dict(haddhisto_params)
    p_hadd_counts = utils.dot_dict(haddhisto_params)
    p_eff_sf      = utils.dot_dict(sf_params)
    p_eff_sf_agg  = utils.dot_dict(sfagg_params)
    p_disc        = utils.dot_dict(discriminator_params)
    p_calc        = utils.dot_dict(calculator_params)
    p_closure     = utils.dot_dict(closure_params)
    
    p_hadd_histo['tprefix'] = main.pref['histos']
    p_eff_sf['tprefix'] = main.pref['histos']
    
    @lutils.WorkflowDebugger(flag=FLAGS.debug_workflow)
    def output(self):
        o1 = dag.dag_outputs(self.params)
 
        target_path = get_target_path(self.__class__.__name__,targets_folder)
        utils.remove( target_path )
        with open( target_path, 'w' ) as f:
            f.write( o1 + '\n' )
 
        _c1 = lutils.convert_to_luigi_local_targets(o1)
        return _c1
 
    @lutils.WorkflowDebugger(flag=FLAGS.debug_workflow)
    def run(self):
        self.p_histos['mode'] = 'histos'
        obj_hdata, obj_hmc, _, _ = processing.processing_outputs(self.p_histos)
        subm_hdata = obj_hdata[1]
        subm_hmc = obj_hmc[1]
 
        self.p_histos['mode'] = 'counts'
        obj_cdata, obj_cmc, _, _ = processing.processing_outputs(self.p_histos)
        subm_cdata = obj_cdata[1]
        subm_cmc = obj_cmc[1]
 
        self.p_hadd_histo['dataset_name']  = data_name
        subm_hadd_hdata = hadd_histo.hadd_histo_outputs(self.p_hadd_histo)[1]
        self.p_hadd_histo['dataset_name']  = mc_name
        subm_hadd_hmc = hadd_histo.hadd_histo_outputs(self.p_hadd_histo)[1]
 
        self.p_hadd_counts['dataset_name'] = data_name
        subm_hadd_cdata = hadd_counts.hadd_counts_outputs(self.p_hadd_counts)[1]
        self.p_hadd_counts['dataset_name'] = mc_name
        subm_hadd_cmc = hadd_counts.hadd_counts_outputs(self.p_hadd_counts)[1]
        
        subm_eff_sf = eff_and_sf.eff_and_sf_outputs(self.p_eff_sf)[1]
        subm_eff_sf_agg = eff_and_sf_aggr.eff_and_sf_aggr_outputs(self.p_eff_sf_agg)[1]
        subm_disc = discriminator.discriminator_outputs(self.p_disc)[1]
        # subm_union = union_calculator.union_calculator_outputs(self.p_calc)[1]
        # subm_closure = closure.closure_outputs(self.p_closure)
 
        jobs = { 'HistosData':     subm_hdata,
                 'HistosMC':       subm_hmc,
                 'CountsData':     subm_cdata,
                 'CountsMC':       subm_cmc,
                 'HaddHistoData':  subm_hadd_hdata,
                 'HaddHistoMC':    subm_hadd_hmc,
                 'HaddCountsData': subm_hadd_cdata,
                 'HaddCountsMC':   subm_hadd_cmc,
                 'EffSF':          [ subm_eff_sf ],
                 'EffSFAgg':       [ subm_eff_sf_agg ],
                 'Discr':          subm_disc,
                 #'Union':          subm_union,
                 #'Closure':        [ subm_closure ],
                }
 
        dag_manager = dag.WriteDAGManager( self.params['localdir'],
                                           self.params['tag'],
                                           jobs,
                                           mode='short' )
        dag_manager.write_all()
        
class SubmitDAG(lutils.ForceRun):
    """Submission class."""
    def edit_condor_submission_file(self, out):
        jw = job_writer.JobWriter()
        with open(out, 'r') as f:
            contents = f.readlines()
        ncontents = len(contents)
        new_content = jw.condor_specific_content(queue=main.queue,
                                                 machine='llrt3condor')
        contents.insert(ncontents-1, new_content + '\n')
        with open(out, 'w') as f:
            contents = ''.join(contents)
            f.write(contents)
        
    @lutils.WorkflowDebugger(flag=FLAGS.debug_workflow)
    def run(self):
        outfile = self.input()[-1][0].path
        com = 'condor_submit_dag -no_submit -f '
        com += '-outfile_dir {} {}'.format(os.path.dirname(outfile), outfile)
 
        os.system(com)
        time.sleep(.5)
        self.edit_condor_submission_file(outfile + '.condor.sub')
        time.sleep(.5)
        os.system('condor_submit {}.condor.sub'.format(outfile))
 
    @lutils.WorkflowDebugger(flag=FLAGS.debug_workflow)
    def output(self):
        # WriteDag dependency is the last one
        target = self.input()[-1][0].path + '.condor.sub'
        return luigi.LocalTarget(target)
 
    @lutils.WorkflowDebugger(flag=FLAGS.debug_workflow)
    def requires(self):
        return [ Processing(mode='histos'),
                 Processing(mode='counts'),
                 HaddHisto(dataset_name=data_name, samples=data_vals),
                 HaddHisto(dataset_name=mc_name, samples=mc_vals),
                 HaddCounts(dataset_name=data_name, samples=data_vals),
                 HaddCounts(dataset_name=mc_name, samples=mc_vals ),
                 EffAndSF(),
                 EffAndSFAggr(),
                 Discriminator(),
                 UnionCalculator(),
                 Closure(),
                 Dag(),
                ]
   
 
utils.create_single_dir( data_storage )
utils.create_single_dir( targets_folder )
    
last_tasks = [ SubmitDAG() ]
 
if FLAGS.scheduler == 'central':
    luigi.build(last_tasks,
                workers=FLAGS.workers, local_scheduler=False, log_level='INFO')
if FLAGS.scheduler == 'local':
    luigi.build(last_tasks,
                local_scheduler=True, log_level='INFO')
 
else:
    raise RuntimeError('This script can only be run directly from the command line.')
