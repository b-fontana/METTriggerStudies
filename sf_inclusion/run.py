# coding: utf-8

_all_ = [ ]
    
import os
import time
import inspect
import re

from utils import utils
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

import luigi
from luigi_conf import luigi_utils as lutils
from luigi_conf.luigi_cfg import cfg, FLAGS
lcfg = cfg() #luigi configuration


### Helper functions 

def convert_to_luigi_local_targets(targets):
    """Converts a list of files into a list of luigi targets."""
    if not isinstance(targets, (tuple,list,set)):
        targets = [ targets ]
    return [ luigi.LocalTarget(t) for t in targets ]

def get_target_path(taskname):
    re_txt = re.compile('\.txt')
    target_path = os.path.join(lcfg.targets_folder,
                               re_txt.sub( '_'+taskname+'.txt', lcfg.targets_default_name) ) 
    return target_path

def luigi_to_raw( param ):
    """
    Converts luigi parameters to their "raw" versions.
    Useful when passing them to functions.
    """
    if isinstance(param, (tuple,list)):
        return [x for x in param]
    else:
        raise NotImplementedError('[' + inspect.stack()[0][3] + ']: ' + 'only tuples/lists implemented so far!')

### Tasks

class DefineBinning(luigi.Task):
    """Calculate the most adequate binning based on data."""
    args = utils.dot_dict(lcfg.bins_params)
    
    @lutils.WorkflowDebugger(flag=FLAGS.debug_workflow)
    def output(self):
        target = def_bins.define_binning_outputs( self.args )

        #write the target files for debugging
        target_path = get_target_path( self.__class__.__name__ )
        utils.remove( target_path )
        with open( target_path, 'w' ) as f:
            f.write( target + '\n' )

        return luigi.LocalTarget(target)

    @lutils.WorkflowDebugger(flag=FLAGS.debug_workflow)
    def run(self):
        def_bins.define_binning( self.args )


class Processing(lutils.ForceRun):
    """Write htcondor files for total and passed trigger histograms."""
    params = utils.dot_dict(lcfg.histos_params)

    mode = luigi.ChoiceParameter(choices=lcfg.modes.keys(),
                                 var_type=str)
    
    @lutils.WorkflowDebugger(flag=FLAGS.debug_workflow)
    def output(self):
        self.params['mode'] = self.mode
        self.params['tprefix'] = lcfg.modes[self.mode]
        obj_data, obj_mc, _, _ = processing.processing_outputs(self.params)
        o1_1, o1_2, _, _ = obj_data
        o2_1, o2_2, _, _ = obj_mc

        target_path = get_target_path( self.__class__.__name__ )
        utils.remove( target_path )

        objs = (o1_1, o1_2, o2_1, o2_2)
        with open( target_path, 'w' ) as f:
            for x in objs:
                for t in x:
                    f.write(t + '\n')

        return sum([convert_to_luigi_local_targets(x) for x in objs], [])
    
    @lutils.WorkflowDebugger(flag=FLAGS.debug_workflow)
    def run(self):
        self.params['mode'] = self.mode
        self.params['tprefix'] = lcfg.modes[self.mode]
        processing.processing(self.params)

    @lutils.WorkflowDebugger(flag=FLAGS.debug_workflow)
    def requires(self):
        return DefineBinning()


class HaddHisto(lutils.ForceRun):
    """Write htcondor files for hadd'ing histograms in root files."""
    samples = luigi.ListParameter()
    dataset_name = luigi.Parameter()
    args = utils.dot_dict(lcfg.haddhisto_params)
    args['tprefix'] = lcfg.modes['histos']
    
    @lutils.WorkflowDebugger(flag=FLAGS.debug_workflow)
    def output(self):
        self.args['samples'] = luigi_to_raw( self.samples )
        self.args['dataset_name'] = self.dataset_name
        o1, o2, _, _ = hadd_histo.hadd_histo_outputs( self.args )
        
        target_path = get_target_path( self.__class__.__name__ )
        utils.remove( target_path )
        with open( target_path, 'w' ) as f:
            for t in o1: f.write( t + '\n' )
            for t in o2: f.write( t + '\n' )

        _c1 = convert_to_luigi_local_targets(o1)
        _c2 = convert_to_luigi_local_targets(o2)
        return _c1 + _c2

    @lutils.WorkflowDebugger(flag=FLAGS.debug_workflow)
    def run(self):
        self.args['samples'] = luigi_to_raw( self.samples )
        self.args['dataset_name'] = self.dataset_name
        hadd_histo.hadd_histo( self.args )


class HaddCounts(lutils.ForceRun):
    """Write htcondor files for hadding txt count files."""
    samples = luigi.ListParameter()
    dataset_name = luigi.Parameter()
    args = utils.dot_dict(lcfg.haddcounts_params)
    args['tprefix'] = lcfg.modes['counts']
    
    @lutils.WorkflowDebugger(flag=FLAGS.debug_workflow)
    def output(self):
        self.args['samples'] = luigi_to_raw( self.samples )
        self.args['dataset_name'] = self.dataset_name
        o1, o2, _, _ = hadd_counts.hadd_counts_outputs( self.args )
        
        target_path = get_target_path( self.__class__.__name__ )
        utils.remove( target_path )
        with open( target_path, 'w' ) as f:
            for t in o1: f.write( t + '\n' )
            for t in o2: f.write( t + '\n' )

        _c1 = convert_to_luigi_local_targets(o1)
        _c2 = convert_to_luigi_local_targets(o2)
        return _c1 + _c2

    @lutils.WorkflowDebugger(flag=FLAGS.debug_workflow)
    def run(self):
        self.args['samples'] = luigi_to_raw( self.samples )
        self.args['dataset_name'] = self.dataset_name
        hadd_counts.hadd_counts( self.args )


class EffAndSF(lutils.ForceRun):
    """Write htcondor files for efficiencies and scale factors."""
    params = utils.dot_dict(lcfg.sf_params)
    params['tprefix'] = lcfg.modes['histos']
    
    @lutils.WorkflowDebugger(flag=FLAGS.debug_workflow)
    def output(self):
        o1, o2, _, _ = eff_and_sf.eff_and_sf_outputs(self.params)

        target_path = get_target_path( self.__class__.__name__ )
        utils.remove( target_path )
        with open( target_path, 'w' ) as f:
            f.write( o1 + '\n' )
            f.write( o2 + '\n' )

        _c1 = convert_to_luigi_local_targets(o1)
        _c2 = convert_to_luigi_local_targets(o2)
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
    params = utils.dot_dict(lcfg.sfagg_params)
    
    @lutils.WorkflowDebugger(flag=FLAGS.debug_workflow)
    def output(self):
        o1, o2, _, _ = eff_and_sf_aggr.eff_and_sf_aggr_outputs(self.params)

        target_path = get_target_path( self.__class__.__name__ )
        utils.remove( target_path )
        with open( target_path, 'w' ) as f:
            f.write( o1 + '\n' )
            f.write( o2 + '\n' )

        _c1 = convert_to_luigi_local_targets(o1)
        _c2 = convert_to_luigi_local_targets(o2)
        return _c1 + _c2

    @lutils.WorkflowDebugger(flag=FLAGS.debug_workflow)
    def run(self):
        eff_and_sf_aggr.eff_and_sf_aggr(self.params)


class Discriminator(lutils.ForceRun):
    """Write htcondor files for the variable discriminator."""
    params = utils.dot_dict(lcfg.discriminator_params)
    
    @lutils.WorkflowDebugger(flag=FLAGS.debug_workflow)
    def output(self):
        o1, o2, _, _ = discriminator.discriminator_outputs(self.params)

        target_path = get_target_path( self.__class__.__name__ )
        utils.remove( target_path )
        with open( target_path, 'w' ) as f:
            for t in o1: f.write( t + '\n' )
            for t in o2: f.write( t + '\n' )

        _c1 = convert_to_luigi_local_targets(o1)
        _c2 = convert_to_luigi_local_targets(o2)
        return _c1 + _c2

    @lutils.WorkflowDebugger(flag=FLAGS.debug_workflow)
    def run(self):
        discriminator.discriminator(self.params)


class UnionCalculator(lutils.ForceRun):
    """Write htcondor files for scale factor calculator."""
    params = utils.dot_dict(lcfg.calculator_params)
    
    @lutils.WorkflowDebugger(flag=FLAGS.debug_workflow)
    def output(self):
        o1, o2, _, _ = union_calculator.union_calculator_outputs(self.params)

        #write the target files for debugging
        target_path = get_target_path( self.__class__.__name__ )
        utils.remove( target_path )
        with open( target_path, 'w' ) as f:
            for t in o1: f.write( t + '\n' )
            for t in o2: f.write( t + '\n' )

        _c1 = convert_to_luigi_local_targets(o1)
        _c2 = convert_to_luigi_local_targets(o2)
        return _c1 + _c2

    @lutils.WorkflowDebugger(flag=FLAGS.debug_workflow)
    def run(self):
        union_calculator.union_calculator(self.params)


class Closure(lutils.ForceRun):
    """Write htcondor files for displaying closure plots."""
    params = utils.dot_dict(lcfg.closure_params)
    
    @lutils.WorkflowDebugger(flag=FLAGS.debug_workflow)
    def output(self):
        o1, o2, _, _ = closure.closure_outputs(self.params)

        target_path = get_target_path( self.__class__.__name__ )
        utils.remove( target_path )
        with open( target_path, 'w' ) as f:
            f.write( o1 + '\n' )
            f.write( o2 + '\n' )

        _c1 = convert_to_luigi_local_targets(o1)
        _c2 = convert_to_luigi_local_targets(o2)
        return _c1 + _c2

    @lutils.WorkflowDebugger(flag=FLAGS.debug_workflow)
    def run(self):
        closure.closure(self.params)


class Dag(lutils.ForceRun):
    """Triggering all htcondor writing classes."""
    params        = utils.dot_dict(lcfg.write_params)
    p_histos      = utils.dot_dict(lcfg.histos_params)
    p_hadd_histo  = utils.dot_dict(lcfg.haddhisto_params)
    p_hadd_counts = utils.dot_dict(lcfg.haddhisto_params)
    p_eff_sf      = utils.dot_dict(lcfg.sf_params)
    p_eff_sf_agg  = utils.dot_dict(lcfg.sfagg_params)
    p_disc        = utils.dot_dict(lcfg.discriminator_params)
    p_calc        = utils.dot_dict(lcfg.calculator_params)
    p_closure     = utils.dot_dict(lcfg.closure_params)
    
    p_hadd_histo['tprefix'] = lcfg.modes['histos']
    p_eff_sf['tprefix'] = lcfg.modes['histos']
    
    @lutils.WorkflowDebugger(flag=FLAGS.debug_workflow)
    def output(self):
        o1 = dag.dag_outputs(self.params)

        target_path = get_target_path( self.__class__.__name__ )
        utils.remove( target_path )
        with open( target_path, 'w' ) as f:
            f.write( o1 + '\n' )

        _c1 = convert_to_luigi_local_targets(o1)
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

        self.p_hadd_histo['dataset_name']  = lcfg.data_name
        subm_hadd_hdata = hadd_histo.hadd_histo_outputs(self.p_hadd_histo)[1]
        self.p_hadd_histo['dataset_name']  = lcfg.mc_name
        subm_hadd_hmc = hadd_histo.hadd_histo_outputs(self.p_hadd_histo)[1]

        self.p_hadd_counts['dataset_name'] = lcfg.data_name
        subm_hadd_cdata = hadd_counts.hadd_counts_outputs(self.p_hadd_counts)[1]
        self.p_hadd_counts['dataset_name'] = lcfg.mc_name
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
                                           self.params['data_name'],
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
        new_content = jw.condor_specific_content(queue='short',
                                                 machine='llrt3condor')
        contents.insert(ncontents-1, new_content + '\n')
        with open(out, 'w') as f:
            contents = "".join(contents)
            f.write(contents)
        
    @lutils.WorkflowDebugger(flag=FLAGS.debug_workflow)
    def run(self):
        outfile = self.input()[-1][0].path
        com = 'condor_submit_dag -no_submit -f '
        com += '-outfile_dir {} {}'.format(os.path.dirname(outfile), outfile)

        os.system(com)
        time.sleep(0.5)
        self.edit_condor_submission_file(outfile + '.condor.sub')
        time.sleep(0.5)
        os.system('condor_submit {}.condor.sub'.format(outfile))
        time.sleep(0.5)
        os.system('condor_q')

    @lutils.WorkflowDebugger(flag=FLAGS.debug_workflow)
    def output(self):
        # WriteDag dependency is the last one
        target = self.input()[-1][0].path + '.condor.sub'
        return luigi.LocalTarget(target)

    @lutils.WorkflowDebugger(flag=FLAGS.debug_workflow)
    def requires(self):
        return [ Processing(mode='histos'),
                 Processing(mode='counts'),
                 HaddHisto(dataset_name=lcfg.data_name, samples=lcfg.data_vals),
                 HaddHisto(dataset_name=lcfg.mc_name, samples=lcfg.mc_vals),
                 HaddCounts(dataset_name=lcfg.data_name, samples=lcfg.data_vals),
                 HaddCounts(dataset_name=lcfg.mc_name, samples=lcfg.mc_vals ),
                 EffAndSF(),
                 EffAndSFAggr(),
                 Discriminator(),
                 UnionCalculator(),
                 Closure(),
                 Dag(),
                ]
   
if __name__ == "__main__":
    utils.create_single_dir( lcfg.data_storage )
    utils.create_single_dir( lcfg.targets_folder )
    
    last_tasks = [ SubmitDAG() ]
                    
    if FLAGS.scheduler == 'central':
        luigi.build(last_tasks,
                    workers=FLAGS.workers, local_scheduler=False, log_level='INFO')
    if FLAGS.scheduler == 'local':
        luigi.build(last_tasks,
                    local_scheduler=True, log_level='INFO')

else:
    raise RuntimeError('This script can only be run directly from the command line.')
