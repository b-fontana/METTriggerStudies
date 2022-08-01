import os
import time
import inspect
import re

import luigi
from luigi_conf import luigi_utils as lutils
from luigi_conf.luigi_cfg import cfg, FLAGS
lcfg = cfg() #luigi configuration

from utils import utils
from scripts import defineBinning
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
    
    @WorkflowDebugger(flag=FLAGS.debug_workflow)
    def output(self):
        target = defineBinning.defineBinning_outputs( self.args )

        #write the target files for debugging
        target_path = get_target_path( self.__class__.__name__ )
        utils.remove( target_path )
        with open( target_path, 'w' ) as f:
            f.write( target + '\n' )

        return luigi.LocalTarget(target)

    @WorkflowDebugger(flag=FLAGS.debug_workflow)
    def run(self):
        defineBinning.defineBinning( self.args )


class Processing(ForceRun):
    """Write htcondor files for total and passed trigger histograms."""
    params = utils.dot_dict(lcfg.histos_params)

    mode = luigi.ChoiceParameter(choices=lcfg.modes.keys(),
                                 var_type=str)
    
    @WorkflowDebugger(flag=FLAGS.debug_workflow)
    def output(self):
        self.params['mode'] = self.mode
        self.params['tprefix'] = lcfg.modes[self.mode]
        objData, objMC, _, _ = processing.processing_outputs(self.params)
        o1, o2, _ = objData
        o3, o4, _ = objMC

        target_path = get_target_path( self.__class__.__name__ )
        utils.remove( target_path )
        with open( target_path, 'w' ) as f:
            for t in o1: f.write(t + '\n')
            for t in o2: f.write(t + '\n')
            for t in o3: f.write(t + '\n')
            for t in o4: f.write(t + '\n')

        _c1 = convert_to_luigi_local_targets(o1)
        _c2 = convert_to_luigi_local_targets(o2)
        return _c1 + _c2
    
    @WorkflowDebugger(flag=FLAGS.debug_workflow)
    def run(self):
        self.params['mode'] = self.mode
        self.params['tprefix'] = lcfg.modes[self.mode]
        processing.processing(self.params)

    @WorkflowDebugger(flag=FLAGS.debug_workflow)
    def requires(self):
        return DefineBinning()


class HaddHisto(ForceRun):
    """Write htcondor files for hadd'ing histograms in root files."""
    samples = luigi.ListParameter()
    dataset_name = luigi.Parameter()
    args = utils.dot_dict(lcfg.haddhisto_params)
    args['tprefix'] = lcfg.modes['histos']
    
    @WorkflowDebugger(flag=FLAGS.debug_workflow)
    def output(self):
        self.args['samples'] = luigi_to_raw( self.samples )
        self.args['dataset_name'] = self.dataset_name
        o1, o2, _ = hadd_histo.hadd_histo_outputs( self.args )
        
        target_path = get_target_path( self.__class__.__name__ )
        utils.remove( target_path )
        with open( target_path, 'w' ) as f:
            for t in o1: f.write( t + '\n' )
            for t in o2: f.write( t + '\n' )

        _c1 = convert_to_luigi_local_targets(o1)
        _c2 = convert_to_luigi_local_targets(o2)
        return _c1 + _c2

    @WorkflowDebugger(flag=FLAGS.debug_workflow)
    def run(self):
        self.args['samples'] = luigi_to_raw( self.samples )
        self.args['dataset_name'] = self.dataset_name
        hadd_histo.hadd_histo( self.args )


class HaddCounts(ForceRun):
    """Write htcondor files for hadding txt count files."""
    samples = luigi.ListParameter()
    dataset_name = luigi.Parameter()
    args = utils.dot_dict(lcfg.haddcounts_params)
    args['tprefix'] = lcfg.modes['counts']
    
    @WorkflowDebugger(flag=FLAGS.debug_workflow)
    def output(self):
        self.args['samples'] = luigi_to_raw( self.samples )
        self.args['dataset_name'] = self.dataset_name
        o1, o2, _ = hadd_counts.hadd_counts_outputs( self.args )
        
        target_path = get_target_path( self.__class__.__name__ )
        utils.remove( target_path )
        with open( target_path, 'w' ) as f:
            for t in o1: f.write( t + '\n' )
            for t in o2: f.write( t + '\n' )

        _c1 = convert_to_luigi_local_targets(o1)
        _c2 = convert_to_luigi_local_targets(o2)
        return _c1 + _c2

    @WorkflowDebugger(flag=FLAGS.debug_workflow)
    def run(self):
        self.args['samples'] = luigi_to_raw( self.samples )
        self.args['dataset_name'] = self.dataset_name
        hadd_counts.hadd_counts( self.args )


class EffAndSF(ForceRun):
    """Write htcondor files for efficiencies and scale factors."""
    params = utils.dot_dict(lcfg.sf_params)
    params['tprefix'] = lcfg.modes['histos']
    
    @WorkflowDebugger(flag=FLAGS.debug_workflow)
    def output(self):
        o1, o2, _ = eff_and_sf.eff_and_sf_outputs(self.params)

        target_path = get_target_path( self.__class__.__name__ )
        utils.remove( target_path )
        with open( target_path, 'w' ) as f:
            f.write( o1 + '\n' )
            f.write( o2 + '\n' )

        _c1 = convert_to_luigi_local_targets(o1)
        _c2 = convert_to_luigi_local_targets(o2)
        return _c1 + _c2

    @WorkflowDebugger(flag=FLAGS.debug_workflow)
    def run(self):
        eff_and_sf.eff_and_sf(self.params)

class EffAndSFAggr(ForceRun):
    """
    Aggregate htcondor files for efficiencies and scale factors.
    Useful for transfering the intersection efficiencies to the KLUB framework.
    Not needed for the following steps.
    """
    params = utils.dot_dict(lcfg.sfagg_params)
    
    @WorkflowDebugger(flag=FLAGS.debug_workflow)
    def output(self):
        o1, o2, _ = eff_and_sf_aggr.eff_and_sf_aggr_outputs(self.params)

        target_path = get_target_path( self.__class__.__name__ )
        utils.remove( target_path )
        with open( target_path, 'w' ) as f:
            f.write( o1 + '\n' )
            f.write( o2 + '\n' )

        _c1 = convert_to_luigi_local_targets(o1)
        _c2 = convert_to_luigi_local_targets(o2)
        return _c1 + _c2

    @WorkflowDebugger(flag=FLAGS.debug_workflow)
    def run(self):
        eff_and_sf_aggr.eff_and_sf_aggr(self.params)


class Discriminator(ForceRun):
    """Write htcondor files for the variable discriminator."""
    params = utils.dot_dict(lcfg.discriminator_params)
    
    @WorkflowDebugger(flag=FLAGS.debug_workflow)
    def output(self):
        o1, o2, _ = discriminator.discriminator_outputs(self.params)

        target_path = get_target_path( self.__class__.__name__ )
        utils.remove( target_path )
        with open( target_path, 'w' ) as f:
            for t in o1: f.write( t + '\n' )
            for t in o2: f.write( t + '\n' )

        _c1 = convert_to_luigi_local_targets(o1)
        _c2 = convert_to_luigi_local_targets(o2)
        return _c1 + _c2

    @WorkflowDebugger(flag=FLAGS.debug_workflow)
    def run(self):
        discriminator.discriminator(self.params)


class UnionCalculator(ForceRun):
    """Write htcondor files for scale factor calculator."""
    params = utils.dot_dict(lcfg.calculator_params)
    
    @WorkflowDebugger(flag=FLAGS.debug_workflow)
    def output(self):
        o1, o2, _ = union_calculator.union_calculator_outputs(self.params)

        #write the target files for debugging
        target_path = get_target_path( self.__class__.__name__ )
        utils.remove( target_path )
        with open( target_path, 'w' ) as f:
            for t in o1: f.write( t + '\n' )
            for t in o2: f.write( t + '\n' )

        _c1 = convert_to_luigi_local_targets(o1)
        _c2 = convert_to_luigi_local_targets(o2)
        return _c1 + _c2

    @WorkflowDebugger(flag=FLAGS.debug_workflow)
    def run(self):
        union_calculator.union_calculator(self.params)


class Closure(ForceRun):
    """Write htcondor files for displaying closure plots."""
    params = utils.dot_dict(lcfg.closure_params)
    
    @WorkflowDebugger(flag=FLAGS.debug_workflow)
    def output(self):
        o1, o2, _ = closure.closure_outputs(self.params)

        target_path = get_target_path( self.__class__.__name__ )
        utils.remove( target_path )
        with open( target_path, 'w' ) as f:
            f.write( o1 + '\n' )
            f.write( o2 + '\n' )

        _c1 = convert_to_luigi_local_targets(o1)
        _c2 = convert_to_luigi_local_targets(o2)
        return _c1 + _c2

    @WorkflowDebugger(flag=FLAGS.debug_workflow)
    def run(self):
        closure.closure(self.params)


class Dag(ForceRun):
    """Triggering all htcondor writing classes."""
    params      = utils.dot_dict(lcfg.write_params)
    pHistos     = utils.dot_dict(lcfg.histos_params)
    pHaddHisto  = utils.dot_dict(lcfg.haddhisto_params)
    pHaddCounts = utils.dot_dict(lcfg.haddhisto_params)
    pEffSF      = utils.dot_dict(lcfg.sf_params)
    pEffSFAgg   = utils.dot_dict(lcfg.sfagg_params)
    pDisc       = utils.dot_dict(lcfg.discriminator_params)
    pSFCalc     = utils.dot_dict(lcfg.calculator_params)
    pClosure    = utils.dot_dict(lcfg.closure_params)
    
    pHaddHisto['tprefix']  = lcfg.modes['histos']
    pEffSF['tprefix'] =  lcfg.modes['histos']
    
    @WorkflowDebugger(flag=FLAGS.debug_workflow)
    def output(self):
        o1 = dag.dag_outputs(self.params)

        target_path = get_target_path( self.__class__.__name__ )
        utils.remove( target_path )
        with open( target_path, 'w' ) as f:
            f.write( o1 + '\n' )

        _c1 = convert_to_luigi_local_targets(o1)
        return _c1

    @WorkflowDebugger(flag=FLAGS.debug_workflow)
    def run(self):
        self.pHistos['mode']             = 'histos'
        objHistosData, objHistosMC, _, _ = processing.processing_outputs(self.pHistos)
        _, submHistosData, _             = objHistosData
        _, submHistosMC, _               = objHistosMC

        self.pHistos['mode']             = 'counts'
        objCountsData, objCountsMC, _, _ = processing.processing_outputs(self.pHistos)
        _, submCountsData, _             = objCountsData
        _, submCountsMC, _               = objCountsMC

        self.pHaddHisto['dataset_name']  = lcfg.data_name
        _, submHaddHistoData, _          = hadd_histo.hadd_histo_outputs(self.pHaddHisto)
        self.pHaddHisto['dataset_name']  = lcfg.mc_name
        _, submHaddHistoMC, _            = hadd_histo.hadd_histo_outputs(self.pHaddHisto)

        self.pHaddCounts['dataset_name'] = lcfg.data_name
        _, submHaddCountsData, _         = hadd_counts.hadd_counts_outputs(self.pHaddCounts)
        self.pHaddCounts['dataset_name'] = lcfg.mc_name
        _, submHaddCountsMC, _           = hadd_counts.hadd_counts_outputs(self.pHaddCounts)
        
        _, submEffSF, _                  = eff_and_sf.eff_and_sf_outputs(self.pEffSF)
        _, submEffSFAgg, _               = eff_and_sf_aggr.eff_and_sf_aggr_outputs(self.pEffSFAgg)
        _, submDisc, _                   = discriminator.discriminator_outputs(self.pDisc)
        # _, submUnion, _                = union_calculator.union_calculator_outputs(self.pSFCalc)
        # _, submClosure, _              = closure.closure_outputs(self.pClosure)

        jobs = { 'HistosData':     submHistosData,
                 'HistosMC':       submHistosMC,
                 'CountsData':     submCountsData,
                 'CountsMC':       submCountsMC,
                 'HaddHistoData':  submHaddHistoData,
                 'HaddHistoMC':    submHaddHistoMC,
                 'HaddCountsData': submHaddCountsData,
                 'HaddCountsMC':   submHaddCountsMC,
                 'EffSF':          [ submEffSF ],
                 'EffSFAgg':       [ submEffSFAgg ],
                 'Discr':          submDisc,
                 #'Union':          submUnion,
                 #'Closure':        [ submClosure ],
                }

        dag_manager = dag.WriteDAGManager( self.params['localdir'],
                                           self.params['tag'],
                                           self.params['data_name'],
                                           jobs,
                                           mode='short' )
        dag_manager.write_all()
        
class SubmitDAG(ForceRun):
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
        
    @WorkflowDebugger(flag=FLAGS.debug_workflow)
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

    @WorkflowDebugger(flag=FLAGS.debug_workflow)
    def output(self):
        # WriteDag dependency is the last one
        target = self.input()[-1][0].path + '.condor.sub'
        return luigi.LocalTarget(target)

    @WorkflowDebugger(flag=FLAGS.debug_workflow)
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
