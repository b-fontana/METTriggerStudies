import os
import time
import utils
import inspect

os.environ['LUIGI_CONFIG_PATH'] = os.environ['PWD']+'/luigi_conf/luigi.cfg'
assert os.path.exists(os.environ['LUIGI_CONFIG_PATH'])
import luigi
from luigi_conf.luigi_utils import  (
    ForceRun,
    WorkflowDebugger,
)
from luigi_conf.luigi_cfg import cfg, FLAGS

from luigi_conf import (
    _placeholder_cuts,
)
lcfg = cfg() #luigi configuration

from scripts.defineBinning import (
    defineBinning,
    defineBinning_outputs,
)
from scripts.writeHTCondorProcessingFiles import (
    writeHTCondorProcessingFiles,
    writeHTCondorProcessingFiles_outputs,
)
from scripts.writeHTCondorHaddHistoFiles import (
    writeHTCondorHaddHistoFiles,
    writeHTCondorHaddHistoFiles_outputs,
)
from scripts.writeHTCondorHaddCountsFiles import (
    writeHTCondorHaddCountsFiles,
    writeHTCondorHaddCountsFiles_outputs,
)
from scripts.writeHTCondorEfficienciesAndScaleFactorsFiles import (
    writeHTCondorEfficienciesAndScaleFactorsFiles,
    writeHTCondorEfficienciesAndScaleFactorsFiles_outputs,
)
from scripts.writeHTCondorEfficienciesAndSFAggregator import (
    writeHTCondorEfficienciesAndSFAggregator,
    writeHTCondorEfficienciesAndSFAggregator_outputs,
)
from scripts.writeHTCondorDiscriminatorFiles import (
    writeHTCondorDiscriminatorFiles,
    writeHTCondorDiscriminatorFiles_outputs,
)
from scripts.writeHTCondorUnionWeightsCalculatorFiles import (
    writeHTCondorUnionWeightsCalculatorFiles,
    writeHTCondorUnionWeightsCalculatorFiles_outputs,
)
# from scripts.writeHTCondorHaddEffFiles import (
#     writeHTCondorHaddEffFiles,
#     writeHTCondorHaddEffFiles_outputs,
# )
from scripts.writeHTCondorClosureFiles import (
    writeHTCondorClosureFiles,
    writeHTCondorClosureFiles_outputs,
)
from scripts.writeHTCondorDAGFiles import (
    WriteDAGManager,
    writeHTCondorDAGFiles_outputs,
)
# from scripts.addTriggerCounts import (
#     addTriggerCounts,
#     addTriggerCounts_outputs
#     )
# from scripts.draw2DTriggerSF import draw2DTriggerSF, draw2DTriggerSF_outputs
# from scripts.drawDistributions import drawDistributions, drawDistributions_outputs

from utils import utils
from scripts.jobWriter import JobWriter

import re
re_txt = re.compile('\.txt')

########################################################################
### HELPER FUNCTIONS ###################################################
########################################################################
def convert_to_luigi_local_targets(targets):
    """Converts a list of files into a list of luigi targets."""
    if not isinstance(targets, (tuple,list,set)):
        targets = [ targets ]
    return [ luigi.LocalTarget(t) for t in targets ]

def get_target_path(taskname):
    target_path = os.path.join(lcfg.targets_folder, re_txt.sub( '_'+taskname+'.txt',
                                                                lcfg.targets_default_name ) ) 
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
    
########################################################################
### CALCULATE THE MOST ADEQUATE BINNING BASED ON DATA ##################
########################################################################
class DefineBinning(luigi.Task):
    args = utils.dot_dict(lcfg.bins_params)
    
    @WorkflowDebugger(flag=FLAGS.debug_workflow)
    def output(self):
        target = defineBinning_outputs( self.args )

        #write the target files for debugging
        target_path = get_target_path( self.__class__.__name__ )
        utils.remove( target_path )
        with open( target_path, 'w' ) as f:
            f.write( target + '\n' )

        return luigi.LocalTarget(target)

    @WorkflowDebugger(flag=FLAGS.debug_workflow)
    def run(self):
        defineBinning( self.args )

########################################################################
### WRITE HTCONDOR FILES FOR TOTAL AND PASSED TRIGGER HISTOGRAMS #######
########################################################################
class WriteHTCondorProcessingFiles(ForceRun):
    params = utils.dot_dict(lcfg.histos_params)

    mode = luigi.ChoiceParameter(choices=lcfg.modes.keys(),
                                 var_type=str)
    
    @WorkflowDebugger(flag=FLAGS.debug_workflow)
    def output(self):
        self.params['mode'] = self.mode
        self.params['tprefix'] = lcfg.modes[self.mode]
        o1, o2, _, _ = writeHTCondorProcessingFiles_outputs(self.params)

        #write the target files for debugging
        target_path = get_target_path( self.__class__.__name__ )
        utils.remove( target_path )
        with open( target_path, 'w' ) as f:
            for t in o1: f.write(t + '\n')
            for t in o2: f.write(t + '\n')

        _c1 = convert_to_luigi_local_targets(o1)
        _c2 = convert_to_luigi_local_targets(o2)
        # for elem in o1:
        #     print(elem, flush=True)
        # print(flush=True)
        # for elem in o2:
        #     print(elem, flush=True)
        # print('----', flush=True)
        # print(flush=True)
        return _c1 + _c2
    
    @WorkflowDebugger(flag=FLAGS.debug_workflow)
    def run(self):
        self.params['mode'] = self.mode
        self.params['tprefix'] = lcfg.modes[self.mode]
        writeHTCondorProcessingFiles(self.params)

    @WorkflowDebugger(flag=FLAGS.debug_workflow)
    def requires(self):
        return DefineBinning()


########################################################################
### WRITE HTCONDOR FILES FOR HADD'ING HISTOGRAMS IN ROOT FILES #########
########################################################################
class WriteHTCondorHaddHistoFiles(ForceRun):
    samples = luigi.ListParameter()
    dataset_name = luigi.Parameter()
    args = utils.dot_dict(lcfg.haddhisto_params)
    args['tprefix'] = lcfg.modes['histos']
    
    @WorkflowDebugger(flag=FLAGS.debug_workflow)
    def output(self):
        self.args['samples'] = luigi_to_raw( self.samples )
        self.args['dataset_name'] = self.dataset_name
        o1, o2, _ = writeHTCondorHaddHistoFiles_outputs( self.args )
        
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
        self.args['samples'] = luigi_to_raw( self.samples )
        self.args['dataset_name'] = self.dataset_name
        writeHTCondorHaddHistoFiles( self.args )

########################################################################
### WRITE HTCONDOR FILES FOR HADDING TXT COUNT FILES ###################
########################################################################
class WriteHTCondorHaddCountsFiles(ForceRun):
    samples = luigi.ListParameter()
    dataset_name = luigi.Parameter()
    args = utils.dot_dict(lcfg.haddcounts_params)
    args['tprefix'] = lcfg.modes['counts']
    
    @WorkflowDebugger(flag=FLAGS.debug_workflow)
    def output(self):
        self.args['samples'] = luigi_to_raw( self.samples )
        self.args['dataset_name'] = self.dataset_name
        o1, o2, _ = writeHTCondorHaddCountsFiles_outputs( self.args )
        
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
        self.args['samples'] = luigi_to_raw( self.samples )
        self.args['dataset_name'] = self.dataset_name
        writeHTCondorHaddCountsFiles( self.args )

########################################################################
### WRITE HTCONDOR FILES FOR EFFICIENCIES AND SCALE FACTORS ############
########################################################################
class WriteHTCondorEfficienciesAndScaleFactorsFiles(ForceRun):
    params = utils.dot_dict(lcfg.sf_params)
    params['tprefix'] = lcfg.modes['histos']
    
    @WorkflowDebugger(flag=FLAGS.debug_workflow)
    def output(self):
        o1, o2, _ = writeHTCondorEfficienciesAndScaleFactorsFiles_outputs(self.params)

        #write the target files for debugging
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
        writeHTCondorEfficienciesAndScaleFactorsFiles(self.params)

########################################################################
### AGGREGATE HTCONDOR FILES FOR EFFICIENCIES AND SCALE FACTORS ########
########################################################################
class WriteHTCondorEfficienciesAndSFAggregator(ForceRun):
    """
    Useful for transfering the intersection efficiencies to the KLUB framework.
    Not needed for the following steps.
    """
    params = utils.dot_dict(lcfg.sfagg_params)
    
    @WorkflowDebugger(flag=FLAGS.debug_workflow)
    def output(self):
        o1, o2, _ = writeHTCondorEfficienciesAndSFAggregator_outputs(self.params)

        #write the target files for debugging
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
        writeHTCondorEfficienciesAndSFAggregator(self.params)

########################################################################
### WRITE HTCONDOR FILES FOR THE VARIABLE DISCRIMINATOR ################
########################################################################
class WriteHTCondorDiscriminatorFiles(ForceRun):
    params = utils.dot_dict(lcfg.discriminator_params)
    
    @WorkflowDebugger(flag=FLAGS.debug_workflow)
    def output(self):
        o1, o2, _ = writeHTCondorDiscriminatorFiles_outputs(self.params)

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
        writeHTCondorDiscriminatorFiles(self.params)

########################################################################
### WRITE HTCONDOR FILES FOR SCALE FACTOR CALCULATOR ###################
########################################################################
class WriteHTCondorUnionWeightsCalculatorFiles(ForceRun):
    params = utils.dot_dict(lcfg.calculator_params)
    
    @WorkflowDebugger(flag=FLAGS.debug_workflow)
    def output(self):
        o1, o2, _ = writeHTCondorUnionWeightsCalculatorFiles_outputs(self.params)

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
        writeHTCondorUnionWeightsCalculatorFiles(self.params)


########################################################################
### WRITE HTCONDOR FILES FOR HADD'ING EFFICIENCIES IN ROOT FILES #######
########################################################################
# class WriteHTCondorHaddEffFiles(ForceRun):
#     samples = luigi.ListParameter()
#     args = utils.dot_dict(lcfg.haddeff_params)
    
#     @WorkflowDebugger(flag=FLAGS.debug_workflow)
#     def output(self):
#         self.args['samples'] = luigi_to_raw( self.samples )
#         o1, o2, _ = writeHTCondorHaddEffFiles_outputs( self.args )
        
#         #write the target files for debugging
#         target_path = get_target_path( self.__class__.__name__ )
#         utils.remove( target_path )
#         with open( target_path, 'w' ) as f:
#             for t in o1: f.write( t + '\n' )
#             for t in o2: f.write( t + '\n' )

#         _c1 = convert_to_luigi_local_targets(o1)
#         _c2 = convert_to_luigi_local_targets(o2)
#         return _c1 + _c2

#     @WorkflowDebugger(flag=FLAGS.debug_workflow)
#     def run(self):
#         self.args['samples'] = luigi_to_raw( self.samples )
#         writeHTCondorHaddEffFiles( self.args )

########################################################################
### WRITE HTCONDOR FILES FOR DISPLAYING CLOSURE PLOTS ##################
########################################################################
class WriteHTCondorClosureFiles(ForceRun):
    params = utils.dot_dict(lcfg.closure_params)
    
    @WorkflowDebugger(flag=FLAGS.debug_workflow)
    def output(self):
        o1, o2, _ = writeHTCondorClosureFiles_outputs(self.params)

        #write the target files for debugging
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
        writeHTCondorClosureFiles(self.params)

########################################################################
### DRAW VARIABLES' DISTRIBUTIONS ######################################
########################################################################
# class DrawDistributions(luigi.Task):
#     args = utils.dot_dict(lcfg.drawcounts_params)
#     args.update( {'tprefix': lcfg.modes['histos'], } )
#     target_path = get_target_path( args.taskname )
    
#     @WorkflowDebugger(flag=FLAGS.debug_workflow)
#     def output(self):
#         targets = []
#         targets_list, _ = drawDistributions_outputs( self.args )

#         #define luigi targets
#         for t in targets_list:
#             targets.append( luigi.LocalTarget(t) )

#         #write the target files for debugging
#         utils.remove( target_path )
#         with open( target_path, 'w' ) as f:
#             for t in targets_list:
#                 f.write( t + '\n' )
                
#         return targets

#     @WorkflowDebugger(flag=FLAGS.debug_workflow)
#     def run(self):
#         drawDistributions( self.args )

#     @WorkflowDebugger(flag=FLAGS.debug_workflow)
#     def requires(self):
#         return [ HaddTriggerEff(samples=self.args.data,
#                                 dataset_name=self.args.data_name),
#                  HaddTriggerEff(samples=self.args.mc_processes,
#                                 dataset_name=self.args.mc_name) ]

########################################################################
### TRIGGERING ALL HTCONDOR WRITING CLASSES ############################
########################################################################
class WriteDAG(ForceRun):
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
        o1 = writeHTCondorDAGFiles_outputs(self.params)

        #write the target files for debugging
        target_path = get_target_path( self.__class__.__name__ )
        utils.remove( target_path )
        with open( target_path, 'w' ) as f:
            f.write( o1 + '\n' )

        _c1 = convert_to_luigi_local_targets(o1)
        return _c1

    @WorkflowDebugger(flag=FLAGS.debug_workflow)
    def run(self):
        self.pHistos['mode'] = 'histos'
        _, submHistos, _, _ = writeHTCondorProcessingFiles_outputs(self.pHistos)
        self.pHistos['mode'] = 'counts'
        _, submCounts, _, _ = writeHTCondorProcessingFiles_outputs(self.pHistos)

        self.pHaddHisto['dataset_name'] = lcfg.data_name
        _, submHaddHistoData, _   = writeHTCondorHaddHistoFiles_outputs(self.pHaddHisto)
        self.pHaddHisto['dataset_name'] = lcfg.mc_name
        _, submHaddHistoMC, _   = writeHTCondorHaddHistoFiles_outputs(self.pHaddHisto)

        self.pHaddCounts['dataset_name'] = lcfg.data_name
        _, submHaddCountsData, _   = writeHTCondorHaddCountsFiles_outputs(self.pHaddCounts)
        self.pHaddCounts['dataset_name'] = lcfg.mc_name
        _, submHaddCountsMC, _   = writeHTCondorHaddCountsFiles_outputs(self.pHaddCounts)

        
        _, submEffSF, _  = writeHTCondorEfficienciesAndScaleFactorsFiles_outputs(self.pEffSF)
        _, submEffSFAgg, _  = writeHTCondorEfficienciesAndSFAggregator_outputs(self.pEffSFAgg)
        _, submDisc, _  = writeHTCondorDiscriminatorFiles_outputs(self.pDisc)
        _, submUnion, _  = writeHTCondorUnionWeightsCalculatorFiles_outputs(self.pSFCalc)
        _, submClosure, _  = writeHTCondorClosureFiles_outputs(self.pClosure)

        jobs = { 'Histos':         submHistos,
                 'Counts':         submCounts,
                 'HaddHistoData':  submHaddHistoData,
                 'HaddHistoMC':    submHaddHistoMC,
                 'HaddCountsData': submHaddCountsData,
                 'HaddCountsMC':   submHaddCountsMC,
                 'EffSF':          [ submEffSF ],
                 'EffSFAgg':       [ submEffSFAgg ],
                 'Discr':          submDisc,
                 'Union':          submUnion,
                 'Closure':        [ submClosure ],
                }

        dag_manager = WriteDAGManager( self.params['localdir'],
                                       self.params['tag'],
                                       self.params['data_name'],
                                       jobs )
        dag_manager.write_all()
        
class SubmitDAG(ForceRun):
    """
    Submission class.
    """
    def edit_condor_submission_file(self, out):
        jw = JobWriter()
        with open(out, 'r') as f:
            contents = f.readlines()
        ncontents = len(contents)
        new_content = jw.llr_condor_specific_content(queue='short')
        contents.insert(ncontents-1, new_content + '\n')
        with open(out, 'w') as f:
            contents = "".join(contents)
            f.write(contents)
        
    @WorkflowDebugger(flag=FLAGS.debug_workflow)
    def run(self):
        out = self.input()[-1][0].path
        os.system('condor_submit_dag -no_submit -f {}'.format(out))
        os.system('sleep 1')
        self.edit_condor_submission_file(out + '.condor.sub')
        os.system('sleep 1')
        os.system('condor_submit {}.condor.sub'.format(out))
        os.system('sleep 1')
        os.system('condor_q')

    @WorkflowDebugger(flag=FLAGS.debug_workflow)
    def output(self):
        # WriteDag dependency is the last one
        target = self.input()[-1][0].path + '.condor.sub'
        return luigi.LocalTarget(target)

    @WorkflowDebugger(flag=FLAGS.debug_workflow)
    def requires(self):
        return [ WriteHTCondorProcessingFiles( mode='histos' ),
                 WriteHTCondorProcessingFiles( mode='counts' ),
                 WriteHTCondorHaddHistoFiles( dataset_name=lcfg.data_name,
                                              samples=lcfg.data_vals ),
                 WriteHTCondorHaddHistoFiles( dataset_name=lcfg.mc_name,
                                              samples=lcfg.mc_vals ),
                 WriteHTCondorHaddCountsFiles( dataset_name=lcfg.data_name,
                                               samples=lcfg.data_vals ),
                 WriteHTCondorHaddCountsFiles( dataset_name=lcfg.mc_name,
                                               samples=lcfg.mc_vals ),
                 WriteHTCondorEfficienciesAndScaleFactorsFiles(),
                 WriteHTCondorEfficienciesAndSFAggregator(),
                 WriteHTCondorDiscriminatorFiles(),
                 WriteHTCondorUnionWeightsCalculatorFiles(),
                 WriteHTCondorClosureFiles(),
                 WriteDAG(),
                ]
   
########################################################################
### MAIN ###############################################################
########################################################################
if __name__ == "__main__":
    utils.create_single_dir( lcfg.data_storage )
    utils.create_single_dir( lcfg.targets_folder )
    
    last_tasks = [ SubmitDAG() ]
        
    # if FLAGS.distributions:
    #     last_tasks += [ DrawDistributions() ]

    # if FLAGS.distributions != 2:
    #     triggercomb = utils.generateTriggerCombinations(FLAGS.triggers)

    #     #one task per trigger combination
    #     for tcomb in triggercomb:
    #         last_tasks += [
    #             Draw1DTriggerScaleFactors(trigger_combination=tcomb)
    #         ]

    # if FLAGS.counts: #overwrites
    #     last_tasks = count_tasks
            
    if FLAGS.scheduler == 'central':
        luigi.build(last_tasks,
                    workers=FLAGS.workers, local_scheduler=False, log_level='INFO')
    if FLAGS.scheduler == 'local':
        luigi.build(last_tasks,
                    local_scheduler=True, log_level='INFO')

else:
    raise RuntimeError('This script can only be run directly from the command line.')
