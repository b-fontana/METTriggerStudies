# coding: utf-8

_all_ = [ ]
    
import os
import sys
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, parent_dir)

import time
import inspect
import re
import argparse
from argparse import RawTextHelpFormatter

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
from luigi.util import inherits
from luigi_conf import luigi_utils as lutils
_extensions = ( 'png',
                'pdf',
                'C',
                'root',
               )
_placeholder_cuts = '_XXX'

_local_home = os.environ['HOME']
_local_cmssw = os.path.join(os.environ['CMSSW_VERSION'], 'src')
_analysis_folders = {'main'    : 'METTriggerStudies/inclusion',
                     'scripts' : 'scripts',
                     'jobs'    : 'condor',
                     'subm'    : 'submission',
                     'outs'    : 'outputs' }

### Channels
#_channels = ( 'all', 'etau', 'mutau', 'tautau', 'mumu' )
_channels = ( 'etau', 'mutau', 'tautau' )
_sel = { 'all':    {'pairType': ('<',  3),},
         'mutau':  {'pairType': ('==', 0),},
         'etau':   {'pairType': ('==', 1),},
         'tautau': {'pairType': ('==', 2),},
         'mumu':   {'pairType': ('==', 3),}, # passMu missing for the mumu channel
         'ee':     {'pairType': ('==', 4),} }

### Variables
# variables considered for calculating and plotting efficiencies
_variables_eff = ['HT20', 'met_et', 'mht_et', 'metnomu_et', 'mhtnomu_et',
                  'dau1_pt', 'dau2_pt', 'dau1_eta', 'dau2_eta']
# variables considered for plotting MC/data comparison distributions
_variables_dist = ['dau1_pt', 'HH_mass']
# joining the two lists above
_variables_join = set(_variables_eff + _variables_dist)

_variables_unionweights = ['dau1_pt', 'dau2_pt', 'dau1_eta', 'dau2_eta']
assert len(set(_variables_unionweights)) == len(_variables_unionweights)
assert set(_variables_unionweights).issubset(set(_variables_eff))

### Triggers
_trigger_linear = lambda x : {'mc': x, 'data': x}
_trigger_shift  = lambda x : {'mc': x, 'data': x+5}
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauTrigger
_triggers_map = {
    'IsoMu24':     _trigger_linear(0),
    #'IsoMu27':     _trigger_linear(1),
    'Ele32':       _trigger_linear(2),
    #'Ele35':       _trigger_linear(3),
    'IsoDoubleTauCustom': {'IsoDoubleTau':    {'mc': 4, 'data': (4,5,6)},
                           'IsoDoubleTauHPS': {'mc': 4, 'data': 7}},
    'IsoMuIsoTauCustom': { 'IsoMuIsoTau':    {'mc': 5, 'data': 9},
                           'IsoMuIsoTauHPS': {'mc': 5, 'data': 8} },
    'EleIsoTauCustom': {'EleIsoTau': {'mc': 6, 'data': 11},
                        'EleIsoTauHPS': {'mc': 6, 'data': 10}},
    'VBFTauCustom':  {'VBFTau':    {'mc': 8, 'data': 12},
                      'VBFTauHPS': _trigger_shift(8)},
    'METNoMu120':  _trigger_shift(9),
    'IsoTau180':   _trigger_shift(11),
}
_triggers_custom = { 'VBFTauCustom',
                     'IsoDoubleTauCustom',
                     'IsoMuIsoTauCustom',
                     'EleIsoTauCustom',
                     }
assert(_triggers_custom.issubset(set(_triggers_map.keys())))

### Cuts
_cuts = {#'METNoMu120': {'metnomu_et': ('>', [120,180]), 'mhtnomu_et': ('>', [100,160])},
         #'IsoTau50':   {'dau1_pt': ('>', [80]), 'dau1_eta': ('<', [2.0]), 'met_et': ('>', [150])},
         }
assert( set(_cuts.keys()).issubset(set(_triggers_map.keys())) )
for x in _cuts.values():
    assert( set(x.keys()).issubset(set(_variables_eff)) )
_cuts_ignored = { 'HT20':       [],
                  'met_et':     ['metnomu_et',],
                  'mht_et':     ['mhtnomu_et',],
                  'metnomu_et': ['met_et',],
                  'mhtnomu_et': ['mht_et',],
                  'dau1_pt':    [],
                  'dau2_pt':    [],
                 }
assert( set(_cuts_ignored.keys()).issubset(_variables_join) )
for x in _cuts_ignored.values():
    assert( set(x).issubset(_variables_join) )
for k,v in _cuts_ignored.items():
    if k in v:
        raise ValueError('[configuration, var={}] It is redundant to specify the same variable: cuts are never applied to variables being displayed. Remove it.'.format(k))

### Correlation Matrix
_corr = {'etau': {},
         'mutau': {},
         'tautau': {} }
    
### 2D Plots
_2Dpairs = {'METNoMu120': (('metnomu_et', 'mhtnomu_et'),
                           ('dau1_pt', 'dau1_eta'),),
            }
assert( set(_2Dpairs.keys()).issubset(set(_triggers_map.keys())) )
for x in _2Dpairs.values():
    for pair in x:
        assert( pair[0] in _variables_eff and pair[1] in _variables_eff )

### Binning
_pog_pt_binedges = [26., 30., 40., 50., 60., 120., 200]
_binedges = { 'dau1_pt': { 'etau':   _pog_pt_binedges,
                           'mutau':  _pog_pt_binedges,
                           'tautau': _pog_pt_binedges },
              'dau2_pt': { 'etau':   _pog_pt_binedges,
                           'mutau':  _pog_pt_binedges,
                           'tautau': _pog_pt_binedges },
             }
assert( set(_binedges.keys()).issubset(_variables_join) )
for x in _binedges.values():
    assert( set(x.keys()).issubset(_channels) )
    assert( len(x) == len(list(_binedges.values())[0]) )

### Data and MC samples
_inputs = [ '/data_CMS/cms/portales/HHresonant_SKIMS/SKIMS_UL18_220420/', ]

# names of the subfolders under '_inputs' above:
# dictionary that maps specific general triggers to datasets 
# both set to MET for framework development phase CHANGE!!!!!!!!!
_data = dict( MET = ['SKIM_MET',],
              EG  = ['SKIM_EGamma',]
              )
#ADD CHECK THAT MAKES SURE THE DATASETS DO NOT REPEAT CHANGE!!!!!!!
_mc_processes = dict( ggfRadions = [],
                      ggfBulkGraviton = [],
                      vbfRadion = [],
                      vbfBulkRadion = [],
                      SingleMuon = [],
                      TT =  ['SKIM_TT_fullyHad',
                             'SKIM_TT_fullyLep',
                             'SKIM_TT_semiLep',],
                      DY = [],
                     )


### Argument parsing
descr = 'Run example: `copython inclusion/run.py --tag abc --data MET EG --mc_process TT`'
parser = argparse.ArgumentParser(description=descr, formatter_class=RawTextHelpFormatter)
parser.add_argument(
    '--nbins',
    type=int,
    default=6,
    help="Number of histogram bins. If fine-grained control is required modify the variable `_bins` in the luigi configuration file."
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
    choices=_data.keys(),
    help='Select the data over which the workflow will be run.'
    )
parser.add_argument(
    '--mc_processes',
    type=str,
    nargs='+',
    required=True,
    choices=_mc_processes.keys(),
    help='Select the MC processes over which the workflow will be run.'
    )
parser.add_argument(
    '--triggers',
    nargs='+', #1 or more arguments
    type=str,
    required=False,
    default=list(_triggers_map.keys()),
    choices=_triggers_map.keys(),
    help='Select the processes over which the workflow will be run.'
    )
parser.add_argument(
    '--triggers_closure',
    nargs='+', #1 or more arguments
    type=str,
    required=False,
    default=list(_triggers_map.keys()),
    choices=_triggers_map.keys(),
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
    default=_channels,
    help='Select the channels over which the workflow will be run.'
    )
parser.add_argument(
    '--variables_for_efficiencies',
    nargs='+', #1 or more arguments
    type=str,
    default=_variables_eff,
    help='Select the variables to be used for the calculation of efficiencies.'
    )
parser.add_argument(
    '--variables_for_distributions',
    nargs='+', #1 or more arguments
    type=str,
    default=_variables_dist,
    help='Select the variables to be used for the display of distributions.'
    )
parser.add_argument(
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
    '--debug_workflow',
    action='store_true',
    help="Explicitly print the functions being run for each task, for workflow debugging purposes."
    )
FLAGS, _ = parser.parse_known_args()
assert set(FLAGS.triggers_closure).issubset(set(FLAGS.triggers))

### Luigi configuration
class cfg(luigi.Config):
     # auxiliary, not used by scripts
    base_name = 'TriggerScaleFactors'
    data_base = os.path.join( '/data_CMS/', 'cms' )
    user = os.environ['USER']
    tag = FLAGS.tag
 
    def flatten_nested_dict(d):
        """
        Splits keys and values.
        Dictionaries are not straightforward to pass as arguments.
        """
        keys, vals = ([] for _ in range(2))
        for k,v in d.items():
            for x in v:
                keys.append(k)
            vals.extend(v)
        return keys, vals
        
    user_data = {k:v for k,v in _data.items() if k in FLAGS.data}
    user_mc   = {k:v for k,v in _mc_processes.items() if k in FLAGS.mc_processes}
    data_keys, data_vals = flatten_nested_dict(user_data)
    assert len(data_keys)==len(data_vals)
    mc_keys, mc_vals = flatten_nested_dict(user_mc)
    assert len(mc_keys)==len(mc_vals)
 
    data_name = 'Data_' + '_'.join(FLAGS.data)
    mc_name   = 'MC_'   + '_'.join(FLAGS.mc_processes)
    
    _storage = os.path.join(data_base, user, base_name, tag)
    data_storage = os.path.join(_storage, 'Data')
    out_storage = os.path.join(_storage, 'Outputs')
 
    # general
    modes = {'histos': 'hist_',
             'counts': 'counts_'}
    _closure_prefix = 'Closure'
    _sf_prefix = 'trigSF_'
 
    subtag = ( FLAGS.subtag if FLAGS.subtag==''
               else ( '_' + FLAGS.subtag if FLAGS.subtag[0] != '_' else FLAGS.subtag ) )
    local_folder = os.path.join(_local_home, _local_cmssw, _analysis_folders['main'])
    targets_folder = os.path.join(data_storage, 'targets')
    targets_default_name = 'DefaultTarget.txt'
    intersection_str = '_PLUS_'
    nocut_dummy_str = 'NoCut'
    
    binedges_filename = os.path.join(data_storage, 'binedges.hdf5')
 
    variables_join = list(set(FLAGS.variables_for_efficiencies + FLAGS.variables_for_distributions))
 
    ####
    #### defineBinning
    ####   
    bins_params = luigi.DictParameter(
        default={ 'nbins'             : FLAGS.nbins,
                  'binedges_filename' : binedges_filename,
                  'indir'             : _inputs,
                  'outdir'            : data_storage,
                  'data_vals'         : data_vals,
                  'variables'         : variables_join,
                  'channels'          : FLAGS.channels,
                  'tag'               : tag,
                  'subtag'            : subtag,
                  'debug'             : FLAGS.debug_workflow} )
 
    ####
    #### dag
    ####
    write_params = { 'data_name' : data_name,
                     'localdir'  : local_folder,
                     'tag'       : tag }
    
    ####
    #### produceTriggerHistograms
    ####
    histos_params = { 'binedges_filename' : binedges_filename,
                      'indir'             : _inputs,
                      'outdir'            : data_storage,
                      'localdir'          : local_folder,
                      'data_keys'         : data_keys,
                      'data_vals'         : data_vals,
                      'mc_keys'           : mc_keys,
                      'mc_vals'           : mc_vals,
                      'triggers'          : FLAGS.triggers,
                      'channels'          : FLAGS.channels,
                      'variables'         : variables_join,
                      'tag'               : tag,
                      'subtag'            : subtag,
                      'intersection_str'  : intersection_str,
                      'nocut_dummy_str'   : nocut_dummy_str,
                      'debug'             : FLAGS.debug_workflow,
                     }
 
    ####
    #### haddHisto
    ####
    haddhisto_params = luigi.DictParameter(
        default={ 'indir'    : data_storage,
                  'localdir' : local_folder,
                  'tag'      : tag,
                  'subtag'   : subtag, } )
 
    ####
    #### haddCounts
    ####
    haddcounts_params = luigi.DictParameter(
        default={ 'indir'    : data_storage,
                  'outdir'   : out_storage,
                  'localdir' : local_folder,
                  'tag'      : tag,
                  'subtag'   : subtag,
                  'channels' : FLAGS.channels, } )
 
    ####
    #### drawTriggerScaleFactors
    ####
    sf_params = luigi.DictParameter(
        default={ 'data_keys'            : data_keys,
                  'mc_keys'              : mc_keys,
                  'data_vals'            : data_vals,
                  'mc_vals'              : mc_vals,
                  'draw_independent_MCs' : False,
                  'indir'                : data_storage,
                  'outdir'               : out_storage,
                  'localdir'             : local_folder,
                  'triggers'             : FLAGS.triggers,
                  'channels'             : FLAGS.channels,
                  'variables'            : FLAGS.variables_for_efficiencies,
                  'binedges_filename'    : binedges_filename,
                  'tag'                  : tag,
                  'subtag'               : subtag,
                  'canvas_prefix'        : 'Canvas1D_',
                  'intersection_str'     : intersection_str,
                  'nocut_dummy_str'      : nocut_dummy_str,
                  'debug'                : FLAGS.debug_workflow,} )
 
    sfagg_params = luigi.DictParameter(
        default={ 'indir'       : out_storage,
                  'outdir'      : out_storage,
                  'localdir'    : local_folder,
                  'channels'    : FLAGS.channels,
                  'variables'   : FLAGS.variables_for_efficiencies,
                  'tag'         : tag,
                  'subtag'      : subtag,
                  'file_prefix' : _sf_prefix,
                  'debug'       : FLAGS.debug_workflow,} )
 
    ####
    #### variableImportanceDiscriminator
    ####
    discriminator_params = luigi.DictParameter(
        default={ 'indir'            : data_storage,
                  'outdir'           : data_storage,
                  'localdir'         : local_folder,
                  'triggers'         : FLAGS.triggers,
                  'channels'         : FLAGS.channels,
                  'variables'        : FLAGS.variables_for_efficiencies,
                  'tag'              : tag,
                  'subtag'           : subtag,
                  'intersection_str' : intersection_str,
                  'debug'            : FLAGS.debug_workflow,} )
 
    ####
    #### scale factor calculator
    ####
    calculator_params = luigi.DictParameter(
        default={ 'binedges_filename'       : binedges_filename,
                  'indir_root'              : _inputs,
                  'indir_json'              : data_storage,
                  'indir_eff'               : out_storage,
                  'outdir'                  : data_storage,
                  'outprefix'               : _closure_prefix,
                  'data_name'               : data_name,
                  'mc_name'                 : mc_name,
                  'mc_processes'            : mc_vals,
                  'localdir'                : local_folder,
                  'triggers'                : FLAGS.triggers,
                  'closure_single_triggers' : FLAGS.triggers_closure,
                  'channels'                : FLAGS.channels,
                  'variables'               : FLAGS.variables_for_efficiencies,
                  'tag'                     : tag,
                  'subtag'                  : subtag,
                  'debug'                   : FLAGS.debug_workflow,} )
 
    ####
    #### draw single efficiencies closure
    ####
    closure_params = luigi.DictParameter(
        default={ 'data_name'               : data_name,
                  'mc_name'                 : mc_name,
                  'binedges_filename'       : binedges_filename,
                  'indir_eff'               : out_storage,
                  'indir_union'             : data_storage,
                  'indir_json'              : data_storage,
                  'outdir'                  : out_storage,
                  'inprefix'                : _closure_prefix,
                  'eff_prefix'              : _sf_prefix,
                  'mc_processes'            : mc_vals,
                  'out_weighted_prefix'     : _closure_prefix,
                  'out_original_prefix'     : modes['histos'],
                  'localdir'                : local_folder,
                  'closure_single_triggers' : FLAGS.triggers_closure,
                  'channels'                : FLAGS.channels,
                  'variables'               : FLAGS.variables_for_efficiencies,
                  'tag'                     : tag,
                  'subtag'                  : subtag,
                  'debug'                   : FLAGS.debug_workflow } )

lcfg = cfg()
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
        new_content = jw.condor_specific_content(queue='long',
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
        time.sleep(1.5)
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
