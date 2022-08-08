import os
import argparse
from argparse import RawTextHelpFormatter
import luigi
from luigi.util import inherits

from . import _inputs, _data, _mc_processes, _triggers_map, _channels
from . import _variables_eff, _variables_dist
from . import _trigger_shift, _triggers_map

######################################################################## 
### ARGUMENT PARSING ###################################################
########################################################################
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

########################################################################
### LUIGI CONFIGURATION ################################################
########################################################################
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

    local_home = os.environ['HOME']
    local_cmssw = os.path.join(os.environ['CMSSW_VERSION'], 'src')
    analysis_folders = {'main'    : 'METTriggerStudies',
                        'scripts' : 'scripts',
                        'jobs'    : 'condor'}

    # general
    modes = {'histos': 'hist_',
             'counts': 'counts_'}
    _closure_prefix = 'Closure'
    _sf_prefix = 'trigSF_'

    subtag = ( FLAGS.subtag if FLAGS.subtag==''
               else ( '_' + FLAGS.subtag if FLAGS.subtag[0] != '_' else FLAGS.subtag ) )
    local_folder = os.path.join(local_home, local_cmssw, analysis_folders['main'])
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
