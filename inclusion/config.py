import os
import argparse
from argparse import RawTextHelpFormatter
import luigi
from luigi.util import inherits

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
if __name__ == '__main__':
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
