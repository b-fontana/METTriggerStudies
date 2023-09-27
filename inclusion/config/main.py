import os

email = 'bruno.alves@cern.ch'
queue = 'short'

storage = '/data_CMS/cms/' + os.environ['USER'] + '/TriggerScaleFactors/'

folders = {'base'    : 'METTriggerStudies',
           'main'    : 'inclusion',
           'scripts' : 'scripts',
           'jobs'    : 'condor',
           'subm'    : 'submission',
           'outs'    : 'outputs' }

base_folder = os.path.join(os.environ['HOME'],
                           'CMSSW_12_5_0_pre1',#os.environ['CMSSW_VERSION'],
                           'src',
                           folders['base'])
local_folder = os.path.join(base_folder, folders['main'])

targ_def = 'DefaultTarget.txt'
inters_str = '_PLUS_'
nocut_dummy = 'NoCut'
 
pref = {'clos': 'Closure',
        'sf': 'trigSF_',
        'canvas': 'Canvas1D_',
        'histos': 'hist_',
        'counts': 'counts_'}
 
extensions = ('png', 'pdf', 'C', 'root')
placeholder_cuts = '_XXX'

#channels = ( 'all', 'etau', 'mutau', 'tautau', 'mumu' )
channels = ('etau', 'mutau', 'tautau', 'mumu')
sel = {'all'    : {'pairType': ('<',  3),},
       'mutau'  : {'pairType': ('==', 0),},
       'etau'   : {'pairType': ('==', 1),},
       'tautau' : {'pairType': ('==', 2),},
       'mumu'   : {'pairType': ('==', 3),}, # passMu missing for the mumu channel
       'ee'     : {'pairType': ('==', 4),},
       'emu'    : {'pairType': ('==', 5),} }

# variables considered for calculating and plotting efficiencies
var_eff = ('HT20', 'met_et', 'mht_et', 'metnomu_et', 'mhtnomu_et',
           'dau1_pt', 'dau2_pt', 'dau1_eta', 'dau2_eta')

# variables considered for plotting MC/data comparison distributions
var_dist = ('dau1_pt', 'HH_mass')
# joining the two lists above
var_join = set(var_eff + var_dist)

var_unionweights = ('dau1_pt', 'dau2_pt', 'dau1_eta', 'dau2_eta')

trig_linear = lambda x : {'mc': x, 'data': x}

# The bit numbering scheme below matches the 'triggerbit' leaf
# It does NOT match the 'pass_triggerbit' leaf, which is a skimmed version of the above that might change more often
# One way to ensure the scheme is still correct is by running the script as shown here:
# https://github.com/bfonta/useful_scripts/commit/f5e4a0096bc74c89176579a336b0f52b74cb3ed2
trig_map = {'IsoMu24':    {'mc': 0,  'data': 0},
            'Ele32':      {'mc': 2,  'data': 2},
            'METNoMu120': {'mc': 40, 'data': 40},
            'IsoTau180':  {'mc': 4,  'data': 4},
            'IsoDoubleTauCustom': {'IsoDoubleTau':    {'mc': 12, 'data': (13,14,15)},
                                   'IsoDoubleTauHPS': {'mc': 12, 'data': 12}},
            'IsoMuIsoTauCustom':  {'IsoMuIsoTau':     {'mc': 8,  'data': 9},
                                   'IsoMuIsoTauHPS':  {'mc': 8,  'data': 8} },
            'EleIsoTauCustom':    {'EleIsoTau':       {'mc': 10, 'data': 11},
                                   'EleIsoTauHPS':    {'mc': 10, 'data': 10}},
            'VBFTauCustom':       {'VBFTau':          {'mc': 16, 'data': 12},
                                   'VBFTauHPS':       {'mc': 16, 'data': 13}}}
    
lep_triggers = {'Ele32', 'EleIsoTauCustom', 'IsoMu24', 'IsoMuIsoTauCustom',
                'IsoDoubleTauCustom'}
assert all(x in trig_map.keys() for x in lep_triggers)

cuts_ignored = {'HT20':       (),
                'met_et':     ('metnomu_et',),
                'mht_et':     ('mhtnomu_et',),
                'metnomu_et': ('met_et',),
                'mhtnomu_et': ('mht_et',),
                'dau1_pt':    (),
                'dau2_pt':    ()}

### Correlation Matrix
corr = {'etau': {},
        'mutau': {},
        'tautau': {},
        'mumu': {} }
    
### Data and MC samples
inputs = ('/data_CMS/cms/alves/HHresonant_SKIMS/SKIMS_UL18_Sep16_Data/',
          '/data_CMS/cms/alves/HHresonant_SKIMS/SKIMS_UL18_Sep16_Background/',)
corrupted_files = () #'/data_CMS/cms/alves/HHresonant_SKIMS/SKIMS_UL18_EOSv5HighPrio_Background/TTTo2L2Nu/output_9.root'

# names of the subfolders under 'inputs' above:
# dictionary that maps specific general triggers to datasets 
data = {'MET'  : ('MET',),
        'EG'   : ('EGamma',),
        'Mu'   : ('SingleMuon',),
        'Tau'  : ('Tau',)
        }

mc_processes = {'ggfRadions': (),
                'ggfBulkGraviton': (),
                'vbfRadion': (),
                'vbfBulkGraviton': (),
                'TT': ('TTToHadronic', 'TTTo2L2Nu', 'TTToSemiLeptonic',),
                'DY': ('DYJetsToLL_M-10to50', 'DYJetsToLL_M-50_TuneCP5_13TeV-amc',
                       'DYJetsToLL_0J', 'DYJetsToLL_1J', 'DYJetsToLL_2J',
                       'DYJetsToLL_LHEFilterPtZ-0To50', 'DYJetsToLL_LHEFilterPtZ-50To100', 'DYJetsToLL_LHEFilterPtZ-100To250',
                       'DYJetsToLL_LHEFilterPtZ-250To400', 'DYJetsToLL_LHEFilterPtZ-400To650', 'DYJetsToLL_LHEFilterPtZ-650ToInf'),
                'WJets': ('WJetsToLNu_TuneCP5_13TeV-madgraph',
                          'WJetsToLNu_HT-70To100',  'WJetsToLNu_HT-100To200', 'WJetsToLNu_HT-200To400',
                          'WJetsToLNu_HT-400To600', 'WJetsToLNu_HT-600To800', 'WJetsToLNu_HT-800To1200',
                          'WJetsToLNu_HT-1200To2500', 'WJetsToLNu_HT-2500ToInf'),
                }


###############################################################################
# Sanity checks
assert len(set(var_unionweights)) == len(var_unionweights)
assert set(var_unionweights).issubset(set(var_eff))

assert( set(cuts_ignored.keys()).issubset(var_join) )
for x in cuts_ignored.values():
    assert( set(x).issubset(var_join) )
for k,v in cuts_ignored.items():
    if k in v:
        mes = ''.join(('[configuration, var={}] It is redundant to specify ',
                       'the same variable: cuts are never applied to ',
                       'variables being displayed. Remove it.'.format(k)))
        raise ValueError(mes)
