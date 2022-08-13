import os

storage = '/data_CMS/cms/' + os.environ['USER'] + '/TriggerScaleFactors/'

folders = {'base'    : 'METTriggerStudies',
           'main'    : 'inclusion',
           'scripts' : 'scripts',
           'jobs'    : 'condor',
           'subm'    : 'submission',
           'outs'    : 'outputs' }

base_folder = os.path.join(os.environ['HOME'],
                           os.environ['CMSSW_VERSION'], 'src',
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

#_channels = ( 'all', 'etau', 'mutau', 'tautau', 'mumu' )
channels = ('etau', 'mutau', 'tautau')
sel = {'all':    {'pairType': ('<',  3),},
       'mutau':  {'pairType': ('==', 0),},
       'etau':   {'pairType': ('==', 1),},
       'tautau': {'pairType': ('==', 2),},
       'mumu':   {'pairType': ('==', 3),}, # passMu missing for the mumu channel
       'ee':     {'pairType': ('==', 4),} }

# variables considered for calculating and plotting efficiencies
var_eff = ('HT20', 'met_et', 'mht_et', 'metnomu_et', 'mhtnomu_et',
           'dau1_pt', 'dau2_pt', 'dau1_eta', 'dau2_eta')

# variables considered for plotting MC/data comparison distributions
var_dist = ('dau1_pt', 'HH_mass')
# joining the two lists above
var_join = set(var_eff + var_dist)

var_unionweights = ('dau1_pt', 'dau2_pt', 'dau1_eta', 'dau2_eta')

trig_linear = lambda x : {'mc': x, 'data': x}
trig_shift  = lambda x : {'mc': x, 'data': x+5}
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/auTrigger
trig_map = {'IsoMu24': trig_linear(0),
            'Ele32': trig_linear(2),
            'IsoDoubleTauCustom': {'IsoDoubleTau':    {'mc': 4, 'data': (4,5,6)},
                                   'IsoDoubleTauHPS': {'mc': 4, 'data': 7}},
            'IsoMuIsoTauCustom': { 'IsoMuIsoTau':    {'mc': 5, 'data': 9},
                                  'IsoMuIsoTauHPS': {'mc': 5, 'data': 8} },
            'EleIsoTauCustom': {'EleIsoTau': {'mc': 6, 'data': 11},
                                'EleIsoTauHPS': {'mc': 6, 'data': 10}},
            'VBFTauCustom':  {'VBFTau':    {'mc': 8, 'data': 12},
                              'VBFTauHPS': trig_shift(8)},
            'METNoMu120': trig_shift(9),
            'IsoTau180': trig_shift(11)}
triggers = tuple(trig_map.keys())
trig_custom = {'VBFTauCustom',
               'IsoDoubleTauCustom',
               'IsoMuIsoTauCustom',
               'EleIsoTauCustom'}

cuts = {#'METNoMu120': {'metnomu_et': ('>', [120,180]), 'mhtnomu_et': ('>', [100,160])},
         #'IsoTau50':   {'dau1_pt': ('>', [80]), 'dau1_eta': ('<', [2.0]), 'met_et': ('>', [150])},
         }
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
        'tautau': {} }
    
### 2D Plots
pairs2D = {'METNoMu120': (('metnomu_et', 'mhtnomu_et'),
                          ('dau1_pt', 'dau1_eta'),)}
assert( set(pairs2D.keys()).issubset(set(triggers)) )
for x in pairs2D.values():
    for pair in x:
        assert( pair[0] in var_eff and pair[1] in var_eff )

### Binning
pog_pt_binedges = (26., 30., 40., 50., 60., 120., 200)
binedges = {'dau1_pt': {'etau':   pog_pt_binedges,
                        'mutau':  pog_pt_binedges,
                        'tautau': pog_pt_binedges },
            'dau2_pt': {'etau':   pog_pt_binedges,
                        'mutau':  pog_pt_binedges,
                        'tautau': pog_pt_binedges },
            }

### Data and MC samples
inputs = ( '/data_CMS/cms/portales/HHresonant_SKIMS/SKIMS_UL18_220420/', )

# names of the subfolders under '_inputs' above:
# dictionary that maps specific general triggers to datasets 
data = {'MET': ('SKIM_MET',),
        'EG': ('SKIM_EGamma',)
        }

mc_processes = {'ggfRadions': (),
                'ggfBulkGraviton': (),
                'vbfRadion': (),
                'vbfBulkRadion': (),
                'TT': ('SKIM_TT_fullyHad', 'SKIM_TT_fullyLep', 'SKIM_TT_semiLep',),
                'DY': (),
                }


###############################################################################
# Sanity checks
assert len(set(var_unionweights)) == len(var_unionweights)
assert set(var_unionweights).issubset(set(var_eff))
assert(trig_custom.issubset(set(triggers)))
assert(set(cuts.keys()).issubset(set(triggers)) )
for x in cuts.values():
    assert( set(x.keys()).issubset(set(var_eff)) )

assert( set(cuts_ignored.keys()).issubset(var_join) )
for x in cuts_ignored.values():
    assert( set(x).issubset(var_join) )
for k,v in cuts_ignored.items():
    if k in v:
        mes = ''.join(('[configuration, var={}] It is redundant to specify ',
                       'the same variable: cuts are never applied to ',
                       'variables being displayed. Remove it.'.format(k)))
        raise ValueError(mes)
assert( set(binedges.keys()).issubset(var_join) )
for x in binedges.values():
    assert( set(x.keys()).issubset(channels) )
    assert( len(x) == len(list(binedges.values())[0]) )
