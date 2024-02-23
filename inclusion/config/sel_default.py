# coding: utf-8

import os
import sys
parent_dir = os.path.abspath(__file__ + 3 * '/..')
sys.path.insert(0, parent_dir)

import inclusion
from inclusion.config import main
from inclusion.utils import utils

bjets_cut = True
mass_cut = None #standard, inverted
custom_cut = None

triggers = tuple(main.trig_map.keys())
trig_custom = {'VBFTauCustom',
               'IsoDoubleTauCustom',
               'IsoMuIsoTauCustom',
               'EleIsoTauCustom'}
assert all(x in main.trig_map['2018'].keys() for x in trig_custom)
assert(trig_custom.issubset(set(triggers)))

# which triggers are exclusive to a particular channel?
exclusive = {'etau'   : ('Ele32', 'EleIsoTauCustom'),
             'mutau'  : ('IsoMu24', 'IsoMuIsoTauCustom'),
             'tautau' : ('IsoDoubleTauCustom', 'VBFTauCustom'),
             'general': ('METNoMu120', 'IsoTau180')}
for year in {'2016', '2017', '2018'}:
    for excl in exclusive.values():
        assert all(x in main.trig_map[year].keys() for x in excl)

cuts = {#'METNoMu120': {'metnomu_et': ('>', [120,180]), 'mhtnomu_et': ('>', [100,160])},
         #'IsoTau50':   {'dau1_pt': ('>', [80]), 'dau1_eta': ('<', [2.0]), 'met_et': ('>', [150])},
         }
assert(set(cuts.keys()).issubset(set(triggers)) )

inters_general = {'MET' : (),
                  'EG'  : (),
                  'Mu'  : (('IsoTau180',),
                           ('METNoMu120',),
                           ('IsoTau180', 'METNoMu120')),
                  'Tau' : ()
                  }

inters = {
    'etau':
    {'MET' : (),
     'EG'  : (),
     'Mu'  : (('Ele32',),
              ('EleIsoTauCustom',),
              ('Ele32', 'EleIsoTauCustom'),
              ('Ele32', 'METNoMu120'),
              ('Ele32', 'IsoTau180'),
              ('EleIsoTauCustom', 'IsoTau180'),
              ('EleIsoTauCustom', 'METNoMu120'),
              ('Ele32', 'EleIsoTauCustom', 'METNoMu120'),
              ('Ele32', 'IsoTau180', 'METNoMu120'),
              ('Ele32', 'EleIsoTauCustom', 'IsoTau180'),
              ('EleIsoTauCustom', 'IsoTau180', 'METNoMu120'),
              ('Ele32', 'EleIsoTauCustom', 'IsoTau180', 'METNoMu120')),
     'Tau' : ()
     },

    'mutau':
    {'MET' : (),
     'EG'  : (('IsoMu24',),
              ('IsoMuIsoTauCustom',),
              ('IsoMu24', 'IsoMuIsoTauCustom'),
              ('IsoMu24', 'IsoTau180'),
              ('IsoMuIsoTauCustom', 'IsoTau180'),
              ('IsoMu24', 'METNoMu120'),
              ('IsoMuIsoTauCustom', 'METNoMu120'),
              ('IsoMu24', 'IsoTau180', 'METNoMu120'),
              ('IsoMu24', 'IsoMuIsoTauCustom', 'METNoMu120'),
              ('IsoMu24', 'IsoMuIsoTauCustom', 'IsoTau180'),
              ('IsoMuIsoTauCustom', 'IsoTau180', 'METNoMu120'),
              ('IsoMu24', 'IsoMuIsoTauCustom', 'IsoTau180', 'METNoMu120')),
     'Mu'  : (),
     'Tau' : ()
     },

    'tautau':
    {'MET' : (),
     'EG'  : (),
     'Mu'  : (('VBFTauCustom',),
              ('IsoDoubleTauCustom',),
              ('IsoTau180', 'VBFTauCustom'),
              ('METNoMu120', 'VBFTauCustom'),
              ('IsoDoubleTauCustom', 'METNoMu120'),
              ('IsoDoubleTauCustom', 'VBFTauCustom'),
              ('IsoDoubleTauCustom', 'IsoTau180'),
              ('IsoDoubleTauCustom', 'IsoTau180', 'VBFTauCustom'),
              ('IsoTau180', 'METNoMu120', 'VBFTauCustom'),
              ('IsoDoubleTauCustom', 'IsoTau180', 'METNoMu120'),
              ('IsoDoubleTauCustom', 'METNoMu120', 'VBFTauCustom'),
              ('IsoDoubleTauCustom', 'IsoTau180', 'METNoMu120', 'VBFTauCustom')),
     'Tau' : ()
     }
}

for x in inters:
    utils.check_inters_correctness(triggers, inters[x], inters_general, channel=x, exclusive=exclusive)

### 2D Plots
pairs2D = {'METNoMu120': (('metnomu_et', 'mhtnomu_et'),
                          ('dau1_pt', 'dau1_eta'),)}
assert( set(pairs2D.keys()).issubset(set(triggers)) )
for x in pairs2D.values():
    for pair in x:
        assert( pair[0] in main.var_eff and pair[1] in main.var_eff )
