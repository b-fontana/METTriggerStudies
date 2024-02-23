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
category = 'baseline'

triggers = ('METNoMu120', 'IsoMu24', 'IsoTau180')
trig_custom = {'EleIsoTauCustom', 'IsoMuIsoTauCustom', 'IsoDoubleTauCustom'}
cuts = {'METNoMu120': {'metnomu_et': ('>', [120,]),
                       'mhtnomu_et': ('>', [100,])},
        }

# which triggers are exclusive to a particular channel?
exclusive = {'mutau'   : (),
             'mumu'    : (),
             'general' : ('IsoMu24','METNoMu120','IsoTau180'),}

inters_general = {'MET' : (),
                  'EG'  : (),
                  'Mu'  : (('METNoMu120',),
                           ('IsoMu24',),
                           ('IsoTau180',),
                           ('IsoMu24', 'METNoMu120'),
                           ('IsoMu24', 'IsoTau180'),
                           ('METNoMu120', 'IsoTau180'),
                           ('IsoMu24', 'METNoMu120','IsoTau180'),
                           ),
                  'Tau' : ()}
inters = {
    'mutau':
    {'MET' : (),
     'EG'  : (),
     'Mu'  : (),
     'Tau' : ()
     },
    'mumu':
    {'MET' : (),
     'EG'  : (),
     'Mu'  : (),
     'Tau' : ()
     }
}
for x in inters:
    utils.check_inters_correctness(triggers, inters[x], inters_general, channel=x, exclusive=exclusive)
