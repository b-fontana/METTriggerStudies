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

triggers = ('METNoMu120', 'IsoMu27')
trig_custom = set()
cuts = {'METNoMu120': {'metnomu_et': ('>', [0.,]),
                       'mhtnomu_et': ('>', [0.,])},
        }

# which triggers are exclusive to a particular channel?
exclusive = {'mutau'   : (),
             'mumu'    : (),
             'general' : ('IsoMu27','METNoMu120'),}

inters_general = {'MET' : (),
                  'EG'  : (),
                  'Mu'  : (('METNoMu120',),
                           ('IsoMu27',),
                           ('IsoMu27', 'METNoMu120'),
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

discr_vars_1D =  {
    'mutau':
    {
        ('IsoMu27',)  : 'dau1_pt',
        ('METNoMu120',) : 'dau1_pt',
        ('IsoMu27', 'METNoMu120') : 'dau1_pt',
    },
    'mumu':
    {
        ('IsoMu27',)  : 'dau1_pt',
        ('METNoMu120',) : 'dau1_pt',
        ('IsoMu27', 'METNoMu120') : 'dau1_pt',
    }

}
for x in discr_vars_1D:
    utils.check_discr_vars_correctness(triggers, discr_vars_1D[x], channel=x, exclusive=exclusive)

### Fit efficiencies for some variables
fit_vars = ["metnomu_et",]

### 2D Plots
pairs2D = {'METNoMu120': (('metnomu_et', 'mhtnomu_et'),)}
assert( set(pairs2D.keys()).issubset(set(triggers)) )
for x in pairs2D.values():
    for pair in x:
        assert( pair[0] in main.var_eff and pair[1] in main.var_eff )

### Binning
#pog_pt_binedges = (26., 30., 40., 50., 60., 120., 200)
metnomu_et_binedges = {
    'mutau': (70., 95., 115., 130., 140., 150., 160., 170., 180., 190., 200., 210., 222.5, 235., 250., 265., 280., 300., 325.),
    'mumu':  (70., 95., 115., 130., 140., 150., 160., 170., 180., 190., 200., 210., 222.5, 235., 250., 265., 280., 300., 325.),
}
binedges = {
    # 'dau1_pt': {'etau':   pog_pt_binedges,
    #             'mutau':  pog_pt_binedges,
    #             'tautau': pog_pt_binedges },
    #("quantiles", 100, 300),
    'metnomu_et': {'mutau' : metnomu_et_binedges['mutau'],
                   'mumu'  : metnomu_et_binedges['mumu'] },
    'mhtnomu_et': {'mutau' : (100,360),
                   'mumu'  : (100,360) },
}
assert( set(binedges.keys()).issubset(main.var_join) )
for x in binedges.values():
    assert( set(x.keys()).issubset(main.channels) )
    assert( len(x) == len(list(binedges.values())[0]) )
