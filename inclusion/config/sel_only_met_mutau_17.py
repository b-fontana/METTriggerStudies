# coding: utf-8

import os
import sys
parent_dir = os.path.abspath(__file__ + 3 * '/..')
sys.path.insert(0, parent_dir)

import inclusion
from inclusion.config import main
from inclusion.utils import utils

bjets_cut = True
category = 'baseline'

custom_cut = ('(self.entries.dau1_pt < 28. and self.entries.dau2_pt < 32.) or ' +
              '(self.entries.dau1_pt < 21. and self.entries.dau2_pt > 32.)')

triggers = ('METNoMu120', 'IsoMu24')
trig_custom = set()
cuts = {'METNoMu120': {'metnomu_et': ('>', [120,]),
                       'mhtnomu_et': ('>', [100,])},
        }

# which triggers are exclusive to a particular channel?
exclusive = {'mutau'   : (),
             'mumu'    : (),
             'general' : ('IsoMu24','METNoMu120'),}

inters_general = {'MET' : (),
                  'EG'  : (),
                  'Mu'  : (('METNoMu120',),
                           ('IsoMu24',),
                           ('IsoMu24', 'METNoMu120'),
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
        ('IsoMu24',)  : 'dau1_pt',
        ('METNoMu120',) : 'dau1_pt',
        ('IsoMu24', 'METNoMu120') : 'dau1_pt',
    },
    'mumu':
    {
        ('IsoMu24',)  : 'dau1_pt',
        ('METNoMu120',) : 'dau1_pt',
        ('IsoMu24', 'METNoMu120') : 'dau1_pt',
    }

}
for x in discr_vars_1D:
    utils.check_discr_vars_correctness(triggers, discr_vars_1D[x], channel=x, exclusive=exclusive)

### Fit efficiencies for some variables
fit_vars = ["metnomu_et",]

### 2D Plots
pairs2D = {'METNoMu120': (('metnomu_et', 'mhtnomu_et'),
                          ('metnomu_et', 'dau1_pt'),
                          ('mhtnomu_et', 'dau1_pt'),
                          ('metnomu_et', 'dau1_eta'),
                          ('mhtnomu_et', 'dau1_eta'),
                          ('dau1_pt', 'dau1_eta'),)}
assert( set(pairs2D.keys()).issubset(set(triggers)) )
for x in pairs2D.values():
    for pair in x:
        assert( pair[0] in main.var_eff and pair[1] in main.var_eff )

### Binning
#pog_pt_binedges = (26., 30., 40., 50., 60., 120., 200)
metnomu_et_binedges = {
    'mutau': (100., 110., 120., 130., 140., 145., 150., 155., 160., 165., 170.,
              175., 180., 185., 190., 195., 200., 210., 220., 230., 240., 250.,
              260., 270., 280.),
    'mumu': (100., 110., 120., 130., 140., 145., 150., 155., 160., 165., 170.,
             175., 180., 185., 190., 195., 200., 210., 220., 230., 240., 250.,
             260., 270., 280.),
}
binedges = {
    # 'dau1_pt': {'etau':   pog_pt_binedges,
    #             'mutau':  pog_pt_binedges,
    #             'tautau': pog_pt_binedges },
    'metnomu_et': {'mutau' : ("quantiles", 100, 300),
                   'mumu'  : ("quantiles", 100, 300) },
    # 'metnomu_et': {'mutau' : metnomu_et_binedges['mutau'],
    #                'mumu'  : metnomu_et_binedges['mumu'] },
    'mhtnomu_et': {'mutau' : (120,360),
                   'mumu'  : (120,360) },
}
assert( set(binedges.keys()).issubset(main.var_join) )
for x in binedges.values():
    assert( set(x.keys()).issubset(main.channels) )
    assert( len(x) == len(list(binedges.values())[0]) )
