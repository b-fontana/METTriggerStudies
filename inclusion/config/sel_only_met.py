# coding: utf-8

import os
import sys
parent_dir = os.path.abspath(__file__ + 3 * '/..')
sys.path.insert(0, parent_dir)

import inclusion
from inclusion.config import main
from inclusion.utils import utils

bjets_cut = False

inters_general = {'MET' : (),
                  'EG'  : (),
                  'Mu'  : (('METNoMu120',),
                           ),
                  'Tau' : ()
                  }
inters = {
    'etau':
    {'MET' : (),
     'EG'  : (),
     'Mu'  : (),
     'Tau' : ()
     },
    'mutau':
    {'MET' : (),
     'EG'  : (),
     'Mu'  : (),
     'Tau' : ()
     },
    'tautau':
    {'MET' : (),
     'EG'  : (),
     'Mu'  : (),
     'Tau' : ()
     }
}
for x in inters:
    utils.check_inters_correctness(inters[x], inters_general, channel=x)

discr_vars_1D =  {
    'etau':
    {
        ('METNoMu120',) : 'dau1_pt',
    },

    'mutau':
    {
        ('METNoMu120',) : 'dau1_pt',
    },

    'tautau':
    {
        ('METNoMu120',) : 'dau1_pt',
    }

}
for x in discr_vars_1D:
    utils.check_discr_vars_correctness(discr_vars_1D[x], channel=x)
