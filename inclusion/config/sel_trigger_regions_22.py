# coding: utf-8

import os
import sys
parent_dir = os.path.abspath(__file__ + 3 * '/..')
sys.path.insert(0, parent_dir)

import inclusion
from inclusion.config import main
from inclusion.utils import utils

bjets_cut = None
mass_cut = None #standard, inverted
custom_cut = None
category = 'baseline'

triggers = (
    # 'METNoMu120', 
    'IsoMu24', 
    'IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1',
    'DoubleMediumDeepTauPFTauHPS35_L2NN_eta2p1',
    # "LooseDeepTauPFTauHPS180_L2NN_eta2p1",
    # "PFMETNoMu120_PFMHTNoMu120_IDTight",
    'QuadPFJet70_50_40_30_PFBTagParticleNet_2BTagSum0p65',
    # "Mu50",
    # "Ele28_eta2p1_WPTight_Gsf_HT150"
    # 'IsoTau180'
    )
trig_custom = {
    # 'EleIsoTauCustom', 'IsoMuIsoTauCustom', 'IsoDoubleTauCustom'
    }
cuts = {
    # 'METNoMu120': {'metnomu_et': ('>', [120,]),
    #                    'mhtnomu_et': ('>', [100,])},
        }

# which triggers are exclusive to a particular channel?
exclusive = {'mutau'   : (
                'IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1',),
             'mumu'    : (),
             'tautau': (
                'DoubleMediumDeepTauPFTauHPS35_L2NN_eta2p1',
             ),
             'general' : (
                # "LooseDeepTauPFTauHPS180_L2NN_eta2p1",
                # "PFMETNoMu120_PFMHTNoMu120_IDTight",
                'IsoMu24',
                'QuadPFJet70_50_40_30_PFBTagParticleNet_2BTagSum0p65',
                # "Mu50",
            #  'METNoMu120','IsoTau180'
             ),}

inters_general = {'MET' : (),
                  'EG'  : (),
                  'Mu'  : (
                    # ('METNoMu120',),
                                                       ),
                  'Tau' : ()}
inters = {
    'mutau':
    {'MET' : (),
     'EG'  : (),
     'Mu'  : (('IsoMu24',),
                            ('IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1',),
                            ('IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1', 'IsoMu24',),
                            ('QuadPFJet70_50_40_30_PFBTagParticleNet_2BTagSum0p65',),
                            ('QuadPFJet70_50_40_30_PFBTagParticleNet_2BTagSum0p65', 'IsoMu24'),
                            ('IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1', 'QuadPFJet70_50_40_30_PFBTagParticleNet_2BTagSum0p65',),
                            ('IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1', 'IsoMu24', 'QuadPFJet70_50_40_30_PFBTagParticleNet_2BTagSum0p65',),
),
     'Tau' : ()
     },
    'tautau':
    {'MET' : (),
     'EG'  : (),
     'Mu'  : (),
     'Tau' : (
            ('IsoMu24',),
            ('IsoMu24', 'QuadPFJet70_50_40_30_PFBTagParticleNet_2BTagSum0p65',),
            ('DoubleMediumDeepTauPFTauHPS35_L2NN_eta2p1', 'IsoMu24',),
            ('DoubleMediumDeepTauPFTauHPS35_L2NN_eta2p1', 'IsoMu24', 'QuadPFJet70_50_40_30_PFBTagParticleNet_2BTagSum0p65'),
            ('DoubleMediumDeepTauPFTauHPS35_L2NN_eta2p1', 'QuadPFJet70_50_40_30_PFBTagParticleNet_2BTagSum0p65',),
            ('QuadPFJet70_50_40_30_PFBTagParticleNet_2BTagSum0p65',),
            ('DoubleMediumDeepTauPFTauHPS35_L2NN_eta2p1',),
    )}
}
for x in inters:
    utils.check_inters_correctness(triggers, inters[x], inters_general, channel=x, exclusive=exclusive)


fit_vars = ["dau1_pt", "dau2_pt"]
pairs2D = {'IsoMu24': (('dau1_pt', 'dau2_pt'),), 'IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1': (('dau1_pt', 'dau2_pt'),),
        "DoubleMediumDeepTauPFTauHPS35_L2NN_eta2p1": (('dau1_pt', 'dau2_pt',), ('dau1_tauIdVSjet', 'dau2_tauIdVSjet'),), "QuadPFJet70_50_40_30_PFBTagParticleNet_2BTagSum0p65": (('dau1_pt', 'dau2_pt'), ('bjet1_btagDeepFlavB', 'bjet2_btagDeepFlavB'),), }
assert( set(pairs2D.keys()).issubset(set(triggers)) )
for x in pairs2D.values():
    for pair in x:
        assert( pair[0] in main.var_eff and pair[1] in main.var_eff )

pog_pt_binedges = (0, 26, 30, 40, 50, 60, 120, 200)
binedges = {
    'dau1_pt': {'etau':   pog_pt_binedges,
                'mutau':  pog_pt_binedges,
                'tautau': pog_pt_binedges },
    'dau2_pt': {
                'etau':   pog_pt_binedges,
                'mutau':  pog_pt_binedges,
                'tautau': pog_pt_binedges 
    },
    'dau1_eta': {
        'mutau': (-2.4, -2.1, -1.2, -0.9, 0, 0.9, 1.2, 2.1, 2.4),
        'etau':  (-2.4, -2.1, -1.2, -0.9, 0, 0.9, 1.2, 2.1, 2.4),
        'tautau':  (-2.4, -2.1, -1.2, -0.9, 0, 0.9, 1.2, 2.1, 2.4),
    },
    'dau1_tauIdVSjet': {
        'tautau': (0, 1, 2, 3, 4, 5, 6, 7),
    }, 
    'dau2_tauIdVSjet': {
        'tautau': (0, 1, 2, 3, 4, 5, 6, 7),
    },
    'bjet1_btagDeepFlavB': {
        'tautau': (0, 0.2, 0.3, 0.4, 0.5, 0.6, 1),
    },
    'bjet2_btagDeepFlavB': {
        'tautau': (0, 0.2, 0.3, 0.4, 0.5, 0.6, 1),
    }
    # ("quantiles", 100, 300),
    # 'metnomu_et': {'mutau' : metnomu_et_binedges['mutau'],
    #                'mumu'  : metnomu_et_binedges['mumu'] },
    # 'mhtnomu_et': {'mutau' : (90,360),
    #                'mumu'  : (90,360) },
}