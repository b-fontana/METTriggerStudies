"""
Configuration file for the Luigi trigger scale factors framework.
Sanity checks included.
"""
_extensions = ( 'png',
                'pdf',
                'C',
                'root',
               )
_placeholder_cuts = '_XXX'

#######################################################################################################
########### CHANNELS ##################################################################################
#######################################################################################################
#_channels = ( 'all', 'etau', 'mutau', 'tautau', 'mumu' )
_channels = ( 'etau', 'mutau', 'tautau' )
_sel = { 'all':    {'pairType': ('<',  3),},
         'mutau':  {'pairType': ('==', 0),},
         'etau':   {'pairType': ('==', 1),},
         'tautau': {'pairType': ('==', 2),},
         'mumu':   {'pairType': ('==', 3),}, # passMu missing for the mumu channel
         'ee':     {'pairType': ('==', 4),} }

#######################################################################################################
########### VARIABLES #################################################################################
#######################################################################################################
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

#######################################################################################################
########### TRIGGERS ##################################################################################
#######################################################################################################
_trigger_linear = lambda x : {'mc': x, 'data': x}
_trigger_shift  = lambda x : {'mc': x, 'data': x+5}
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauTrigger
_triggers_map = {
    'IsoMu24':     _trigger_linear(0),
    'IsoMu27':     _trigger_linear(1),
    'Ele32':       _trigger_linear(2),
    'Ele35':       _trigger_linear(3),
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

#######################################################################################################
########### CUTS ######################################################################################
#######################################################################################################
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

#######################################################################################################
########### CORRELATION MATRIX ########################################################################
#######################################################################################################
_corr = {'etau': {},
         'mutau': {},
         'tautau': {}
         }
    
#######################################################################################################
########### 2D PLOTS ##################################################################################
#######################################################################################################
_2Dpairs = {#'METNoMu120':      (('metnomu_et', 'mhtnomu_et'),),
            }
assert( set(_2Dpairs.keys()).issubset(set(_triggers_map.keys())) )
for x in _2Dpairs.values():
    for pair in x:
        assert( pair[0] in _variables_eff and pair[1] in _variables_eff )

#######################################################################################################
########### BINNING ###################################################################################
#######################################################################################################
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

#######################################################################################################
########### DATA AND MC SAMPLES #######################################################################
#######################################################################################################
_inputs = [ '/data_CMS/cms/portales/HHresonant_SKIMS/SKIMS_UL18_220420/', ] #data, MC signal and MC backgrounds
# _inputs = [ '/data_CMS/cms/portales/HHresonant_SKIMS/SKIMS_UL18_220325_MiniAODv2/',
#             '/data_CMS/cms/portales/HHresonant_SKIMS/SKIMS_UL18_220331_data_MiniAODv2/' ] 

# names of the subfolders under '_inputs' above
_data = dict( MET = ['SKIM_MET',] )
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

"""
'pass_triggerbit' leaf

data:
0 - HLT_IsoMu24_v
1 - HLT_IsoMu27_v
2 - HLT_Ele32_WPTight_Gsf_v
3 - HLT_Ele35_WPTight_Gsf_v
4 - HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_v
5 - HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_v
6 - HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg_v
7 - HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_v
8 - HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1_v
9 - HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1_v
10 - HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_CrossL1_v
11 - HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1_v
12 - HLT_VBF_DoubleLooseChargedIsoPFTau20_Trk1_eta2p1_v
13 - HLT_VBF_DoubleLooseChargedIsoPFTauHPS20_Trk1_eta2p1_v
14 - HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v
15 - HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_v
16 - HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr_v

MC:
0 - HLT_IsoMu24_v
1 - HLT_IsoMu27_v
2 - HLT_Ele32_WPTight_Gsf_v
3 - HLT_Ele35_WPTight_Gsf_v
4 - HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_v
5 - HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1_v
6 - HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_CrossL1_v
7 - HLT_VBF_DoubleLooseChargedIsoPFTau20_Trk1_eta2p1_v
8 - HLT_VBF_DoubleLooseChargedIsoPFTauHPS20_Trk1_eta2p1_v
9 - HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_vgithgith
10 - HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100_v
11 - HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_v

Leg splitting:

4: one tau, one tau
5: one muon, one tau
6: one electron, one tau
7/8: remove (dedicated space region)
9: MET, MHT
10: remove
11: one tau






















'triggerbit' leaf

data and MC:
0 - HLT_IsoMu24_v
1 - HLT_IsoMu27_v
2 - HLT_Ele32_WPTight_Gsf_v
3 - HLT_Ele35_WPTight_Gsf_v
4 - HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_v
5 - HLT_MediumChargedIsoPFTau200HighPtRelaxedIso_Trk50_eta2p1_v
6 - HLT_MediumChargedIsoPFTau220HighPtRelaxedIso_Trk50_eta2p1_v
7 - HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr_v
8 - HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1_v
9 - HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1_v
10 - HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_CrossL1_v
11 - HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1_v
12 - HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_v
13 - HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_v
14 - HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_v
15 - HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg_v
16 - HLT_VBF_DoubleLooseChargedIsoPFTauHPS20_Trk1_eta2p1_v
17 - HLT_VBF_DoubleLooseChargedIsoPFTau20_Trk1_eta2p1_v
18 - HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET90_v
19 - HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100_v
20 - HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET110_v
21 - HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET130_v
22 - HLT_Ele28_eta2p1_WPTight_Gsf_HT150_v
23 - HLT_Ele32_WPTight_Gsf_L1DoubleEG_v
24 - HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto_v
25 - HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165_v
26 - HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5_v
27 - HLT_Mu50_v
28 - HLT_TkMu100_v
29 - HLT_OldMu100_v
30 - HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight_v
31 - HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v
32 - HLT_Mu17_Photon30_IsoCaloId_v
33 - HLT_DoubleMu4_Mass3p8_DZ_PFHT350_v
34 - HLT_DoubleMu3_DCA_PFMET50_PFMHT60_v
35 - HLT_AK8PFJet330_TrimMass30_PFAK8BoostedDoubleB_np4_v
36 - HLT_QuadPFJet103_88_75_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_v
37 - HLT_Photon35_TwoProngs35_v
38 - HLT_PFHT500_PFMET100_PFMHT100_IDTight_v
39 - HLT_AK8PFJet400_TrimMass30_v
40 - HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v
41 - HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v

all remaining bits are defined as zero.
"""
