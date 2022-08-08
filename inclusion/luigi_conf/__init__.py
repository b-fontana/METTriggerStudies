
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
