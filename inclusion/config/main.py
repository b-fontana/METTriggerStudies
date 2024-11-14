import os

email = 'bruno.alves@cern.ch'
queue = 'short'
machine = 'slurm' #'lxplus' or 'llrt3condor'

storage = {
           '2022':    os.path.join('/t3home/', os.environ['USER'], 'TriggerScaleFactors'),
           '2018':    os.path.join('/data_CMS/cms/', os.environ['USER'], 'TriggerScaleFactors'),
           '2017':    os.path.join('/data_CMS/cms/', os.environ['USER'], 'TriggerScaleFactors'),
           '2016':    os.path.join('/data_CMS/cms/', os.environ['USER'], 'TriggerScaleFactors'),
           '2016APV': os.path.join('/data_CMS/cms/', os.environ['USER'], 'TriggerScaleFactors')}

folders = {'base'    : 'METTriggerStudies',
           'main'    : 'inclusion',
           'scripts' : 'scripts',
           'jobs'    : 'condor',
           'subm'    : 'submission',
           'outs'    : 'outputs'}

base_folder = {
              'slurm':
               os.path.join(os.environ['HOME'], 'HHbbtautau/inclusionRun3'),
              'llrt3condor':
               os.path.join(os.environ['HOME'], 'CMSSW_14_1_0_pre0', 'src', folders['base']),
               'llrt3condor7':
               os.path.join(os.environ['HOME'], 'CMSSW_14_1_0_pre0', 'src', folders['base']),
               'lxplus':
               os.path.join('/afs/cern.ch/work/', os.environ['USER'][0], os.environ['USER'],
                            'CMSSW_14_1_0_pre0', 'src', folders['base']),}
local_folder = os.path.join(base_folder[machine], folders['main'])

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

channels = ('etau', 'mutau', 'tautau', 'mumu')
sel = {'all'    : {'pairType': ('<',  3),},
       'mutau'  : {'pairType': ('==', 0),},
       'etau'   : {'pairType': ('==', 1),},
       'tautau' : {'pairType': ('==', 2),},
       'mumu'   : {'pairType': ('==', 3),}, # passMu missing for the mumu channel
       'ee'     : {'pairType': ('==', 4),},
       'emu'    : {'pairType': ('==', 5),} }

# variables considered for calculating and plotting efficiencies
var_eff = (
    # 'met_et', 'mht_et', 'metnomu_et', 'mhtnomu_et',
           'dau1_pt', 'dau2_pt', 'dau1_eta', 
           'dau1_tauIdVSjet', 'dau2_tauIdVSjet',
           'bjet1_btagDeepFlavB', 'bjet2_btagDeepFlavB'
        #    'dau2_eta'
           )

# variables considered for plotting MC/data comparison distributions
var_dist = ('dau1_pt', 
# 'HH_mass'
)
# joining the two lists above
var_join = set(var_eff + var_dist)

var_unionweights = ('dau1_pt', 'dau2_pt', 'dau1_eta', 
# 'dau2_eta'
)

trig_linear = lambda x : {'mc': x, 'data': x}

# The bit numbering scheme below matches the 'triggerbit' leaf
# It does NOT match the 'pass_triggerbit' leaf, which is a skimmed version of the above that might change more often
# One way to ensure the scheme is still correct is by running the script as shown here:
# https://github.com/bfonta/useful_scripts/commit/f5e4a0096bc74c89176579a336b0f52b74cb3ed2
filterbit_map = {
    '2022': {
        'IsoMu24': {'mc': {13: (1, 3)},  'data': {13: (1, 3)}},
        'IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1': {'mc': {13: (2, ), 15: (3, 5, 9)}, 'data': {13: (2, ), 15: (3, 5, 9)}},
        "LooseDeepTauPFTauHPS180_L2NN_eta2p1": {'mc': {}, 'data': {}}, 
        "PFMETNoMu120_PFMHTNoMu120_IDTight": {'mc': {}, 'data': {}}, 
    }
}

trig_map = {
       '2022': {
            'IsoMu24': {'mc': 0,  'data': 0},
            'IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1': {'mc': 1, 'data': 1},
            "LooseDeepTauPFTauHPS180_L2NN_eta2p1": {'mc': 11, 'data': 11}, 
            "PFMETNoMu120_PFMHTNoMu120_IDTight": {'mc': 12, 'data': 12}, 
       },
       '2018':
            {'IsoMu24':    {'mc': 0,  'data': 0},
             'Ele32':      {'mc': 2,  'data': 2},
             'METNoMu120': {'mc': 40, 'data': 40},
             'IsoTau180':  {'mc': (4,7), 'data': (4,7)},
             'IsoDoubleTauCustom': {'IsoDoubleTau':    {'mc': 12, 'data': (13,14,15)},
                                    'IsoDoubleTauHPS': {'mc': 12, 'data': 12}},
             'IsoMuIsoTauCustom':  {'IsoMuIsoTau':     {'mc': 8,  'data': 9},
                                    'IsoMuIsoTauHPS':  {'mc': 8,  'data': 8} },
             'EleIsoTauCustom':    {'EleIsoTau':       {'mc': 10, 'data': 11},
                                    'EleIsoTauHPS':    {'mc': 10, 'data': 10}},
             'VBFTauCustom':       {'VBFTau':          {'mc': 16, 'data': 12},
                                    'VBFTauHPS':       {'mc': 16, 'data': 13}}}
            ,

            '2017':
            {'IsoMu27':      {'mc': 0,  'data': 0},
             'Ele32':        {'mc': 2,  'data': 2},
             'METNoMu120':   {'mc': 40, 'data': 40},
             'IsoTau180':    {'mc': (31,32), 'data': (31,32)},
             'IsoDoubleTau': {'mc': (25,27,29), 'data': (25,27,29)},
             'IsoMuIsoTau':  {'mc': 5, 'data': 5},
             'EleIsoTau':    {'mc': 17, 'data': 17}}
            ,
            '2016':
            {'IsoMu24':      {'mc': 5,  'data': 5},
             'Ele25':        {'mc': 13, 'data': 13},
             'METNoMu90':    {'mc': 48, 'data': 48},
             'IsoTau130':    {'mc': (46,47),  'data': (46,47)},
             'IsoDoubleTau': {'mc': (27,29), 'data': (27,29)},
             'IsoMuIsoTau':  {'mc': (19,20), 'data': (19,20)}}
            }
trig_map['2016APV'] = trig_map['2016']

lep_triggers = {
                '2022': {
                    'IsoMu24', 'IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1'
                },
                '2018':
                {'Ele32', 'EleIsoTauCustom', 'IsoMu24', 'IsoMuIsoTauCustom', 'IsoDoubleTauCustom'},
                '2017':
                {'Ele32', 'EleIsoTau', 'IsoMu27', 'IsoMuIsoTau', 'IsoDoubleTau'},
                '2016':
                {'Ele25', 'IsoMu24', 'IsoMuIsoTau', 'IsoDoubleTau'}
                }
for year in ('2016', '2017', '2018', '2022'):
    assert all(x in trig_map[year].keys() for x in lep_triggers[year])

cuts_ignored = {#'HT20':       (),
                # 'met_et':     ('metnomu_et',),
                # 'mht_et':     ('mhtnomu_et',),
                # 'metnomu_et': ('met_et',),
                # 'mhtnomu_et': ('mht_et',),
                'dau1_pt':    (),
                'dau2_pt':    ()}

### Correlation Matrix
corr = {'etau': {},
        'mutau': {},
        'tautau': {},
        'mumu': {} }
    
### Data and MC samples
inputs = {
        '2022': ('/t3home/fbilandz/HHbbtautau/Run3/MC/PreprocessOutputs'),
        '2018':
          ('/data_CMS/cms/portales/HHresonant_SKIMS/SKIMS_UL18_OpenCADI_Data/',
           '/data_CMS/cms/alves/HHresonant_SKIMS/SKIMS_UL18_OpenCADI_MC/'),
          '2017':
          ('/eos/home-t/tokramer/hhbbtautau/skims/SKIMS_UL17/',),
          '2016APV':
          ('/eos/home-s/spalluot/HHbbtautau/SKIM/SKIMS_UL2016APV_14Feb2024/',),
          '2016':
          ('/eos/home-d/dzuolo/SKIMS_UL2016_12Feb2024_newMET_newBoosted/',)
          }
corrupted_files = () #'/data_CMS/cms/alves/HHresonant_SKIMS/SKIMS_UL18_EOSv5HighPrio_Background/TTTo2L2Nu/output_9.root'

# names of the subfolders under 'inputs' above:
# dictionary that maps specific general triggers to datasets 
data = {
    "2022":
        {'MET'  : ('MET',),
         'EG'   : ('EGamma',),
         'Mu'   : (
            #  'SingleMuon',
               'Muon',),
         'Tau'  : ('Tau',)
         },
    "2018":
        {'MET'  : ('MET',),
         'EG'   : ('EGamma',),
         'Mu'   : ('SingleMuon',),
         'Tau'  : ('Tau',)
         },
        "2017":
        {'MET'  : ('MET',),
         'EG'   : ('EGamma',),
         'Mu'   : ('SingleMuon',),
         'Tau'  : ('Tau',)
         },
        "2016APV":
        {'MET'  : ('SKIM_MET',),
         'EG'   : ('SKIM_Electron',),
         'Mu'   : ('SKIM_SingleMuon',),
         'Tau'  : ('SKIM_Tau',)
         },
        "2016":
        {'MET'  : ('SKIM_MET',),
         'EG'   : ('SKIM_Ele',),
         'Mu'   : ('SKIM_Mu',),
         'Tau'  : ('SKIM_Tau',)
         },
}

mc_processes = {
    "2022": {
        'TT': ('TTToHadronic', 'TTToFullyLeptonic', 
        'TTToSemiLeptonic',),
        'DY': ('DYJetsToLL_M-50')
    },
    "2018":
                {'ggfRadions': (), 'ggfBulkGraviton': (),
                 'TT': ('TTToHadronic', 'TTTo2L2Nu', 'TTToSemiLeptonic',),
                 'DY': ('DYJetsToLL_M-50_TuneCP5_13TeV-amc',
                        'DYJetsToLL_0J', 'DYJetsToLL_1J',   'DYJetsToLL_2J',
                        'DYJetsToLL_LHEFilterPtZ-0To50',    'DYJetsToLL_LHEFilterPtZ-50To100',
                        'DYJetsToLL_LHEFilterPtZ-100To250', 'DYJetsToLL_LHEFilterPtZ-250To400',
                        'DYJetsToLL_LHEFilterPtZ-400To650', 'DYJetsToLL_LHEFilterPtZ-650ToInf'),
                 'WJets': ('WJetsToLNu_TuneCP5_13TeV-madgraph',
                           'WJetsToLNu_HT-70To100',    'WJetsToLNu_HT-100To200', 'WJetsToLNu_HT-200To400',
                           'WJetsToLNu_HT-400To600',   'WJetsToLNu_HT-600To800', 'WJetsToLNu_HT-800To1200',
                           'WJetsToLNu_HT-1200To2500', 'WJetsToLNu_HT-2500ToInf'),
                 },
                
"2017":
                {'ggfRadions': (), 'ggfBulkGraviton': (),
                 'TT': ('TTToHadronic', 'TTTo2L2Nu', 'TTToSemiLeptonic',),
                 'DY': ('DYJetsToLL_M-50_TuneCP5_13TeV-amc',
                        'DYJetsToLL_0J', 'DYJetsToLL_1J',   'DYJetsToLL_2J',
                        'DYJetsToLL_LHEFilterPtZ-0To50',    'DYJetsToLL_LHEFilterPtZ-50To100',
                        'DYJetsToLL_LHEFilterPtZ-100To250', 'DYJetsToLL_LHEFilterPtZ-250To400',
                        'DYJetsToLL_LHEFilterPtZ-400To650', 'DYJetsToLL_LHEFilterPtZ-650ToInf'),
                 'WJets': ('WJetsToLNu_TuneCP5_13TeV-madgraph',
                           'WJetsToLNu_HT-70To100',    'WJetsToLNu_HT-100To200', 'WJetsToLNu_HT-200To400',
                           'WJetsToLNu_HT-400To600',   'WJetsToLNu_HT-600To800', 'WJetsToLNu_HT-800To1200',
                           'WJetsToLNu_HT-1200To2500', 'WJetsToLNu_HT-2500ToInf'),
                 },

"2016APV":
                {'ggfRadions': (), 'ggfBulkGraviton': (),
                 'TT': ('SKIM_TTToHadronic', 'SKIM_TTTo2L2Nu', 'SKIM_TTToSemiLeptonic',),
                 'DY': ('SKIM_DYJetsToLL_M-50',
                        'SKIM_DYJetsToLL_0J',  'SKIM_DYJetsToLL_1J', 'SKIM_DYJetsToLL_2J',
                        'SKIM_DYJetsToLL_LHEFilterPtZ-0To50',    'SKIM_DYJetsToLL_LHEFilterPtZ-50To100',
                        'SKIM_DYJetsToLL_LHEFilterPtZ-100To250', 'SKIM_DYJetsToLL_LHEFilterPtZ-250To400',
                        'SKIM_DYJetsToLL_LHEFilterPtZ-400To650', 'SKIM_DYJetsToLL_LHEFilterPtZ-650ToInf'),
                 'WJets': ('SKIM_WJetsToLNu',
                           'SKIM_WJetsToLNu_HT-70To100',    'SKIM_WJetsToLNu_HT-100To200',
                           'SKIM_WJetsToLNu_HT-200To400',   'SKIM_WJetsToLNu_HT-400To600',
                           'SKIM_WJetsToLNu_HT-600To800',   'SKIM_WJetsToLNu_HT-800To1200',
                           'SKIM_WJetsToLNu_HT-1200To2500', 'SKIM_WJetsToLNu_HT-2500ToInf'),
                 },

"2016":
                {'ggfRadions': (), 'ggfBulkGraviton': (),
                 'TT': ('SKIM_TT_fullyHad', 'SKIM_TT_fullyLep', 'SKIM_TT_semiLep',),
                 'DY': ('SKIM_DY_amc_incl',
                        'SKIM_DY_amc_0j',  'SKIM_DY_amc_1j', 'SKIM_DY_amc_2j',
                        'SKIM_DY_amc_PtZ_0To50',    'SKIM_DY_amc_PtZ_50To100',
                        'SKIM_DY_amc_PtZ_100To250', 'SKIM_DY_amc_PtZ_250to400',
                        'SKIM_DY_amc_PtZ_400to650', 'SKIM_DY_amc_PtZ_650toInf'),
                 'WJets': ('SKIM_WJets_HT0To70',
                           'SKIM_WJets_HT70To100',    'SKIM_WJets_HT100To200',
                           'SKIM_WJets_HT200To400',   'SKIM_WJets_HT400To600',
                           'SKIM_WJets_HT600To800',   'SKIM_WJets_HT800To1200',
                           'SKIM_WJets_HT1200To2500', 'SKIM_WJets_HT2500ToInf'),
                 },

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
