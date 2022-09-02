# coding: utf-8

inters_general = {'MET' : (),
                  'EG'  : (),
                  'Mu'  : (('IsoTau180',),
                           ('METNoMu120',),
                           ('IsoTau180', 'METNoMu120')),
                  'Tau' : ()
                  }

inters_etau = {'MET' : (),
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
               }

inters_mutau = {'MET' : (),
                'EG'  : (('IsoMu24',),
                         ('IsoMuIsoTauCustom',),
                         ('IsoMu24', 'IsoMuIsoTauCustom'),
                         ('IsoMu24', 'IsoTau180'),
                         ('IsoMuIsoTauCustom', 'IsoTau180'),
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
                }

inters_tautau = {'MET' : (),
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

# inters_with_veto = (('IsoTau180', 'METNoMu120'),
#                     ('METNoMu120', 'VBFTauCustom'),
#                     # mutau channel
#                     ('IsoMu24', 'METNoMu120'),
#                     ('IsoMuIsoTauCustom', 'METNoMu120'),
#                     ('IsoMu24', 'IsoTau180', 'METNoMu120'),
#                     ('IsoMu24', 'IsoMuIsoTauCustom', 'METNoMu120'),
#                     ('IsoMuIsoTauCustom', 'IsoTau180', 'METNoMu120'),
#                     # etau channel
#                     ('Ele32', 'METNoMu120'),
#                     ('Ele32', 'EleIsoTauCustom', 'METNoMu120'),
#                     ('Ele32', 'IsoTau180', 'METNoMu120'),
#                     ('EleIsoTauCustom', 'IsoTau180', 'METNoMu120'),
#                     # tautau channel
#                     ('IsoDoubleTauCustom', 'METNoMu120'),
#                     ('IsoDoubleTauCustom', 'IsoTau180', 'METNoMu120'),
#                     ('IsoDoubleTauCustom', 'METNoMu120', 'VBFTauCustom'),
#                     # large intersections
#                     ('IsoMu24', 'IsoMuIsoTauCustom', 'IsoTau180', 'METNoMu120'),
#                     ('IsoDoubleTauCustom', 'IsoTau180', 'METNoMu120', 'VBFTauCustom'),
#                     ('Ele32', 'EleIsoTauCustom', 'IsoTau180', 'METNoMu120'))
