# coding: utf-8

inters_general = {'MET': (('IsoTau180',)),
                  'EG': ( ),
                  'Mu': (('METNoMu120',),
                         ('IsoTau180', 'METNoMu120')),
                  'Tau': ( )
                  }

inters_etau = {'MET': (('Ele32',),
                       ('EleIsoTauCustom',),
                       ('Ele32', 'EleIsoTauCustom'),
                       ('Ele32', 'IsoTau180'),
                       ('EleIsoTauCustom', 'IsoTau180'),
                       ('EleIsoTauCustom', 'IsoTau180'),
                       ('Ele32', 'EleIsoTauCustom', 'IsoTau180')),
               'EG': ( ),
               'Mu': (('Ele32', 'METNoMu120'),
                      ('EleIsoTauCustom', 'METNoMu120'),
                      ('Ele32', 'EleIsoTauCustom', 'METNoMu120'),
                      ('Ele32', 'IsoTau180', 'METNoMu120'),
                      ('EleIsoTauCustom', 'IsoTau180', 'METNoMu120'),
                      ('Ele32', 'EleIsoTauCustom', 'IsoTau180', 'METNoMu120')),
               'Tau': ( )
               }

inters_mutau = {'MET': (('IsoMu24',),
                        ('IsoMuIsoTauCustom',),
                        ('IsoMu24', 'IsoMuIsoTauCustom'),
                        ('IsoMu24', 'IsoTau180'),
                        ('IsoMuIsoTauCustom', 'IsoTau180'),
                        ('IsoMuIsoTauCustom', 'IsoTau180'),
                        ('IsoMu24', 'IsoMuIsoTauCustom', 'IsoTau180')),
                'EG': (('IsoMu24', 'METNoMu120'),
                       ('IsoMuIsoTauCustom', 'METNoMu120'),
                       ('IsoMu24', 'IsoTau180', 'METNoMu120'),
                       ('IsoMu24', 'IsoMuIsoTauCustom', 'METNoMu120'),
                       ('IsoMuIsoTauCustom', 'IsoTau180', 'METNoMu120'),
                       ('IsoMu24', 'IsoMuIsoTauCustom', 'IsoTau180', 'METNoMu120')),
                'Mu': ( ),
                'Tau': ( )}

inters_tautau = {'MET': (('IsoDoubleTauCustom',),
                         ('IsoDoubleTauCustom', 'VBFTauCustom'),
                         ('IsoDoubleTauCustom', 'IsoTau180'),
                         ('IsoDoubleTauCustom', 'IsoTau180', 'VBFTauCustom')),
                 'EG': ( ),
                 'Mu': (('VBFTauCustom',),
                        ('IsoTau180', 'VBFTauCustom')
                        ('METNoMu120', 'VBFTauCustom'),
                        ('IsoTau180', 'METNoMu120', 'VBFTauCustom'),
                        ('IsoDoubleTauCustom', 'METNoMu120'),
                        ('IsoDoubleTauCustom', 'IsoTau180', 'METNoMu120'),
                        ('IsoDoubleTauCustom', 'METNoMu120', 'VBFTauCustom'),
                        ('IsoDoubleTauCustom', 'IsoTau180', 'METNoMu120', 'VBFTauCustom')),
                 'Tau': ( )}
