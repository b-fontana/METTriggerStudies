# coding: utf-8

_all_ = [ 'EventSelection' ]

from collections import defaultdict

import os
import sys
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, parent_dir)

import inclusion
from inclusion import config

class EventSelection:
    def __init__(self, entries, dataset, isdata, debug=False):
        self.entries = entries
        self.bit = self.entries.pass_triggerbit
        self.run = self.entries.RunNumber
        self.isdata = isdata
        self.debug = debug

        self.datasets = ('MET', 'EG', 'Mu', 'Tau')
        self.prefix = 'Data_' if self.isdata else 'TT_'
        self.ds_name = lambda ds : self.prefix + ds
        
        self.this_processed_dataset = self.ds_name(dataset)
        
        for d in self.datasets:
            assert( d in config.data )

        self.ref_trigs = tuple(self.ds_name(x) for x in self.datasets)
        
        if dataset not in config.data and dataset not in config.mc_processes:
            mes = 'Dataset {} is not supported '.format(dataset)
            mes += '(prefix `{}` added).'.format(self.prefix)
            raise ValueError(mes)

        self.noref_str = 'NoReference'

    def any_trigger(self, trigs):
        """
        Checks at least one trigger was fired.
        Considers all framework triggers.
        """
        return self._pass_triggers(trigs)

    def check_bit(self, bitpos):
        bitdigit = 1
        res = bool(self.bit&(bitdigit<<bitpos))
        return res

    def dataset_cuts(self, tcomb, channel):
        """
        Applies selection depending on the reference trigger being considered.
        Reference triggers tend to be associated with a certain dataset.
        For instance, the 'MET' dataset is connected to the MET Trigger.
        """
        reference, lepton_veto = self.find_inters_for_reference(tcomb, channel)

        if reference == self.noref_str:
            return False

        if not any(x in reference for x in self.datasets):
            m = "Only datasets {} are supported. You tried using '{}'.".format(self.datasets, reference)
            raise ValueError(m)
                    
        return self.selection_cuts(lepton_veto=lepton_veto)

    def dataset_triggers(self, trigs, tcomb, channel):
        """
        Checks at least one trigger was fired.
        Considers framework triggers for a specific dataset.
        """
        dataset_ref_trigs = {
            self.ds_name('MET') : ('METNoMu120',),
            self.ds_name('EG')  : ('Ele32',),
            self.ds_name('Mu')  : ('IsoMu24',),
            self.ds_name('Tau') : ('IsoTau180',),
            }
        for k in dataset_ref_trigs:
            if k not in self.ref_trigs:
                mes = 'Reference trigger {} is being set but not defined.'
                raise ValueError(mes.format(k))

        for k in self.ref_trigs:
            if k not in dataset_ref_trigs:
                mes = 'Specify the reference trigger {} for dataset {}!'
                raise ValueError(mes.format(k, self.this_processed_dataset))
        for vals in dataset_ref_trigs.values():
            for v in vals:
                if v not in trigs:
                    mes = 'Reference trigger {} is not part of triggers {}.'
                    raise ValueError(mes.format(v,trigs))

        reference, _ = self.find_inters_for_reference(tcomb, channel)
        if reference == self.noref_str:
            raise OverflowError('Intersection is too long.')

        return ( self._pass_triggers(dataset_ref_trigs[reference]),
                 dataset_ref_trigs[reference] )

    def get_trigger_bit(self, trigger_name):
        """
        Returns the trigger bit corresponding to 'config.trig_map'
        """
        s = 'data' if self.isdata else 'mc'
        res = config.trig_map[trigger_name]
        try:
            res = res[s]
        except KeyError:
            print('You likely forgot to add your custom trigger to `config.trig_custom`.')
            raise
        return res

    def find_inters_for_reference(self, tcomb, channel):
        wrong_comb = 'Combination {} is not supported for channel {}.'

        # Ignore long intersections for simplicity
        # Besides, long intersections tend to have lower statistics
        if len(tcomb) > 3:
            return self.noref_str, None

        # whether to apply 3rd lepton veto
        veto = False
        if tcomb not in (# general triggers
                         ('IsoTau180', 'METNoMu120'),
                         ('METNoMu120', 'VBFTauCustom'),
                         # mutau channel
                         ('IsoMu24', 'METNoMu120'),
                         ('IsoMuIsoTauCustom', 'METNoMu120'),
                         ('IsoMu24', 'IsoTau180', 'METNoMu120'),
                         ('IsoMu24', 'METNoMu120', 'VBFTauCustom'),
                         ('IsoMu24', 'IsoMuIsoTauCustom', 'METNoMu120'),
                         ('IsoMuIsoTauCustom', 'IsoTau180', 'METNoMu120'),
                         # etau channel
                         ('Ele32', 'METNoMu120'),
                         ('Ele32', 'EleIsoTauCustom', 'METNoMu120'),
                         ('Ele32', 'IsoTau180', 'METNoMu120'),
                         ('EleIsoTauCustom', 'IsoTau180', 'METNoMu120'),
                         ('EleIsoTauCustom', 'METNoMu120', 'VBFTauCustom'),
                         ('Ele32', 'METNoMu120', 'VBFTauCustom'),
                         # tautau channel
                         ('IsoDoubleTauCustom', 'METNoMu120'),
                         ('IsoDoubleTauCustom', 'IsoTau180', 'METNoMu120'),
                         ('IsoDoubleTauCustom', 'METNoMu120', 'VBFTauCustom'),
                         # large intersections
                         ('IsoMu24', 'IsoMuIsoTauCustom', 'IsoTau180', 'METNoMu120'),
                         ('IsoDoubleTauCustom', 'IsoTau180', 'METNoMu120', 'VBFTauCustom'),
                         ('Ele32', 'EleIsoTauCustom', 'IsoTau180', 'METNoMu120'),
                         ):
            veto = True
        
        #general triggers
        if tcomb in ( ('IsoTau180',),
                      ('VBFTauCustom',),
                      ('IsoTau180', 'VBFTauCustom'),                     
                     ):
            return self.ds_name('MET'), veto
        elif tcomb in ( ):
            return self.ds_name('EG'), veto
        elif tcomb in ( ):
            return self.ds_name('Mu'), veto
        elif tcomb in ( ('METNoMu120',),
                        ('IsoTau180', 'METNoMu120'),
                        ('METNoMu120', 'VBFTauCustom'),
                        ('IsoTau180', 'METNoMu120', 'VBFTauCustom'),):
            return self.ds_name('Tau'), veto
     
        # channel-specific triggers
        if channel == 'etau':
            if tcomb in ( ('Ele32',),
                          ('EleIsoTauCustom',),
                          ('Ele32', 'EleIsoTauCustom'),
                          ('Ele32', 'VBFTauCustom'),
                          ('Ele32', 'IsoTau180'),
                          ('EleIsoTauCustom', 'VBFTauCustom'),
                          ('EleIsoTauCustom', 'IsoTau180'),
                          ('EleIsoTauCustom', 'IsoTau180'),
                          ('Ele32', 'EleIsoTauCustom', 'VBFTauCustom'),
                          ('Ele32', 'EleIsoTauCustom', 'IsoTau180'),
                          ('Ele32', 'IsoTau180', 'VBFTauCustom'),
                          ('EleIsoTauCustom', 'IsoTau180', 'VBFTauCustom'),
                         ):
                return self.ds_name('MET'), veto
            elif tcomb in ( ) :
                return self.ds_name('EG'), veto
            elif tcomb in ( ):
                return self.ds_name('Mu'), veto
            elif tcomb in (('Ele32', 'METNoMu120'),
                           ('EleIsoTauCustom', 'METNoMu120'),
                           ('Ele32', 'EleIsoTauCustom', 'METNoMu120'),
                           ('Ele32', 'IsoTau180', 'METNoMu120'),
                           ('Ele32', 'METNoMu120', 'VBFTauCustom'),
                           ('EleIsoTauCustom', 'IsoTau180', 'METNoMu120'),
                           ('EleIsoTauCustom', 'METNoMu120', 'VBFTauCustom'),
                           ('Ele32', 'EleIsoTauCustom', 'IsoTau180', 'METNoMu120')):
                return self.ds_name('Tau'), veto
            else:
                raise ValueError(wrong_comb.format(tcomb, channel))
            
        elif channel == 'mutau':
            if tcomb in ( ('IsoMu24',),
                          ('IsoMuIsoTauCustom',),
                          ('IsoMu24', 'IsoMuIsoTauCustom'),
                          ('IsoMu24', 'VBFTauCustom'),
                          ('IsoMu24', 'IsoTau180'),
                          ('IsoMuIsoTauCustom', 'VBFTauCustom'),
                          ('IsoMuIsoTauCustom', 'IsoTau180'),
                          ('IsoMuIsoTauCustom', 'IsoTau180'),
                          ('IsoMu24', 'IsoMuIsoTauCustom', 'VBFTauCustom'),
                          ('IsoMu24', 'IsoMuIsoTauCustom', 'IsoTau180'),
                          ('IsoMu24', 'IsoTau180', 'VBFTauCustom'),
                          ('IsoMuIsoTauCustom', 'IsoTau180', 'VBFTauCustom'),
                         ):
                return self.ds_name('MET'), veto
            elif tcomb in ( ) :
                return self.ds_name('EG'), veto
            elif tcomb in (('IsoMu24', 'METNoMu120'),
                           ('IsoMuIsoTauCustom', 'METNoMu120'),
                           ('IsoMu24', 'IsoTau180', 'METNoMu120'),
                           ('IsoMu24', 'METNoMu120', 'VBFTauCustom'),
                           ('IsoMu24', 'IsoMuIsoTauCustom', 'METNoMu120'),
                           ('IsoMuIsoTauCustom', 'IsoTau180', 'METNoMu120'),
                           ('IsoMuIsoTauCustom', 'METNoMu120', 'VBFTauCustom'),
                           ('IsoMu24', 'IsoMuIsoTauCustom', 'IsoTau180', 'METNoMu120'),):
                return self.ds_name('Mu'), veto
            elif tcomb in ( ):
                return self.ds_name('Tau'), veto
            else:
                raise ValueError(wrong_comb.format(tcomb, channel))
            
        elif channel == 'tautau':
            if tcomb in (('IsoDoubleTauCustom',),
                         ('IsoDoubleTauCustom', 'VBFTauCustom'),
                         ('IsoDoubleTauCustom', 'IsoTau180'),
                         ('IsoDoubleTauCustom', 'IsoTau180', 'VBFTauCustom'),
                         ):
                return self.ds_name('MET'), veto
            elif tcomb in ( ) :
                return self.ds_name('EG'), veto
            elif tcomb in ( ):
                return self.ds_name('Mu'), veto
            elif tcomb in (('IsoDoubleTauCustom', 'METNoMu120'),
                           ('IsoDoubleTauCustom', 'IsoTau180', 'METNoMu120'),
                           ('IsoDoubleTauCustom', 'METNoMu120', 'VBFTauCustom'),
                           ('IsoDoubleTauCustom', 'IsoTau180', 'METNoMu120', 'VBFTauCustom'),):
                return self.ds_name('Tau'), veto
            else:
                raise ValueError(wrong_comb.format(tcomb, channel))
            
        else:
            raise ValueError('Channel {} is not supported.'.format(channel))

    def check_inters_with_dataset(self, tcomb, channel):
        """
        All input files on which the selection is applied correspond
        to a different dataset.  This function makes sure there is a
        match between the trigger intersection and the dataset being
        used. It is used to skip the processing of trigger
        intersections which are calculated with other datasets.
        
        Example:
        The 'IsoMu24' trigger intersection will use the MET
        dataset in some analysis. This functions skips the processing
        whenever we are running the selection over other datasets,
        such as EGamma or SingleMuon.
        """
        reference, _ = self.find_inters_for_reference(tcomb, channel)
        if reference == self.noref_str:
            return False

        return True if self.this_processed_dataset == reference else False

    def _pass_triggers(self, trigs):
        """
        Internal only function.
        Checks at least one trigger was fired.
        """
        flag = False
        for trig in trigs:
            if trig in config.trig_custom:
                flag = set_custom_trigger_bit(trig, self.bit, self.run)
            else:
                flag = self.check_bit(self.get_trigger_bit(trig))
            if flag:
                return True
        return False    

    def selection_cuts(self, iso_cuts=dict(), lepton_veto=True,
                       bjets_cut=True, invert_mass_cut=True):
        """
        Applies selection cut to one event.
        Returns `True` only if all selection cuts pass.
        """
        mhh = self.entries.HHKin_mass
        if mhh < 1:
            return False

        pairtype    = self.entries.pairType
        dau1_eleiso = self.entries.dau1_eleMVAiso
        dau1_muiso  = self.entries.dau1_iso
        dau1_tauiso = self.entries.dau1_deepTauVsJet
        dau2_tauiso = self.entries.dau2_deepTauVsJet

        # third lepton veto
        nleps = self.entries.nleps
        if nleps > 0 and lepton_veto:
            return False

        # require at least two b jet candidates
        nbjetscand = self.entries.nbjetscand
        if nbjetscand <= 1 and bjet_cuts:
            return False

        # Loose / Medium / Tight
        iso_allowed = { 'dau1_ele': 1., 'dau1_mu': 0.15,
                        'dau1_tau': 5., 'dau2_tau': 5. }
        if any(x not in iso_allowed for x in iso_cuts.keys()):
            mes = 'At least one of the keys is not allowed. '
            mes += 'Keys introduced: {}.'.format(iso_cuts.keys())
            raise ValueError(mes)

        # setting to the defaults in case the user did not specify the values
        for k, v in iso_allowed.items():
            if k not in iso_cuts: iso_cuts[k] = v
        
        bool0 = pairtype==0 and (dau1_muiso >= iso_cuts['dau1_mu'] or
                                 dau2_tauiso < iso_cuts['dau2_tau'])
        bool1 = pairtype==1 and (dau1_eleiso != iso_cuts['dau1_ele'] or
                                 dau2_tauiso < iso_cuts['dau2_tau'])
        bool2 = pairtype==2 and (dau1_tauiso < iso_cuts['dau1_tau'] or
                                 dau2_tauiso < iso_cuts['dau2_tau'])
        if bool0 or bool1 or bool2:
            return False

        #((tauH_SVFIT_mass-116.)*(tauH_SVFIT_mass-116.))/(35.*35.) + ((bH_mass_raw-111.)*(bH_mass_raw-111.))/(45.*45.) <  1.0
        svfit_mass = self.entries.tauH_SVFIT_mass
        bh_mass    = self.entries.bH_mass_raw

        mcut = ( (svfit_mass-129.)*(svfit_mass-129.) / (53.*53.) +
                 (bh_mass-169.)*(bh_mass-169.) / (145.*145.) ) < 1.0
        if mcut and invert_mass_cut: # inverted elliptical mass cut (-> ttCR)
            return False

        return True

    def set_custom_trigger_bit(self, trigger):
        """
        The VBF trigger was updated during data taking, adding HPS
        https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauTrigger
        """
        if trigger not in config.trig_custom:
            import inspect
            currentFunction = inspect.getframeinfo(frame).function
            raise ValueError('[{}] option not supported.'.format(currentFunction))

        if self.run < 317509 and self.isdata:
            if trigger == 'VBFTauCustom':
                bits = self.check_bit(config.trig_map[trigger]['VBFTau']['data'])
            elif trigger == 'IsoDoubleTauCustom':
                bits = ( self.check_bit(config.trig_map[trigger]['IsoDoubleTau']['data'][0]) or
                         self.check_bit(config.trig_map[trigger]['IsoDoubleTau']['data'][1]) or
                         self.check_bit(config.trig_map[trigger]['IsoDoubleTau']['data'][2]) )
            elif trigger == 'IsoMuIsoTauCustom':
                bits = self.check_bit(config.trig_map[trigger]['IsoMuIsoTau']['data'])
            elif trigger == 'EleIsoTauCustom':
                bits = self.check_bit(config.trig_map[trigger]['EleIsoTau']['data'])

        else:
            s = 'data' if self.isdata else 'mc'
            if trigger == 'VBFTauCustom':
                bits = self.check_bit(config.trig_map[trigger]['VBFTauHPS'][s])
            elif trigger == 'IsoDoubleTauCustom':
                bits = self.check_bit(config.trig_map[trigger]['IsoDoubleTauHPS'][s])
            elif trigger == 'IsoMuIsoTauCustom':
                bits = self.check_bit(config.trig_map[trigger]['IsoMuIsoTauHPS'][s])
            elif trigger == 'EleIsoTauCustom':
                bits = self.check_bit(config.trig_map[trigger]['EleIsoTauHPS'][s])

        return bits

    def trigger_bits(self, trig):
        if trig in config.trig_custom:
            return self.set_custom_trigger_bit(trig)
        else:
            return self.check_bit(self.get_trigger_bit(trig))

    def var_cuts(self, trig, variables, nocut_dummy_str):
        """
        Handles cuts on trigger variables that enter the histograms. 
        Variables being displayed are not cut (i.e., they pass the cut).
        Checks all combinations of cuts specified in '_cuts':
            example: _cuts = {'A': ('>', [10,20]), 'B': ('<', [50,40]))}
            `passes_cuts` will check 4 combinations and return a dict of length 4
            (unless some cuts are ignored according to '_cuts_ignored') 
        Works for both 1D and 2D efficiencies.
        """
        if self.debug:
            print('Trigger={}; Variables={}'.format(trig, variables))

        flagnameJoin = lambda var,sign,val: ('_'.join([str(x) for x in [var,sign,val]])).replace('.','p')
    
        dflags = defaultdict(lambda: [])
    
        try:
            trig_cuts = config.cuts[trig]
        except KeyError: # the trigger has no cut associated
            if self.debug:
                print('KeyError')            
            return {nocut_dummy_str: True}

        for avar,acut in trig_cuts.items():
            # ignore cuts according to the user's definition in '_cuts_ignored'
            # example: do not cut on 'met_et' when displaying 'metnomu_et'
            ignore = functools.reduce( lambda x, y: x or y,
                                       [ avar in config.cuts_ignored[k] for k in variables
                                         if k in config.cuts_ignored ],
                                      False
                                      )

            # additionally, by default do not cut on the variable(s) being plotted
            if avar not in variables and not ignore:
                value = self.entries[avar]

                for c in acut[1]:
                    flagname = flagnameJoin(avar, acut[0], c)

                    if self.debug:
                        print('Cut: {} {} {}'.format(avar, acut[0], c))

                    if acut[0]=='>':
                        dflags[avar].append( (flagname, value > c) )
                    elif acut[0]=='<':
                        dflags[avar].append( (flagname, value < c) )
                    else:
                        mes = 'The operator for the cut is currently '
                        mess += 'not supported: Use `>` or `<`.' 
                        raise ValueError(mes)

        def apply_cuts_combinations(dflags):
            tmp = {}
            allNames = sorted(dflags)
            combinations = it.product(*(dflags[name] for name in allNames))
            combinations = list(combinations)

            if self.debug:
                print('--------- [applyCutsCombinations] ------------------')
                print('Variables being cut: {}'.format(allNames))
                print('All cut combinations: {}'.format(combinations))
                print('----------------------------------------------------')
            
            for comb in combinations:
                joinFlag = functools.reduce( lambda x,y: x and y, [k[1] for k in comb] )
                tmp[ '_AND_'.join([k[0] for k in comb]) ] = joinFlag

            return tmp

        if dflags:
            res = apply_cuts_combinations(dflags)
        else:
            res = {nocut_dummy_str: True}

        return res
