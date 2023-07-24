# coding: utf-8

_all_ = [ 'EventSelection' ]

from collections import defaultdict

import os
import sys
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, parent_dir)

import inclusion
from inclusion import config
from inclusion.config import main

class EventSelection:
    def __init__(self, entries, isdata, configuration=None, debug=False):
        self.entries = entries
        self.bit = self.entries['triggerbit']
        self.run = self.entries['RunNumber']
        self.isdata = isdata
        self.debug = debug
        self.prefix = 'Data_' if self.isdata else 'MC_'
        
        # dependency injection
        self.cfg = configuration

        self.datasets, self.dataset_ref_trigs = self._deduce_datasets(self.cfg.inters_general,
                                                                      self.cfg.inters)
        
        self.categories = ('baseline', 's1b1jresolvedMcut', 's2b0jresolvedMcut', 'sboostedLLMcut')
        
        for d in self.datasets:
            assert d in main.data

        
    def any_trigger(self, trigs):
        """
        Checks at least one trigger was fired.
        Considers all framework triggers.
        """
        return self.pass_triggers(trigs)

    def check_bit(self, bitpos):
        bitdigit = 1
        res = bool(self.bit&(bitdigit<<bitpos))
        return res
    
    def dataset_cuts(self, tcomb, channel):
        """
        Applies selection depending on the reference trigger being considered.
        Datasets are defined according to the applied trigger(s).
        For instance, the 'MET' dataset is by construction a group of events to which
        the MET Trigger was applied.
        Currently all datasets have the same selection modulos the lepton veto.
        """
        reference = self.find_inters_for_reference(tcomb, channel)
        if reference is None:
            return False

        if not any(x in reference for x in self.datasets):
            m = "Only datasets {} are supported. You tried using '{}'.".format(self.datasets, reference)
            raise ValueError(m)

        lepton_veto = self.should_apply_lepton_veto(tcomb)
        
        return self.selection_cuts(lepton_veto=lepton_veto,
                                   bjets_cut=self.cfg.bjets_cut)

    def dataset_name(self, dataset):
        if dataset not in main.data and dataset not in main.mc_processes:
            mes = 'Dataset {} is not supported '.format(dataset)
            mes += '(prefix `{}` added).'.format(self.prefix)
            raise ValueError(mes)

        return self.prefix + dataset
    
    def dataset_triggers(self, tcomb, channel, trigs, dataset):
        """
        Checks at least one trigger was fired.
        Considers framework triggers for a specific dataset.
        """
        this_processed_dataset = self.dataset_name(dataset)
        for vals in self.dataset_ref_trigs.values():
            for v in vals:
                if v not in trigs:
                    mes = 'Reference trigger {} is not part of triggers {}.'
                    raise ValueError(mes.format(v,trigs))

        reference = self.find_inters_for_reference(tcomb, channel)
        if reference is None:
            raise OverflowError('Intersection is too long.')

        in_lep = all(x in main.lep_triggers for x in self.dataset_ref_trigs[reference])
        lept = self.entries['isLeptrigger'] if in_lep else True
        pass_trg = lept and self.pass_triggers(self.dataset_ref_trigs[reference])       
        return pass_trg, self.dataset_ref_trigs[reference]

    def _deduce_datasets(self, int_gen, int_chn):
        """
        Deduce the required datasets to be looped over based on the triggers
        defined in the configuration. Only datasets with at least
        one trigger are used.
        Reference triggers for each dataset are also defined.
        Example: If no trigger is defined on the MET dataset in the
        cnofiguration file, the MET dataset is ignored.
        """
        chns = int_chn.keys()
        for chn in chns:
            assert int_gen.keys() == int_chn[chn].keys()

        _ref_trigs = {'MET' : ('METNoMu120',),
                      'EG'  : ('Ele32',),
                      'Mu'  : ('IsoMu24',),
                      'Tau' : ('IsoTau180',)}

        ds = []
        ds_ref_trigs = {}
        for key,val in int_gen.items():
            if len(val) + sum([len(int_chn[chn][key]) for chn in chns]) > 0:
                ds.append(key)
                ds_ref_trigs.update({self.dataset_name(key): _ref_trigs[key]})
        return tuple(ds), ds_ref_trigs
        
    def get_trigger_bit(self, trigger_name):
        """
        Returns the trigger bit corresponding to 'main.trig_map'
        """
        s = 'data' if self.isdata else 'mc'
        res = main.trig_map[trigger_name]
        try:
            res = res[s]
        except KeyError:
            print('You likely forgot to add your custom trigger to `self.cfg.trig_custom`.')
            raise
        return res

    def find_inters_for_reference(self, tcomb, channel):
        wrong_comb = 'Combination {} is not supported for channel {}.'

        # Ignore long intersections for simplicity
        # Besides, long intersections tend to have lower statistics
        if len(tcomb) > 4:
            return None
        
        # general triggers
        for k in main.data:
            if tcomb in self.cfg.inters_general[k]:
                return self.dataset_name(k)
                 
        # channel-specific triggers
        for k in main.data:
            if tcomb in self.cfg.inters[channel][k]:
                return self.dataset_name(k)
        raise ValueError(wrong_comb.format(tcomb, channel))

    def check_inters_with_dataset(self, tcomb, channel, dataset):
        """
        All input files on which the selection is applied correspond
        to a different dataset.  This function makes sure there is a
        match between the trigger intersection and the dataset being
        used. It is used to skip the processing of trigger
        intersections which are calculated with other datasets.
        No MC datasets are skipped.
        
        Example:
        The 'IsoMu24' trigger intersection will use the MET
        dataset in some analysis. This functions skips the processing
        whenever we are running the selection over other datasets,
        such as EGamma or SingleMuon.
        """
        if not self.isdata:
            return True

        this_processed_dataset = self.dataset_name(dataset)
        reference = self.find_inters_for_reference(tcomb, channel)
        if reference is None:
            return False

        return True if this_processed_dataset == reference else False

    def pass_triggers(self, trigs):
        """
        Checks at least one trigger was fired.
        """
        flag = False
        for trig in trigs:
            if trig in self.cfg.trig_custom:
                flag = self.set_custom_trigger_bit(trig)
            else:
                flag = self.check_bit(self.get_trigger_bit(trig))
            if flag:
                return True
        return False    

    def sel_category(self, category):
        assert category in self.categories
        btagLL = self.entries['bjet1_bID_deepFlavor'] > 0.0490 and self.entries['bjet2_bID_deepFlavor'] > 0.0490
        btagM  = ((self.entries['bjet1_bID_deepFlavor'] > 0.2783 and self.entries['bjet2_bID_deepFlavor'] < 0.2783) or
                  (self.entries['bjet1_bID_deepFlavor'] < 0.2783 and self.entries['bjet2_bID_deepFlavor'] > 0.2783))
        btagMM = self.entries['bjet1_bID_deepFlavor'] > 0.2783 and self.entries['bjet2_bID_deepFlavor'] > 0.2783
        
        common = not (self.entries['isVBF'] == 1 and self.entries['VBFjj_mass'] > 500 and self.entries['VBFjj_deltaEta'] > 3 and
                      (self.entries['bjet1_bID_deepFlavor'] > 0.2783 or self.entries['bjet2_bID_deepFlavor'] > 0.2783))

        if category == 'baseline':
            specific = True
        elif category == 's1b1jresolvedMcut':
            specific = self.isBoosted != 1 and btagM
        elif category == 's2b0jresolvedMcut':
            specific = self.isBoosted != 1 and btagMM
        elif category == 'sboostedLLMcut':
            specific = self.isBoosted == 1 and btagLL

        return common and specific

    def selection_cuts(self, iso_cuts=dict(), lepton_veto=True, bjets_cut=True,
                       invert_mass_cut=True, standard_mass_cut=False):
        """
        Applies selection cut to one event.
        Returns `True` only if all selection cuts pass.
        """
        assert invert_mass_cut is not standard_mass_cut
        
        # When one only has 0 or 1 bjet th HH mass is not well defined,
        # and a value of -1 is assigned. One thus has to remove the cut below
        # when considering events with less than 2 b-jets.
        mhh = self.entries['HHKin_mass']
        if mhh < 1 and bjets_cut:
            return False

        pairtype    = self.entries['pairType']
        dau1_eleiso = self.entries['dau1_eleMVAiso']
        dau1_muiso  = self.entries['dau1_iso']
        dau1_tauiso = self.entries['dau1_deepTauVsJet']
        dau2_tauiso = self.entries['dau2_deepTauVsJet']

        # third lepton veto
        nleps = self.entries['nleps']
        if nleps > 0 and lepton_veto:
            return False

        # require at least two b jet candidates
        nbjetscand = self.entries['nbjetscand']
        if nbjetscand <= 1 and bjets_cut:
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
        svfit_mass = self.entries['tauH_SVFIT_mass']
        bh_mass    = self.entries['bH_mass_raw']

        mpoint = ((svfit_mass-129.)*(svfit_mass-129.) / (53.*53.) +
                  (bh_mass-169.)*(bh_mass-169.) / (145.*145.))
        if mpoint < 1.0 and invert_mass_cut: # inverted elliptical mass cut
            return False
        if mpoint > 1.0 and standard_mass_cut: # standard elliptical mass cut
            return False

        return True

    def set_custom_trigger_bit(self, trigger):
        """
        The VBF trigger was updated during data taking, adding HPS
        https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauTrigger
        """
        if trigger not in self.cfg.trig_custom:
            import inspect
            currentFunction = inspect.getframeinfo(frame).function
            raise ValueError('[{}] option not supported.'.format(currentFunction))

        if self.run < 317509 and self.isdata:
            if trigger == 'VBFTauCustom':
                bits = self.check_bit(main.trig_map[trigger]['VBFTau']['data'])
            elif trigger == 'IsoDoubleTauCustom':
                bits = ( self.check_bit(main.trig_map[trigger]['IsoDoubleTau']['data'][0]) or
                         self.check_bit(main.trig_map[trigger]['IsoDoubleTau']['data'][1]) or
                         self.check_bit(main.trig_map[trigger]['IsoDoubleTau']['data'][2]) )
            elif trigger == 'IsoMuIsoTauCustom':
                bits = self.check_bit(main.trig_map[trigger]['IsoMuIsoTau']['data'])
            elif trigger == 'EleIsoTauCustom':
                bits = self.check_bit(main.trig_map[trigger]['EleIsoTau']['data'])

        else:
            s = 'data' if self.isdata else 'mc'
            if trigger == 'VBFTauCustom':
                bits = self.check_bit(main.trig_map[trigger]['VBFTauHPS'][s])
            elif trigger == 'IsoDoubleTauCustom':
                bits = self.check_bit(main.trig_map[trigger]['IsoDoubleTauHPS'][s])
            elif trigger == 'IsoMuIsoTauCustom':
                bits = self.check_bit(main.trig_map[trigger]['IsoMuIsoTauHPS'][s])
            elif trigger == 'EleIsoTauCustom':
                bits = self.check_bit(main.trig_map[trigger]['EleIsoTauHPS'][s])

        return bits

    def should_apply_lepton_veto(self, tcomb):
        """Whether to apply 3rd lepton veto. The veto is always applied to MC."""
        if not self.isdata:
            return True
        
        if (tcomb in self.cfg.inters_general['MET'] or
            tcomb in self.cfg.inters['etau']['MET'] or
            tcomb in self.cfg.inters['mutau']['MET'] or
            tcomb in self.cfg.inters['tautau']['MET']):
            return True
        return False
    
    def trigger_bits(self, trig):
        if trig in self.cfg.trig_custom:
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
            trig_cuts = self.cfg.cuts[trig]
        except KeyError: # the trigger has no cut associated
            if self.debug:
                print('KeyError')            
            return {nocut_dummy_str: True}

        for avar,acut in trig_cuts.items():
            # ignore cuts according to the user's definition in '_cuts_ignored'
            # example: do not cut on 'met_et' when displaying 'metnomu_et'
            ignore = functools.reduce( lambda x, y: x or y,
                                       [ avar in main.cuts_ignored[k] for k in variables
                                         if k in main.cuts_ignored ],
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
                tmp[ (main.inters_str).join([k[0] for k in comb]) ] = joinFlag

            return tmp

        if dflags:
            res = apply_cuts_combinations(dflags)
        else:
            res = {nocut_dummy_str: True}

        return res
