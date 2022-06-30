from collections import defaultdict
from luigi_conf import (
    _cuts,
    _cuts_ignored,
    _data,
    _triggers_map as tmap,
    _triggers_custom as tcust,
)

class EventSelection:
    def __init__(self, leaf_manager, dataset, isdata, debug=False):
        self.bit = leaf_manager.get_leaf('pass_triggerbit')
        self.run = leaf_manager.get_leaf('RunNumber')
        self.lm = leaf_manager
        self.isdata = isdata
        self.debug = debug

        self.dataset = dataset
        if self.dataset not in _data:
            mes = 'Dataset {} is not supported.'.format(self.dataset)
            raise ValueError(mes)
        
    def check_bit(self, bitpos):
        bitdigit = 1
        res = bool(self.bit&(bitdigit<<bitpos))
        return res

    def _pass_triggers(self, trigs):
        """
        Internal only function.
        Checks at least one trigger was fired.
        """
        flag = False
        for trig in trigs:
            if trig in tcust:
                flag = set_custom_trigger_bit(trig, self.bit, self.run)
            else:
                flag = self.check_bit(self.get_trigger_bit(trig))
            if flag:
                return True
        return False
    
    def any_trigger(self, trigs):
        """
        Checks at least one trigger was fired.
        Considers all framework triggers.
        """
        return self._pass_triggers(trigs)

    def dataset_triggers(self, trigs):
        """
        Checks at least one trigger was fired.
        Considers framework triggers for a specific dataset.
        """
        dataset_ref_trigs = { 'MET': ('METNoMu120',), # triggers for the MET dataset
                              'EG':  ('Ele32',) } # triggers for the EGamma dataset
        for k in _data:
            if k not in dataset_ref_trigs:
                mes = 'Specify the reference trigger for dataset {}!'
                raise ValueError(mes.format(self.dataset))
        for vals in dataset_ref_trigs.values():
            for v in vals:
                if v not in trigs:
                    mes = 'Reference trigger {} is not part of triggers {}.'
                    raise ValueError(mes.format(v,trigs))

        return self._pass_triggers(dataset_ref_trigs[self.dataset])

    def trigger_bits(self, trig):
        if trig in tcust:
            return self.set_custom_trigger_bit(trig)
        else:
            return self.check_bit(self.get_trigger_bit(trig))

    def selection_cuts(self, iso_cuts=dict(), lepton_veto=True,
                       bjets_cut=True, invert_mass_cut=True):
        """
        Applies selection cut to one event.
        Returns `True` only if all selection cuts pass.
        """
        mhh = self.lm.get_leaf('HHKin_mass')
        if mhh < 1:
            return False

        pairtype = self.lm.get_leaf('pairType')
        dau1_eleiso = self.lm.get_leaf('dau1_eleMVAiso')
        dau1_muiso  = self.lm.get_leaf('dau1_iso')
        dau1_tauiso = self.lm.get_leaf('dau1_deepTauVsJet')
        dau2_tauiso = self.lm.get_leaf('dau2_deepTauVsJet')

        # third lepton veto
        nleps = self.lm.get_leaf('nleps')
        if nleps > 0 and lepton_veto:
            return False

        # require at least two b jet candidates
        nbjetscand = self.lm.get_leaf('nbjetscand')
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
        svfit_mass = self.lm.get_leaf('tauH_SVFIT_mass')
        bh_mass    = self.lm.get_leaf('bH_mass_raw')

        mcut = ( (svfit_mass-129.)*(svfit_mass-129.) / (53.*53.) +
                 (bh_mass-169.)*(bh_mass-169.) / (145.*145.) ) <  1.0
        if mcut and invert_mass_cut: # inverted elliptical mass cut (-> ttCR)
            return False

        return True

    def dataset_cuts(self):
        """
        Applies selection cuts depending on the reference trigger being considered.
        Reference triggers tend to be associated with a certain dataset.
        For instance, the 'MET' dataset is connected to the MET Trigger.
        """
        if self.dataset == 'MET':
            return self.selection_cuts(lepton_veto=True)
        elif self.dataset == 'EG':
            return self.selection_cuts(lepton_veto=True)
        else:
            if self.dataset in _data.keys():
                mes = 'You forgot to include dataset {} in EventSelection!'
            else:
                mes = 'Dataset {} is not supported.'
            raise ValueError(mes.format(self.dataset))

    def set_custom_trigger_bit(self, trigger):
        """
        The VBF trigger was updated during data taking, adding HPS
        https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauTrigger
        """
        if trigger not in tcust:
            import inspect
            currentFunction = inspect.getframeinfo(frame).function
            raise ValueError('[{}] option not supported.'.format(currentFunction))

        if self.run < 317509 and self.isdata:
            if trigger == 'VBFTauCustom':
                bits = self.check_bit(tmap[trigger]['VBFTau']['data'])
            elif trigger == 'IsoDoubleTauCustom':
                bits = ( self.check_bit(tmap[trigger]['IsoDoubleTau']['data'][0]) or
                         self.check_bit(tmap[trigger]['IsoDoubleTau']['data'][1]) or
                         self.check_bit(tmap[trigger]['IsoDoubleTau']['data'][2]) )
            elif trigger == 'IsoMuIsoTauCustom':
                bits = self.check_bit(tmap[trigger]['IsoMuIsoTau']['data'])
            elif trigger == 'EleIsoTauCustom':
                bits = self.check_bit(tmap[trigger]['EleIsoTau']['data'])

        else:
            s = 'data' if self.isdata else 'mc'
            if trigger == 'VBFTauCustom':
                bits = self.check_bit(tmap[trigger]['VBFTauHPS'][s])
            elif trigger == 'IsoDoubleTauCustom':
                bits = self.check_bit(tmap[trigger]['IsoDoubleTauHPS'][s])
            elif trigger == 'IsoMuIsoTauCustom':
                bits = self.check_bit(tmap[trigger]['IsoMuIsoTauHPS'][s])
            elif trigger == 'EleIsoTauCustom':
                bits = self.check_bit(tmap[trigger]['EleIsoTauHPS'][s])

        return bits

    def get_trigger_bit(self, trigger_name):
        """
        Returns the trigger bit corresponding to 'tmap'
        """
        s = 'data' if self.isdata else 'mc'
        res = tmap[trigger_name]
        try:
            res = res[s]
        except KeyError:
            print('You likely forgot to add your custom trigger to tcust.')
            raise
        return res

    def match_inters_with_dataset(self, tcomb, channel):
        """Matches a trigger intersection with a reference trigger."""
        if self.dataset not in _data.keys():
            raise ValueError('Dataset {} is not supported.'.format(self.dataset))
        
        # Ignore long intersections for simplicity CHANGE !!!!!!! ??
        if len(tcomb) > 3:
            return False
     
        wrong_comb = 'Combination {} is not supported for channel {}.'
     
        # one single reference per trigger combination
        decision = lambda ds : True if self.dataset == ds else False
     
        #general triggers
        if tcomb in ( ('IsoTau180',),
                      ('VBFTauCustom'),
                      ('METNoMu120'),
                      ('IsoTau180', 'METNoMu120'),
                      ('IsoTau180', 'VBFTauCustom'),
                      ('METNoMu120', 'VBFTauCustom'),
                      ('IsoTau180', 'METNoMu120', 'VBFTauCustom'),
                     ):
            return decision('MET')
        elif tcomb in ( ):
            return decision('EG')
     
        # channel-specific triggers
        if channel == 'etau':
            if tcomb in ( ('Ele32',),
                          ('EleIsoTauCustom',),
                          ('Ele32', 'EleIsoTauCustom'),
                          ('Ele32', 'VBFTauCustom'),
                          ('Ele32', 'METNoMu120'),
                          ('Ele32', 'IsoTau180'),
                          ('EleIsoTauCustom', 'VBFTauCustom'),
                          ('EleIsoTauCustom', 'METNoMu120'),
                          ('EleIsoTauCustom', 'IsoTau180'),
                          ('EleIsoTauCustom', 'IsoTau180'),
                          ('Ele32', 'EleIsoTauCustom', 'VBFTauCustom'),
                          ('Ele32', 'EleIsoTauCustom', 'METNoMu120'),
                          ('Ele32', 'EleIsoTauCustom', 'IsoTau180'),
                          ('Ele32', 'IsoTau180', 'VBFTauCustom'),
                         ):
                return decision('MET')
            elif tcomb in ( ) :
                return decision('EG')
            else:
                raise ValueError(wrong_comb.format(tcomb, channel))
            
        elif channel == 'mutau':
            if tcomb in ( ('IsoMu24',),
                          ('IsoMuIsoTauCustom',),
                          ('IsoMu24', 'IsoMuIsoTauCustom'),
                          ('IsoMu24', 'VBFTauCustom'),
                          ('IsoMu24', 'METNoMu120'),
                          ('IsoMu24', 'IsoTau180'),
                          ('IsoMuIsoTauCustom', 'VBFTauCustom'),
                          ('IsoMuIsoTauCustom', 'METNoMu120'),
                          ('IsoMuIsoTauCustom', 'IsoTau180'),
                          ('IsoMuIsoTauCustom', 'IsoTau180'),
                          ('IsoMu24', 'IsoMuIsoTauCustom', 'VBFTauCustom'),
                          ('IsoMu24', 'IsoMuIsoTauCustom', 'METNoMu120'),
                          ('IsoMu24', 'IsoMuIsoTauCustom', 'IsoTau180'),
                          ('IsoMu24', 'IsoTau180', 'VBFTauCustom'),
                         ):
                return decision('MET')
            elif tcomb in ( ) :
                return decision('EG')
            else:
                raise ValueError(wrong_comb.format(tcomb, channel))
            
        elif channel == 'tautau':
            if tcomb in ( ('IsoDoubleTauCustom',),
                          ('IsoDoubleTauCustom', 'VBFTauCustom'),
                          ('IsoDoubleTauCustom', 'METNoMu120'),
                          ('IsoDoubleTauCustom', 'IsoTau180'),
                          ('IsoDoubleTauCustom', 'IsoTau180', 'VBFTauCustom'),
                         ):
                return decision('MET')
            elif tcomb in ( ) :
                return decision('EG')
            else:
                raise ValueError(wrong_comb.format(tcomb, channel))
            
        else:
            raise ValueError('Channel {} is not supported.'.format(channel))

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
            trig_cuts = _cuts[trig]
        except KeyError: # the trigger has no cut associated
            if self.debug:
                print('KeyError')            
            return {nocut_dummy_str: True}

        for avar,acut in trig_cuts.items():
            # ignore cuts according to the user's definition in '_cuts_ignored'
            # example: do not cut on 'met_et' when displaying 'metnomu_et'
            ignore = functools.reduce( lambda x, y: x or y,
                                       [ avar in _cuts_ignored[k] for k in variables
                                         if k in _cuts_ignored ],
                                      False
                                      )

            # additionally, by default do not cut on the variable(s) being plotted
            if avar not in variables and not ignore:
                value = self.lm.get_leaf(avar) 

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
