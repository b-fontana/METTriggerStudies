# coding: utf-8

__all__ = ['TriggerBitsNumber']

import unittest

import os
import sys
parent_dir = os.path.abspath(__file__ + 2 * '/..')
sys.path.insert(0, parent_dir)

import inclusion
from inclusion import selection
from inclusion import config

import ROOT

class TriggerBitsNumber(unittest.TestCase):
    def setUp(self):
        self.isdata = True
        self.entry_names = ('triggerbit', 'RunNumber', 'PUReweight', 'lumi',
                            'IdAndIsoSF_deep_pt', 'HHKin_mass', 'pairType',
                            'dau1_eleMVAiso', 'dau1_iso', 'dau1_deepTauVsJet',
                            'dau2_deepTauVsJet', 'nleps', 'nbjetscand',
                            'tauH_SVFIT_mass', 'bH_mass_raw',)
        self.entry_names += tuple(config.var_eff)
        self.dummy_dataset = 'MET' # not used but required by the EventSelection class

    def my_print(self, i, n):
        print('\r{} / {}'.format(i, n), flush=True, end='')
        
    def test_met_bits(self):
        infile = '/data_CMS/cms/alves/HHresonant_SKIMS/SKIMS_UL18_Aug15Evening/MET__Run2018A/output_10.root'
        f_in = ROOT.TFile(infile)
        t_in = f_in.Get('HTauTauTree')
     
        t_in.SetBranchStatus('*', 0)
        for ientry in self.entry_names:
            t_in.SetBranchStatus(ientry, 1)
            
        cmet = 0
        nentries = t_in.GetEntriesFast()
        for ientry,entry in enumerate(t_in):
            if ientry%1000==0:
                self.my_print(ientry, nentries)

            entries = {x: getattr(entry, x) for x in self.entry_names}
            sel = selection.EventSelection(entries, self.dummy_dataset, self.isdata)
     
            if sel.check_bit(sel.get_trigger_bit('METNoMu120')):
                cmet += 1

        return cmet > 50 # some tunable cut

    def test_muon_bits(self):
        infile = '/data_CMS/cms/alves/HHresonant_SKIMS/SKIMS_UL18_Aug15Evening/SingleMuon__Run2018A/output_10.root'
        f_in = ROOT.TFile(infile)
        t_in = f_in.Get('HTauTauTree')
     
        t_in.SetBranchStatus('*', 0)
        for ientry in self.entry_names:
            t_in.SetBranchStatus(ientry, 1)
            
        cmuon = 0
        nentries = t_in.GetEntriesFast()
        for ientry,entry in enumerate(t_in):
            if ientry%1000==0:
                self.my_print(ientry, nentries)

            entries = {x: getattr(entry, x) for x in self.entry_names}
            sel = selection.EventSelection(entries, self.dummy_dataset, self.isdata)
     
            if sel.check_bit(sel.get_trigger_bit('IsoMu24')):
                cmuon += 1

        return cmuon > 50 # some tunable cut
