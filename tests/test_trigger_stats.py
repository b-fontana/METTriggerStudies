# coding: utf-8

_all_ = [ 'test_trigger_stats' ]

import os
import sys
parent_dir = os.path.abspath(__file__ + 2 * '/..')
sys.path.insert(0, parent_dir)
import argparse
import glob
import multiprocessing
import itertools as it

import inclusion
from inclusion import selection
from inclusion.config import main
from inclusion.utils import utils

import ROOT

def get_outname(suffix, mode, cut, ext):
    pref = {'met': 'MET', 'tau': 'Tau', 'met_tau': 'MET_and_Tau'}
    assert mode in list(pref.keys())
    if ext == 'csv':
        s = 'counts_{}_{}_{}/table.csv'.format(suffix, pref[mode], cut)
    elif ext == 'root':
        utils.create_single_dir('data')
        s = 'data/data_{}_{}_{}.{}'.format(suffix, pref[mode], cut, ext)
    else:
        raise ValueError('The {} extension is not supported.'.format(ext))
    return s

def set_plot_definitions():
    ROOT.gROOT.SetBatch(ROOT.kTRUE)
    ROOT.gStyle.SetOptStat(ROOT.kFALSE)
    ret = {'BoxTextSize'    : 50,
           'BoxTextFont'    : 43,
           'BoxTextColor'   : ROOT.kBlack,
           'XTitleSize'     : 0.045,
           'YTitleSize'     : 0.045,
           'LineWidth'      : 2,
           'FrameLineWidth' : 1,
           }
    return ret
    
def plot(mode, hbase, htrg, htrgcut, cut_strs, var, channel, sample, category, directory):
    legends = {'met': ['MET', 'MET + cut', 'met'],
               'tau': ['Tau', 'Tau + cut', 'tau'],
               'met_tau': ['MET', 'MET + cuts', 'met_and_tau', 'Tau', 'Tau + cuts']}

    assert mode in list(legends.keys())
    
    cat_dir = os.path.join(directory, sample, category)
    utils.create_single_dir(cat_dir)

    if isinstance(htrg, (tuple,list)):
        plot_three_histos(legends[mode], cut_strs,
                          hbase, htrg, htrgcut, var, channel, directory, cat_dir)
    else:
        plot_two_histos(legends[mode], cut_strs,
                        hbase, htrg, htrgcut, var, channel, directory, cat_dir)


def plot_three_histos(legends, cut_strs, hbase, htrg, htrgcut, var, channel, directory, cat_folder):
    defs = set_plot_definitions()
    
    hbase1 = hbase.Clone('hbase1')
    htrgcut1 = [ h.Clone('htrgcut1_' + str(ih)) for ih,h in enumerate(htrgcut) ]

    hbase2 = hbase.Clone('hbase2')
    htrg2 = [ h.Clone('htrg2_' + str(ih)) for ih,h in enumerate(htrg) ]
    htrgcut2 = [ h.Clone('htrg2_' + str(ih)) for ih,h in enumerate(htrgcut) ]

    c1 = ROOT.TCanvas('c1', '', 600, 400)
    c1.cd()
        
    # Absolute shapes    
    max_base = hbase1.GetMaximum() + (hbase1.GetMaximum()-hbase1.GetMinimum())/5.
    max_cut   = max(htrgcut1[0].GetMaximum(),htrgcut1[1].GetMaximum()) + (htrgcut1[0].GetMaximum()-htrgcut1[0].GetMinimum())/5.
    htrgcut1[0].SetMaximum( max(max_base, max_cut) )

    htrgcut1[1].GetXaxis().SetTitleSize(defs['XTitleSize']);
    htrgcut1[1].GetXaxis().SetTitle(var + ' [GeV]');
    htrgcut1[1].GetYaxis().SetTitleSize(defs['YTitleSize']);
    htrgcut1[1].GetYaxis().SetTitle('a. u.');
    htrgcut1[1].SetLineWidth(defs['LineWidth']);
    htrgcut1[1].SetLineColor(8);

    hbase1.SetLineWidth(1);
    hbase1.SetLineColor(1);
    htrgcut1[0].SetLineWidth(1);
    htrgcut1[0].SetLineColor(1);
    htrgcut1[1].SetLineWidth(1);
    htrgcut1[1].SetLineColor(1);

    htrgcut1[1].Add(hbase)
    htrgcut1[1].Add(htrgcut1[0])
    htrgcut1[1].SetFillColor(8)
    htrgcut1[1].Draw('hist')
    htrgcut1[0].Add(hbase)
    htrgcut1[0].SetFillColor(2)
    htrgcut1[0].Draw('histsame')
    hbase1.SetFillColor(4)
    hbase1.Draw('histsame')

    leg1 = ROOT.TLegend(0.65, 0.75, 0.90, 0.9)
    leg1.SetNColumns(1)
    leg1.SetFillStyle(0)
    leg1.SetBorderSize(0)
    leg1.SetTextFont(43)
    leg1.SetTextSize(10)
    leg1.AddEntry(htrgcut1[1], legends[3], 'F')
    leg1.AddEntry(htrgcut1[0], legends[0], 'F')
    leg1.AddEntry(hbase1, '+'.join(triggers[channel]), 'F')
    leg1.Draw('same')

    c1.Update();
    for ext in ('png', 'pdf'):
        c1.SaveAs( os.path.join(cat_folder, legends[2] + '_abs_' + var + '_' + cut_strs + '.' + ext) )
    c1.Close()

    # Normalized shapes
    c2 = ROOT.TCanvas('c2', '', 600, 400)
    c2.cd()
    try:
        for h in htrg2:
            h.Scale(1/h.Integral())
    except ZeroDivisionError:
        pass
    try:
        hbase2.Scale(1/hbase2.Integral())
    except ZeroDivisionError:
        pass
    try:
        for h in htrgcut2:
            h.Scale(1/h.Integral())
    except ZeroDivisionError:
        pass

    max_base2 = hbase2.GetMaximum() + (hbase2.GetMaximum()-hbase2.GetMinimum())/5.
    max_cut2   = (max(htrgcut2[0].GetMaximum(),htrgcut2[1].GetMaximum()) +
                  max((htrgcut2[0].GetMaximum()-htrgcut2[0].GetMinimum())/5., (htrgcut2[1].GetMaximum()-htrgcut2[1].GetMinimum())/5.))
    hbase2.SetMaximum( max(max_base2, max_cut2) )
    
    hbase2.GetXaxis().SetTitleSize(defs['XTitleSize']);
    hbase2.GetXaxis().SetTitle(var + ' [GeV]');
    hbase2.GetYaxis().SetTitleSize(defs['YTitleSize']);
    hbase2.GetYaxis().SetTitle('Normalized to 1');
    hbase2.SetLineWidth(defs['LineWidth']);
    hbase2.SetLineColor(4);
    htrgcut2[0].SetLineWidth(defs['LineWidth']);
    htrgcut2[0].SetLineColor(2);
    htrgcut2[1].SetLineWidth(defs['LineWidth']);
    htrgcut2[1].SetLineColor(8);

    #htrg2.Draw('hist')
    hbase2.Draw('hist')
    htrgcut2[0].Draw('histsame')
    htrgcut2[1].Draw('histsame')

    leg2 = ROOT.TLegend(0.69, 0.77, 0.90, 0.9)
    leg2.SetNColumns(1)
    leg2.SetFillStyle(0)
    leg2.SetBorderSize(0)
    leg2.SetTextFont(43)
    leg2.SetTextSize(10)
    #leg2.AddEntry(htrg2, legends[0])
    leg2.AddEntry(htrgcut2[0], legends[1])
    leg2.AddEntry(htrgcut2[1], legends[4])
    leg2.AddEntry(hbase2, '+\n'.join(triggers[channel]))
    leg2.Draw('same')
    
    c2.Update();
    for ext in ('png', 'pdf'):
        c2.SaveAs( os.path.join(cat_folder, legends[2] + '_norm_' + var + '_' + cut_strs + '.' + ext) )
    c2.Close()
    
def plot_two_histos(legends, cut_strs, hbase, htrg, htrgcut, var, channel, directory, cat_folder):
    defs = set_plot_definitions()
    
    hbase1 = hbase.Clone('hbase1')
    htrgcut1 = htrgcut.Clone('htrgcut1')

    htrg2 = htrg.Clone('htrg2')
    hbase2 = hbase.Clone('hbase2')
    htrgcut2 = htrgcut.Clone('htrgcut2')

    c1 = ROOT.TCanvas('c1', '', 600, 400)
    c1.cd()
        
    # Absolute shapes    
    max_base = hbase1.GetMaximum() + (hbase1.GetMaximum()-hbase1.GetMinimum())/5.
    max_cut   = htrgcut1.GetMaximum() + (htrgcut1.GetMaximum()-htrgcut1.GetMinimum())/5.
    htrgcut1.SetMaximum( max(max_base, max_cut) )

    htrgcut1.GetXaxis().SetTitleSize(defs['XTitleSize']);
    htrgcut1.GetXaxis().SetTitle(var + ' [GeV]');
    htrgcut1.GetYaxis().SetTitleSize(defs['YTitleSize']);
    htrgcut1.GetYaxis().SetTitle('a. u.');
    htrgcut1.SetLineWidth(defs['LineWidth']);
    htrgcut1.SetLineColor(8);

    hbase1.SetLineWidth(2);
    hbase1.SetLineColor(4);
    htrgcut1.SetLineWidth(2);
    htrgcut1.SetLineColor(2);

    #htrg.Draw('hist')
    htrgcut1.Add(hbase)
    htrgcut1.SetFillColor(2)
    htrgcut1.Draw('hist')
    hbase1.SetFillColor(4)
    hbase1.Draw('histsame')

    leg1 = ROOT.TLegend(0.69, 0.77, 0.90, 0.9)
    leg1.SetNColumns(1)
    leg1.SetFillStyle(0)
    leg1.SetBorderSize(0)
    leg1.SetTextFont(43)
    leg1.SetTextSize(10)
    leg1.AddEntry(htrgcut1, legends[0])
    leg1.AddEntry(hbase1, '+'.join(triggers[channel]))
    leg1.Draw('same')

    c1.Update();
    for ext in ('png', 'pdf'):
        c1.SaveAs( os.path.join(cat_folder, legends[2] + '_abs_' + var + '_' + cut_strs + '.' + ext) )
    c1.Close()

    # Normalized shapes
    c2 = ROOT.TCanvas('c2', '', 600, 400)
    c2.cd()
    try:
        htrg2.Scale(1/htrg2.Integral())
    except ZeroDivisionError:
        pass
    try:
        hbase2.Scale(1/hbase2.Integral())
    except ZeroDivisionError:
        pass
    try:
        htrgcut2.Scale(1/htrgcut2.Integral())
    except ZeroDivisionError:
        pass

    max_met2   = htrg2.GetMaximum() + (htrg2.GetMaximum()-htrg2.GetMinimum())/5.
    max_base2 = hbase2.GetMaximum() + (hbase2.GetMaximum()-hbase2.GetMinimum())/5.
    max_cut2   = htrgcut2.GetMaximum() + (htrgcut2.GetMaximum()-htrgcut2.GetMinimum())/5.
    htrg2.SetMaximum( max(max_met2, max_base2, max_cut2) )

    htrg2.GetXaxis().SetTitleSize(defs['XTitleSize']);
    htrg2.GetXaxis().SetTitle(var + ' [GeV]');
    htrg2.GetYaxis().SetTitleSize(defs['YTitleSize']);
    htrg2.GetYaxis().SetTitle('Normalized to 1');
    htrg2.SetLineWidth(defs['LineWidth']);
    htrg2.SetLineColor(8);

    hbase2.SetLineWidth(2);
    hbase2.SetLineColor(4);
    htrgcut2.SetLineWidth(2);
    htrgcut2.SetLineColor(2);

    htrg2.Draw('hist')
    hbase2.Draw('histsame')
    htrgcut2.Draw('histsame')

    leg2 = ROOT.TLegend(0.69, 0.77, 0.90, 0.9)
    leg2.SetNColumns(1)
    leg2.SetFillStyle(0)
    leg2.SetBorderSize(0)
    leg2.SetTextFont(43)
    leg2.SetTextSize(10)
    leg2.AddEntry(htrg2, legends[0])
    leg2.AddEntry(htrgcut2, legends[1])
    leg2.AddEntry(hbase2, '+\n'.join(triggers[channel]))
    leg2.Draw('same')
    
    c2.Update();
    for ext in ('png', 'pdf'):
        c2.SaveAs( os.path.join(cat_folder, legends[2] + '_norm_' + var + '_' + cut_strs + '.' + ext) )
    c2.Close()

        
def plot2D(mode, hbase, htrg, htrgcut, cut_strs, two_vars, channel, sample, category, directory):
    legends = {'met': ['MET', 'MET + cut', 'met'],
               'tau': ['Tau', 'Tau + cut', 'tau'],
               'met_tau': ['MET', 'MET + cuts', 'met_and_tau', 'Tau', 'Tau + cuts']}

    htrg1 = htrg.Clone('htrg1')
    hbase1 = hbase.Clone('hbase1')
    htrgcut1 = htrgcut.Clone('htrgcut1')

    defs = set_plot_definitions()    
    c = ROOT.TCanvas('c', '', 600, 400)

    c.cd()
    pad1 = ROOT.TPad('pad1', 'pad1', 0., 0., 0.333, 1.)
    pad1.SetFrameLineWidth(defs['FrameLineWidth'])
    pad1.SetLeftMargin(0.15);
    pad1.SetRightMargin(0.0);
    pad1.SetBottomMargin(0.08);
    pad1.SetTopMargin(0.055);
    pad1.Draw()
    pad1.cd()

    hbase.GetXaxis().SetTitle('');
    hbase.GetYaxis().SetTitle('');
    try:
        hbase.Scale(1/hbase.Integral())
    except ZeroDivisionError:
        pass
    hbase.Draw('colz');

    c.cd()
    pad2 = ROOT.TPad('pad2', 'pad2', 0.333, 0.0, 0.665, 1.0)
    pad2.SetFrameLineWidth(defs['FrameLineWidth'])
    pad2.SetLeftMargin(0.0);
    pad2.SetRightMargin(0.0);
    pad2.SetBottomMargin(0.08);
    pad2.SetTopMargin(0.055);
    pad2.Draw()
    pad2.cd()

    htrg1.GetXaxis().SetTitle('')
    htrg1.GetYaxis().SetTitle(two_vars[1])
    htrgcut1.GetYaxis().SetTitleSize(0.045)
    htrgcut1.GetYaxis().SetTitle(two_vars[1])
    try:
        htrg1.Scale(1/htrg1.Integral())
    except ZeroDivisionError:
        pass
    htrg1.Draw('colz');
 
    c.cd()
    pad3 = ROOT.TPad('pad3', 'pad3', 0.665, 0.0, 1.0, 1.0)
    pad3.SetFrameLineWidth(defs['FrameLineWidth'])
    pad3.SetLeftMargin(0.0);
    pad3.SetRightMargin(0.15);
    pad3.SetBottomMargin(0.08);
    pad3.SetTopMargin(0.055);
    pad3.Draw()
    pad3.cd()

    htrgcut1.GetXaxis().SetTitle(two_vars[0])
    htrgcut1.GetXaxis().SetTitleSize(0.045)
    try:
        htrgcut1.Scale(1/htrgcut1.Integral())
    except ZeroDivisionError:
        pass
    htrgcut1.Draw('colz')

    cat_folder = os.path.join(directory, sample, category)
    utils.create_single_dir(cat_folder)
    c.Update();
    for ext in ('png', 'pdf'):
        c.SaveAs( os.path.join(cat_folder, legends[mode][2] + '_' + '_VS_'.join(two_vars) + '_' + cut_strs + '.' + ext) )
    c.Close()

def count(mode, hbase, htrg, htrgcut, cut_strings, var, channel, sample, category, directory):
    cat_folder = os.path.join(directory, sample, category)
    titles = {'met': ['MET', 'MET + cut', 'Trigger baseline (no MET)', 'Fraction [%]: {[MET + Cut] / [Trigger baseline]} + 1\n'],
              'tau': ['Tau', 'Tau + cut', 'Trigger baseline (no Tau)', 'Fraction [%]: {[Tau + Cut] / [Trigger baseline]} + 1\n'],
              'met_tau': ['MET + Tau', 'MET + Tau + cut', 'Trigger baseline (no MET + Tau)', 'Fraction [%]: {[MET + Tau + Cut] / [Trigger baseline]} + 1\n'],}

    def calc_frac(c1, c2):
        try:
            frac = c2 / (c1 + c2)
        except ZeroDivisionError:
            frac = 0
        return frac * 100

    name_met = os.path.join(cat_folder,
                            get_outname(suffix=var, mode=mode, cut=cut_strings, ext='csv'))
    utils.create_single_dir(os.path.dirname(name_met))
    with open(name_met, 'w') as f:
        f.write(','.join(('Bin label', titles[mode][0], titles[mode][1], titles[mode][2], titles[mode][3])))
        for ibin in range(1, htrg.GetNbinsX()+1):
            label = str(round(htrg.GetXaxis().GetBinLowEdge(ibin),2)) + ' / ' + str(round(htrg.GetXaxis().GetBinLowEdge(ibin+1),2))
            cmet = htrg.GetBinContent(ibin)
            cnomet = hbase.GetBinContent(ibin)
            cmetcut = htrgcut.GetBinContent(ibin)
            assert cmet >= cmetcut

            frac  = calc_frac(cnomet, cmetcut)
            f.write(','.join((label, str(round(cmet,2)), str(round(cmetcut,2)), str(round(cnomet,2)), str(round(frac,2)))) + '\n')

        totmet = htrg.Integral(0,htrg.GetNbinsX()+1)
        totnomet = hbase.Integral(0,hbase.GetNbinsX()+1)
        totmetcut = htrgcut.Integral(0,htrgcut.GetNbinsX()+1)

        totfrac = calc_frac(totnomet, totmetcut)
        f.write(','.join(('Total', str(round(totmet,2)), str(round(totmetcut,2)), str(round(totnomet,2)), str(round(totfrac,2)))) + '\n')
    return totfrac

def counts_total(mode, totarr, channel, category, directory):
    titles = {'met': ['MET', 'MET + cut', 'Trigger baseline (no MET)', 'Fraction [%]: {[MET + Cut] / [Trigger baseline]} + 1\n'],
              'tau': ['Tau', 'Tau + cut', 'Trigger baseline (no Tau)', 'Fraction [%]: {[Tau + Cut] / [Trigger baseline]} + 1\n'],
              'met_tau': ['MET + Tau', 'MET + Tau + cut', 'Trigger baseline (no MET + Tau)', 'Fraction [%]: {[MET + Tau + Cut] / [Trigger baseline]} + 1\n'],}

    name = os.path.join(directory, 'counts_total_' + category + '_' + mode,  'table.csv')
    utils.create_single_dir(os.path.dirname(name))
    with open(name, 'w') as f:
        f.write(','.join(('Sample', 'Fraction')) + '\n')
        for sample, frac in totarr:
            f.write(','.join((sample, str(round(frac,3)))) + '\n')

def test_triger_stats(indir, sample, channel, plot_only, cut_strings):
    outname = get_outname(suffix=sample+'_'+channel, mode='met',
                          cut=cut_strings['met_tau'], ext='root')

    if channel == 'etau' or channel == 'mutau':
        iso1 = (24, 0, 8)
    elif channel == 'tautau':
        iso1 = binning['dau2_iso']
    binning.update({'HHKin_mass': (20, float(sample)-300, float(sample)+300),
                    'dau1_iso': iso1})
    
    full_sample = 'GluGluToBulkGravitonToHHTo2B2Tau_M-' + sample + '_'
    
    t_in = ROOT.TChain('HTauTauTree')
    glob_files = glob.glob( os.path.join(indir, full_sample, 'output_*.root') )
    for f in glob_files:
        t_in.Add(f)

    hBaseline = {}
    hMET, hMETWithCut = ({} for _ in range(2))
    hTau, hTauWithCut = ({} for _ in range(2))
    hTauNoMET, hTauNoMETWithCut = ({} for _ in range(2))
    hOR, hORWithCut = ({} for _ in range(2))
    for v in tuple(variables):
        hBaseline[v] = {}
        hMET[v], hMETWithCut[v] = ({} for _ in range(2))
        hTau[v], hTauWithCut[v] = ({} for _ in range(2))
        hTauNoMET[v], hTauNoMETWithCut[v] = ({} for _ in range(2))
        hOR[v], hORWithCut[v] = ({} for _ in range(2))
        for cat in categories:
            hBaseline[v][cat] = ROOT.TH1D('hBaseline_'+v+'_'+cat, '', *binning[v])
            hMET[v][cat] = ROOT.TH1D('hMET_'+v+'_'+cat, '', *binning[v])
            hMETWithCut[v][cat] = ROOT.TH1D('hMETWithCut_'+v+'_'+cat, '', *binning[v])
            hTau[v][cat] = ROOT.TH1D('hTau_'+v+'_'+cat, '', *binning[v])
            hTauWithCut[v][cat] = ROOT.TH1D('hTauWithCut_'+v+'_'+cat, '', *binning[v])
            hTauNoMET[v][cat] = ROOT.TH1D('hTauNoMET_'+v+'_'+cat, '', *binning[v])
            hTauNoMETWithCut[v][cat] = ROOT.TH1D('hTauNoMETWithCut_'+v+'_'+cat, '', *binning[v])
            hOR[v][cat] = ROOT.TH1D('hOR_'+v+'_'+cat, '', *binning[v])
            hORWithCut[v][cat] = ROOT.TH1D('hORWithCut_'+v+'_'+cat, '', *binning[v])

    hBaseline_2D = {}
    hMET_2D, hMETWithCut_2D = ({} for _ in range(2))
    hTau_2D, hTauWithCut_2D = ({} for _ in range(2))
    hOR_2D, hORWithCut_2D = ({} for _ in range(2))
    for v in variables_2D:
        hBaseline_2D[v] = {}
        hMET_2D[v], hMETWithCut_2D[v] = ({} for _ in range(2))
        hTau_2D[v], hTauWithCut_2D[v] = ({} for _ in range(2))
        hOR_2D[v], hORWithCut_2D[v] = ({} for _ in range(2))
        for cat in categories:
            hBaseline_2D[v][cat] = ROOT.TH2D('hBaseline_2D_'+'_'.join(v)+'_'+cat, '', *binning[v[0]], *binning[v[1]])
            hMET_2D[v][cat] = ROOT.TH2D('hMET_2D_'+'_'.join(v)+'_'+cat, '', *binning[v[0]], *binning[v[1]])
            hMETWithCut_2D[v][cat] = ROOT.TH2D('hMETWithCut_2D_'+'_'.join(v)+'_'+cat, '', *binning[v[0]], *binning[v[1]])        
            hTau_2D[v][cat] = ROOT.TH2D('hTau_2D_'+'_'.join(v)+'_'+cat, '', *binning[v[0]], *binning[v[1]])
            hTauWithCut_2D[v][cat] = ROOT.TH2D('hTauWithCut_2D_'+'_'.join(v)+'_'+cat, '', *binning[v[0]], *binning[v[1]])        
            hOR_2D[v][cat] = ROOT.TH2D('hOR_2D_'+'_'.join(v)+'_'+cat, '', *binning[v[0]], *binning[v[1]])
            hORWithCut_2D[v][cat] = ROOT.TH2D('hORWithCut_2D_'+'_'.join(v)+'_'+cat, '', *binning[v[0]], *binning[v[1]])        

    t_in.SetBranchStatus('*', 0)
    _entries = ('triggerbit', 'RunNumber', 'isLeptrigger',
                'bjet1_bID_deepFlavor', 'bjet2_bID_deepFlavor', 'isBoosted',
                'isVBF', 'VBFjj_mass', 'VBFjj_deltaEta', 'PUReweight', 'lumi', 'IdAndIsoSF_deep_pt',
                'pairType', 'dau1_eleMVAiso', 'dau1_iso', 'dau1_deepTauVsJet', 'dau2_deepTauVsJet',
                'nleps', 'nbjetscand', 'tauH_SVFIT_mass', 'bH_mass_raw',)
    _entries += tuple(variables)
    for ientry in _entries:
        t_in.SetBranchStatus(ientry, 1)
  
    for entry in t_in:
        # this is slow: do it once only
        entries = utils.dot_dict({x: getattr(entry, x) for x in _entries})
        sel = selection.EventSelection(entries, isdata=False, configuration=None)
        
        # mcweight   = entries.MC_weight
        pureweight = entries.PUReweight
        lumi       = entries.lumi
        idandiso   = entries.IdAndIsoSF_deep_pt
        
        #if utils.is_nan(mcweight)  : mcweight=1
        if utils.is_nan(pureweight) : pureweight=1
        if utils.is_nan(lumi)       : lumi=1
        if utils.is_nan(idandiso)   : idandiso=1
  
        evt_weight = pureweight*lumi*idandiso
        if utils.is_nan(evt_weight):
            evt_weight = 1
  
        if utils.is_channel_consistent(channel, entries.pairType):
            if not sel.selection_cuts(lepton_veto=True, bjets_cut=True,
                                      standard_mass_cut=True, invert_mass_cut=False):
                continue
  
            for v in variables:
                for cat in categories:
                    if sel.sel_category(cat):
  
                        # passes the OR of the trigger baseline (not including METNoMu120 trigger)
                        pass_trg = sel.pass_triggers(triggers[channel]) and entries.isLeptrigger
                        if pass_trg:
                            hBaseline[v][cat].Fill(entries[v], evt_weight)

                        met_cut_expr = entries.metnomu_et > met_cut
                        tau_cut_expr = ((entries.dau1_pt > tau_cut and args.channel=='tautau') or
                                        (entries.dau2_pt > tau_cut and args.channel!='tautau'))

                        # passes the METNoMu120 trigger and does *not* pass the OR of the baseline
                        if not pass_trg and eval(' '.join(args.custom_cut)):
                            if sel.pass_triggers(('METNoMu120',)):
                                hMET[v][cat].Fill(entries[v], evt_weight)
                                if met_cut_expr:
                                    hMETWithCut[v][cat].Fill(entries[v], evt_weight)
                            
                            # passes the IsoTau180 trigger and does *not* pass the OR of the baseline
                            if sel.pass_triggers(('IsoTau180',)):
                                hTau[v][cat].Fill(entries[v], evt_weight)
                                if tau_cut_expr:
                                    hTauWithCut[v][cat].Fill(entries[v], evt_weight)

                            # passes the IsoTau180 trigger and does *not* pass the OR of the baseline and METNoMu120
                            if sel.pass_triggers(('IsoTau180',)) and not sel.pass_triggers(('METNoMu120',)):
                                hTauNoMET[v][cat].Fill(entries[v], evt_weight)
                                if tau_cut_expr:
                                    hTauNoMETWithCut[v][cat].Fill(entries[v], evt_weight)

                            # passes the METNoMu120 or the IsoTau180 triggers and does *not* pass the OR of the baseline
                            if sel.pass_triggers(('METNoMu120', 'IsoTau180',)):
                                hOR[v][cat].Fill(entries[v], evt_weight)
                                if ((sel.pass_triggers(('METNoMu120',)) and met_cut_expr) or
                                    (sel.pass_triggers(('IsoTau180',)) and tau_cut_expr)):
                                    hORWithCut[v][cat].Fill(entries[v], evt_weight)

            for v in variables_2D:
                for cat in categories:
                    if sel.sel_category(cat):
  
                        # passes the OR of the trigger baseline (not including METNoMu120 trigger)
                        pass_trg = sel.pass_triggers(triggers[channel]) and entries.isLeptrigger
                        if pass_trg:
                            hBaseline_2D[v][cat].Fill(entries[v[0]], entries[v[1]], evt_weight)
  
                        # passes the METNoMu120 trigger and does *not* pass the OR of the baseline
                        if not pass_trg and eval(' '.join(args.custom_cut)):
                            if sel.pass_triggers(('METNoMu120',)):
                                hMET_2D[v][cat].Fill(entries[v[0]], entries[v[1]], evt_weight)
                                if met_cut_expr:
                                    hMETWithCut_2D[v][cat].Fill(entries[v[0]], entries[v[1]], evt_weight)

                            if sel.pass_triggers(('IsoTau180',)):
                                hTau_2D[v][cat].Fill(entries[v[0]], entries[v[1]], evt_weight)
                                if met_cut_expr:
                                    hTauWithCut_2D[v][cat].Fill(entries[v[0]], entries[v[1]], evt_weight)

                            if sel.pass_triggers(('METNoMu120', 'IsoTau180',)):
                                hOR_2D[v][cat].Fill(entries[v[0]], entries[v[1]], evt_weight)
                                if ((sel.pass_triggers(('METNoMu120',)) and met_cut_expr) or
                                    (sel.pass_triggers(('IsoTau180',)) and tau_cut_expr)):
                                    hORWithCut_2D[v][cat].Fill(entries[v[0]], entries[v[1]], evt_weight)


    f_out = ROOT.TFile(outname, 'RECREATE')
    f_out.cd()
    for cat in categories:
        for v in variables:
            hBaseline[v][cat].Write('hBaseline_' + v + '_' + cat)
            hMET[v][cat].Write('hMET_' + v + '_' + cat)
            hMETWithCut[v][cat].Write('hMETWithCut_' + v + '_' + cat)
            hTau[v][cat].Write('hTau_' + v + '_' + cat)
            hTauWithCut[v][cat].Write('hTauWithCut_' + v + '_' + cat)
            hTauNoMET[v][cat].Write('hTauNoMET_' + v + '_' + cat)
            hTauNoMETWithCut[v][cat].Write('hTauNoMETWithCut_' + v + '_' + cat)
            hOR[v][cat].Write('hOR_' + v + '_' + cat)
            hORWithCut[v][cat].Write('hORWithCut_' + v + '_' + cat)
        for v in variables_2D:
            hBaseline_2D[v][cat].Write('hBaseline_2D_' + '_'.join(v)+'_'+ cat)
            hMET_2D[v][cat].Write('hMET_2D_' + '_'.join(v)+'_'+ cat)
            hMETWithCut_2D[v][cat].Write('hMETWithCut_2D_' + '_'.join(v)+'_'+ cat)
            hTau_2D[v][cat].Write('hTau_2D_' + '_'.join(v)+'_'+ cat)
            hTauWithCut_2D[v][cat].Write('hTauWithCut_2D_' + '_'.join(v)+'_'+ cat)
            hOR_2D[v][cat].Write('hOR_2D_' + '_'.join(v)+'_'+ cat)
            hORWithCut_2D[v][cat].Write('hORWithCut_2D_' + '_'.join(v)+'_'+ cat)
    f_out.Close()
    print('Raw histograms saved in {}.'.format(outname), flush=True)

if __name__ == '__main__':
    triggers = {'etau': ('Ele32', 'EleIsoTauCustom'),
                'mutau': ('IsoMu24', 'IsoMuIsoTauCustom'),
                'tautau': ('IsoDoubleTauCustom',)}
    binning = {'metnomu_et': (20, 0, 450),
               'dau1_pt': (30, 0, 450),
               'dau1_eta': (20, -2.5, 2.5),
               'dau2_iso': (20, 0.88, 1.01),
               'dau2_pt': (30, 0, 350),
               'dau2_eta': (20, -2.5, 2.5),
               'ditau_deltaR': (30, 0.3, 1.3),
               'dib_deltaR': (25, 0, 2.5),
               'bH_pt': (20, 70, 600),
               'bH_mass': (30, 0, 280),
               'tauH_mass': (30, 0, 170),
               'tauH_pt': (30, 0, 500),
               'tauH_SVFIT_mass': (30, 0, 250),
               'tauH_SVFIT_pt': (20, 200, 650),
               'bjet1_pt': (25, 10, 600),
               'bjet2_pt': (25, 10, 550),
               'bjet1_eta': (20, -2.5, 2.5),
               'bjet2_eta': (20, -2.5, 2.5),
               }
    variables = tuple(binning.keys()) + ('HHKin_mass', 'dau1_iso')
    #variables = tuple(('dau1_pt', 'dau2_pt', 'metnomu_et', 'HHKin_mass', 'dau1_iso', 'dau2_iso'))
    variables_2D = (('dau1_pt', 'dau2_pt'),)#('dau1_iso', 'dau2_iso'))

    #categories = ('baseline', 's1b1jresolvedMcut', 's2b0jresolvedMcut', 'sboostedLLMcut')
    categories = ('baseline',)
    met_cut = 200
    tau_cut = 190
    
    # Parse input arguments
    parser = argparse.ArgumentParser(description='Producer trigger histograms.')

    parser.add_argument('--indir', required=True, type=str,
                        help='Full path of ROOT input file')
    parser.add_argument('--samples', required=True, nargs='+', type=str,
                        help='Full path of ROOT input file')
    parser.add_argument('--channel', required=True, type=str,  
                        help='Select the channel over which the workflow will be run.' )
    parser.add_argument('--plot_only', action='store_true',
                        help='Reuse previously produced data for quick plot changes.')
    parser.add_argument('--plot_2D_only', action='store_true',
                        help='Reuse previously produced data for quick plot changes.')
    parser.add_argument('--custom_cut', nargs='+', default=['True'],
                        help='Customisable cut provided by the user.')
    parser.add_argument('--no_copy', action='store_true',
                        help='Do not copy the outputs to EOS at the end.')
    parser.add_argument('--sequential', action='store_true',
                        help='Do not use the multiprocess package.')
    args = utils.parse_args(parser)
            
    main_dir = 'TriggerStudy_MET'+str(met_cut)+'_SingleTau'+str(tau_cut)
    cut_expr = ''
    if '_'.join(args.custom_cut) != 'True':
        cut_expr = '_CUT_' + '_'.join(args.custom_cut)
        cut_expr = cut_expr.replace('.','_').replace(' ','_').replace('>','GT').replace('<','ST')
        main_dir += cut_expr

    cut_strings = {'met': str(met_cut) + cut_expr,
                   'tau': str(tau_cut) + cut_expr,
                   'met_tau': str(met_cut) + '_' + str(tau_cut) + cut_expr}

    #### run major function ###
    if args.sequential:
        if not args.plot_only and not args.plot_2D_only:
            for sample in args.samples:
                test_triger_stats(args.indir, sample, args.channel, args.plot_only, cut_strings)
    else:
        if not args.plot_only and not args.plot_2D_only:
            pool = multiprocessing.Pool(processes=6)
            pool.starmap(test_triger_stats, zip(it.repeat(args.indir), args.samples,
                                                it.repeat(args.channel), it.repeat(args.plot_only),
                                                it.repeat(cut_strings)))
    ###########################

    from_directory = os.path.join(main_dir, args.channel)
    totcounts = {'met': {}, 'tau': {}, 'met_tau': {}}
    for cat in categories:
        totcounts['met'][cat] = []
        totcounts['tau'][cat] = []
        totcounts['met_tau'][cat] = []
        
    for sample in args.samples:
        outname = get_outname(suffix=sample+'_'+args.channel, mode='met',
                              cut=cut_strings['met_tau'], ext='root')
        f_in = ROOT.TFile(outname, 'READ')
        f_in.cd()
        for cat in categories:
            if not args.plot_2D_only:
                for v in variables:
                    suff = lambda x : x + v + '_' + cat
                    
                    hBaseline = f_in.Get(suff('hBaseline_'))
                    hMET = f_in.Get(suff('hMET_'))
                    hMETWithCut = f_in.Get(suff('hMETWithCut_'))
                    hTau = f_in.Get(suff('hTau_'))
                    hTauWithCut = f_in.Get(suff('hTauWithCut_'))
                    hTauNoMET = f_in.Get(suff('hTauNoMET_'))
                    hTauNoMETWithCut = f_in.Get(suff('hTauNoMETWithCut_'))
                    hOR = f_in.Get(suff('hOR_'))
                    hORWithCut = f_in.Get(suff('hORWithCut_'))

                    hBaseline_c = hBaseline.Clone(suff('hBaseline_') + '_c')
                    hMET_c = hMET.Clone(suff('hMET_') + '_c')
                    hMETWithCut_c = hMETWithCut.Clone(suff('hMETWithCut_') + '_c')
                    hTau_c = hTau.Clone(suff('hTau_') + '_c')
                    hTauWithCut_c = hTauWithCut.Clone(suff('hTauWithCut_') + '_c')
                    hTauNoMET_c = hTauNoMET.Clone(suff('hTauNoMET_') + '_c')
                    hTauNoMETWithCut_c = hTauNoMETWithCut.Clone(suff('hTauNoMETWithCut_') + '_c')
                    # hOR_c = hOR.Clone(suff('hOR_') + '_c')
                    # hORWithCut_c = hORWithCut.Clone(suff('hORWithCut_') + '_c')

                    hOverlayBaseline_c = hBaseline.Clone(suff('hOverlayBaseline_') + '_c')
                    hOverlayMET_c = hMET.Clone(suff('hOverlayMET_') + '_c')
                    hOverlayBaseline_c.Add(hOverlayMET_c)
                    
                    opt = (v, args.channel, sample, cat, from_directory)

                    plot('met', hBaseline, hMET, hMETWithCut, cut_strings['met'], *opt)
                    plot('tau', hBaseline, hTau, hTauWithCut, cut_strings['tau'], *opt)
                    plot('met_tau', hBaseline, [hMET,hTauNoMET], [hMETWithCut,hTauNoMETWithCut],
                         cut_strings['met_tau'], *opt)

                    c1 = count('met', hBaseline_c, hMET_c, hMETWithCut_c, cut_strings['met'], *opt)
                    c2 = count('tau', hBaseline_c, hTau_c, hTauWithCut_c, cut_strings['tau'], *opt)
                    c3 = count('met_tau',  hBaseline_c, hTauNoMET_c, hTauNoMETWithCut_c,
                               cut_strings['met_tau'],*opt)
                    c4 = count('met_tau', hOverlayBaseline_c, hTauNoMET_c, hTauNoMETWithCut_c,
                               cut_strings['met_tau'], *opt)
                    assert c3 <= c2
                    assert c4 <= c2
                    assert c4 <= c3

                totcounts['met'][cat].append((sample, c1))
                totcounts['tau'][cat].append((sample, c2))
                totcounts['met_tau'][cat].append((sample, c4))
                    
            for v in variables_2D:
                opt_2D = (v, args.channel, sample, cat, from_directory)
                hBaseline_2D = f_in.Get('hBaseline_2D_' + '_'.join(v)+'_'+ cat)
                hMET_2D = f_in.Get('hMET_2D_' + '_'.join(v)+'_'+ cat)
                hMETWithCut_2D = f_in.Get('hMETWithCut_2D_' + '_'.join(v)+'_'+ cat)
                plot2D('met', hBaseline_2D, hMET_2D, hMETWithCut_2D, cut_strings['met'], *opt_2D)
                hTau_2D = f_in.Get('hTau_2D_' + '_'.join(v)+'_'+ cat)
                hTauWithCut_2D = f_in.Get('hTauWithCut_2D_' + '_'.join(v)+'_'+ cat)
                plot2D('tau', hBaseline_2D, hTau_2D, hTauWithCut_2D, cut_strings['tau'], *opt_2D)
                hOR_2D = f_in.Get('hOR_2D_' + '_'.join(v)+'_'+ cat)
                hORWithCut_2D = f_in.Get('hORWithCut_2D_' + '_'.join(v)+'_'+ cat)
                plot2D('met_tau', hBaseline_2D, hOR_2D, hORWithCut_2D, cut_strings['met_tau'], *opt_2D)
                    
        f_in.Close()

    for cat in categories:
        opt2 = (args.channel, cat, from_directory)
        counts_total('met', totcounts['met'][cat], *opt2)
        counts_total('tau', totcounts['tau'][cat], *opt2)
        counts_total('met_tau', totcounts['met_tau'][cat], *opt2)

    if not args.no_copy:
        import subprocess
        to_directory = os.path.join('/eos/user/b/bfontana/www/TriggerScaleFactors', main_dir)
        to_directory = os.path.join(to_directory, args.channel)

        for m in ('met', 'tau', 'met_tau'):
            for cat in categories:
                folder_name = 'counts_total_' + cat + '_' + m
                folder_to = os.path.join(to_directory, folder_name)
                utils.create_single_dir(folder_to)
                counts_files = os.path.join(from_directory, folder_name, 'table.csv')
                print('Copying: {}\t\t--->\t{}'.format(counts_files, folder_to), flush=True)
                subprocess.run(['rsync', '-ah', counts_files, os.path.join(folder_to, 'table.csv')])

        for sample in args.samples:
            sample_from = os.path.join(from_directory, sample)
            print('Copying: {}\t\t--->\t{}'.format(sample_from, to_directory), flush=True)
            subprocess.run(['rsync', '-ah', sample_from, to_directory])

    print('Done.')
