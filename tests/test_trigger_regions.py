# coding: utf-8

_all_ = [ 'test_trigger_regions' ]

import os
import sys
parent_dir = os.path.abspath(__file__ + 2 * '/..')
sys.path.insert(0, parent_dir)
import argparse
import glob
import multiprocessing
import itertools as it
import csv
from collections import defaultdict

import inclusion
from inclusion import selection
from inclusion.config import main
from inclusion.utils import utils

import ROOT

import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn3_circles

def square_diagram(mass, c_ditau_trg, c_met_trg, c_tau_trg,
                   region_cuts, pt_cuts, output_name):
    from bokeh.plotting import figure, output_file, save
    from bokeh.models import Range1d, Label

    tau = '\u03C4'
    ditau = tau+tau

    output_file(output_name)
    
    p = figure(title='m(X)={}GeV'.format(mass), width=600, height=400)
    p.x_range = Range1d(0, 12)
    p.y_range = Range1d(0, 10)
    p.outline_line_color = None
    p.xgrid.grid_line_color = None
    p.ygrid.grid_line_color = None
     
    # add a square renderer with a size, color, and alpha
    polyg_opt = dict(alpha=0.3)
    p.multi_polygons(xs=[[[[2.6, 9.9, 9.9, 2.6]]]],
                     ys=[[[[9.9, 9.9, 2.6, 2.6]]]], color='green', legend_label=ditau, **polyg_opt)
    p.multi_polygons(xs=[[[[0.1, 1.9, 1.9, 5.4, 5.4, 0.1]]]],
                     ys=[[[[5.4, 5.4, 1.9, 1.9, 0.1, 0.1]]]], color='blue', legend_label='MET', **polyg_opt)
    p.multi_polygons(xs=[[[[0.1, 0.1, 1.9, 1.9], [5.6, 5.6, 9.9, 9.9]]]],
                     ys=[[[[5.6, 9.9, 9.9, 5.6], [0.1, 1.9, 1.9, 0.1]]]], legend_label=tau, color='red', **polyg_opt)
    p.legend.title = 'Regions'
    p.legend.title_text_font_style = 'bold'
    p.legend.border_line_color = None
    p.legend.background_fill_color = 'white'
     
    label_opt = dict(x_units='data', y_units='data', text_font_size='10pt', render_mode='canvas')

    stats_ditau_fraction = str(round(100*float(c_met_trg['ditau']+c_tau_trg['ditau'])/(c_ditau_trg['ditau']+c_met_trg['ditau']+c_tau_trg['ditau']),2))
    stats_ditau = {'ditau': Label(x=6., y=9.3, text=ditau+': '+str(c_ditau_trg['ditau']),
                                  text_color='green', **label_opt),
                   'met':   Label(x=6., y=8.9, text='met && !'+ditau+': '+str(c_met_trg['ditau']),
                                  text_color='blue', **label_opt),
                   'tau':   Label(x=6., y=8.5, text=tau+' && !met && !'+ditau+': '+str(c_tau_trg['ditau']),
                                  text_color='red', **label_opt),
                   'gain':  Label(x=6., y=8.1, text='Gain: '+stats_ditau_fraction+'%',
                                  text_color='black', **label_opt),}
    for elem in stats_ditau.values():
        p.add_layout(elem)

    try:
        contam_tau = str(round(100*float(c_tau_trg['met'])/(c_met_trg['met']+c_tau_trg['met']),2))
    except ZeroDivisionError:
        contam_tau = '0'
    stats_met = {'met':   Label(x=2.6, y=1.3, text='met: '+str(c_met_trg['met']),
                                text_color='blue', **label_opt),
                 'tau':   Label(x=2.6, y=0.9, text=tau+' && !met: '+str(c_tau_trg['met']),
                                text_color='red', **label_opt),
                 'ditau': Label(x=2.6, y=0.5, text=ditau+': '+str(c_ditau_trg['met']),
                                text_color='green', **label_opt),
                 'contamination': Label(x=2.6, y=0.1, text='Contam.: '+contam_tau+'%',
                                        text_color='black', **label_opt),}
    for elem in stats_met.values():
        p.add_layout(elem)

    try:
        contam_met = str(round(100*float(c_met_trg['tau'])/(c_met_trg['tau']+c_tau_trg['tau']),2))
    except ZeroDivisionError:
        contam_met = '0'
        
    stats_tau = {'tau':   Label(x=6.5, y=1.3, text=tau+': '+str(c_tau_trg['tau']),
                                 text_color='red', **label_opt),
                 'met':   Label(x=6.5, y=0.9, text='met && !'+tau+': '+str(c_met_trg['tau']),
                                text_color='blue', **label_opt),
                 'ditau': Label(x=6.5, y=0.5, text=ditau+': '+str(c_ditau_trg['tau']),
                                text_color='green', **label_opt),
                 'contamination': Label(x=6.5, y=0.1, text='Contam.: '+contam_met+'%',
                                        text_color='black', **label_opt),}

    for elem in stats_tau.values():
        p.add_layout(elem)
     
    line_opt = dict(color='black', line_dash='dashed', line_width=2)
    p.line(x=[2.0, 2.0], y=[0.0, 10.0], **line_opt)
    p.line(x=[0.0, 10.0], y=[2.0, 2.0], **line_opt)
    if pt_cuts[0] != '40':
        p.line(x=[2.5, 2.5], y=[0.0, 10.0], **line_opt)
    if pt_cuts[1] != '40':
        p.line(x=[0.0, 10.0], y=[2.5, 2.5], **line_opt)
     
    p.xaxis.ticker = [2.0, 2.5, 5.5]
    p.xaxis.major_label_overrides = {2: '40', 2.5: pt_cuts[0], 5.5: region_cuts[0]}
    p.xaxis.axis_label = 'dau1_pT [GeV]'
     
    p.yaxis.ticker = [2.0, 2.5, 5.5]
    p.yaxis.major_label_overrides = {2: '40', 2.5: pt_cuts[1], 5.5: region_cuts[1]}
    p.yaxis.axis_label = 'dau2_pT [GeV]'
     
    p.output_backend = "svg"
    save(p)
    
def get_outname(suffix, mode, cut='', ext=''):
    pref = {'met': 'MET', 'tau': 'Tau', 'met_tau': 'MET_and_Tau'}
    assert mode in list(pref.keys())
    if ext == 'csv':
        s = 'counts_{}_{}'.format(suffix, pref[mode])
        if cut:
            s += '_{}'.format(cut)
        s += '/table.csv'
    elif ext == 'root':
        utils.create_single_dir('data')
        s = 'data/regions_{}_{}'.format(suffix, pref[mode])
        if cut:
            s += '_{}'
        s += '.{}'.format(ext)
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

def plot(hbase, hmet, htau, var, channel, sample, region, category, directory):
    defs = set_plot_definitions()
    
    hbase2 = hbase.Clone('hbase2')
    hmet2  = hmet.Clone('hmet2')
    htau2  = htau.Clone('htau2')

    # Normalized shapes
    c2 = ROOT.TCanvas('c2', '', 600, 400)
    c2.cd()
    try:
        hbase2.Scale(1/hbase2.Integral())
    except ZeroDivisionError:
        pass
    try:
        hmet2.Scale(1/hmet2.Integral())
    except ZeroDivisionError:
        pass
    try:
        htau.Scale(1/htau.Integral())
    except ZeroDivisionError:
        pass

    shift_scale = 5.
    amax = hbase2.GetMaximum() + (hbase2.GetMaximum()-hbase2.GetMinimum()) /  shift_scale
    amax = max(amax, hmet2.GetMaximum() + (hmet2.GetMaximum()-hmet2.GetMinimum()) /  shift_scale)
    amax = max(amax, htau2.GetMaximum() + (htau2.GetMaximum()-htau2.GetMinimum()) /  shift_scale)
    hbase2.SetMaximum(amax)
    
    hbase2.GetXaxis().SetTitleSize(defs['XTitleSize']);
    hbase2.GetXaxis().SetTitle(var + ' [GeV]');
    hbase2.GetYaxis().SetTitleSize(defs['YTitleSize']);
    hbase2.GetYaxis().SetTitle('Normalized to 1');
    hbase2.SetLineWidth(defs['LineWidth']);
    hbase2.SetLineColor(4);

    hmet2.SetLineWidth(defs['LineWidth']);
    hmet2.SetLineColor(2);
    hmet2.SetLineWidth(defs['LineWidth']);
    hmet2.SetLineColor(8);

    htau2.SetLineWidth(defs['LineWidth']);
    htau2.SetLineColor(2);
    htau2.SetLineWidth(defs['LineWidth']);
    htau2.SetLineColor(8);

    hbase2.Draw('hist')
    hmet.Draw('histsame')
    htau.Draw('histsame')

    leg2 = ROOT.TLegend(0.69, 0.77, 0.90, 0.9)
    leg2.SetNColumns(1)
    leg2.SetFillStyle(0)
    leg2.SetBorderSize(0)
    leg2.SetTextFont(43)
    leg2.SetTextSize(10)
    leg2.AddEntry(hmet2, '!ditau && MET')
    leg2.AddEntry(htau2, '!ditau && !MET && tau')
    leg2.AddEntry(hbase2, '+\n'.join(triggers[channel]))
    leg2.Draw('same')
    
    c2.Update();
    cat_dir = os.path.join(directory, sample, category)
    utils.create_single_dir(cat_dir)
    for ext in ('png', 'pdf'):
        c2.SaveAs( os.path.join(cat_dir, 'norm_' + reg + '_' + var + '.' + ext) )
    c2.Close()

def plot2D(histo, two_vars, channel, sample, suffix, category, directory, region):
    histo = histo.Clone('histo_clone')
    defs = set_plot_definitions()    

    c = ROOT.TCanvas('c', '', 800, 800)
    c.cd()
    
    pad1 = ROOT.TPad('pad1', 'pad1', 0., 0., 1., 1.)
    pad1.SetFrameLineWidth(defs['FrameLineWidth'])
    pad1.SetLeftMargin(0.12);
    pad1.SetRightMargin(0.1);
    pad1.SetBottomMargin(0.12);
    pad1.SetTopMargin(0.04);
    pad1.Draw()
    pad1.cd()

    histo.GetXaxis().SetTitle(two_vars[0])
    histo.GetXaxis().SetTitleSize(0.04)
    histo.GetYaxis().SetTitle(two_vars[1])
    histo.GetYaxis().SetTitleSize(0.04)
    histo.Draw('colz')

    c.cd()
    cat_folder = os.path.join(directory, sample, category)
    utils.create_single_dir(cat_folder)
    c.Update();
    for ext in ('png', 'pdf'):
        c.SaveAs( os.path.join(cat_folder, suffix + '_' + 'reg' + region + '_' + '_VS_'.join(two_vars) + '.' + ext) )
    c.Close()

def test_trigger_regions(indir, sample, channel):
    outname = get_outname(suffix=sample+'_'+channel, mode='met', ext='root')

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

    jBase, jNoBase_MET, jNoBase_NoMET_Tau = (defaultdict(lambda: defaultdict(dict)) for _ in range(3)) # one trigger, 1-dimensional
    hBase, hMET, hTau = (defaultdict(lambda: defaultdict(dict)) for _ in range(3)) # one trigger
    hBase_MET, hBase_Tau, hMET_Tau = (defaultdict(lambda: defaultdict(dict)) for _ in range(3)) # intersection of two triggers
    hNoBase_MET, hNoBase_NoMET_Tau = (defaultdict(lambda: defaultdict(dict)) for _ in range(2)) # negations
    hNoBase_Tau, hNoBase_MET_NoTau = (defaultdict(lambda: defaultdict(dict)) for _ in range(2)) # negations
    hBase_MET_Tau = defaultdict(lambda: defaultdict(dict)) # intersection of three triggers
    norphans, ntotal = ({k:0 for k in categories} for _ in range(2))
    for reg in regions:
        for v in tuple(variables):
            for cat in categories:
                suff = lambda x : x + '_' + reg + '_' + v + '_' + cat
                jopt = ('', *binning[v])
                jBase[reg][v][cat] = ROOT.TH1D(suff('jBase'), *jopt)
                jNoBase_MET[reg][v][cat] = ROOT.TH1D(suff('jNoBase_MET'), *jopt)
                jNoBase_NoMET_Tau[reg][v][cat] = ROOT.TH1D(suff('jNoBase_NoMET_Tau'), *jopt)
                
        for v in tuple(variables_2D):
            for cat in categories:
                suff = lambda x: x + '_' + reg + '_' + v[0] + '_VS_' + v[1] + '_' + cat
                hopt = ('', *binning[v[0]], *binning[v[1]])

                hBase[reg][v][cat] = ROOT.TH2D(suff('hBase'), *hopt)
                hMET[reg][v][cat]  = ROOT.TH2D(suff('hMET'), *hopt)
                hTau[reg][v][cat]  = ROOT.TH2D(suff('hTau'), *hopt)
                
                hBase_MET[reg][v][cat] = ROOT.TH2D(suff('hBase_MET'), *hopt)
                hBase_Tau[reg][v][cat] = ROOT.TH2D(suff('hBase_Tau'), *hopt)
                hMET_Tau[reg][v][cat]  = ROOT.TH2D(suff('hMET_Tau') , *hopt)

                hNoBase_MET[reg][v][cat]       = ROOT.TH2D(suff('hNoBase_MET'), *hopt)
                hNoBase_NoMET_Tau[reg][v][cat] = ROOT.TH2D(suff('hNoBase_NoMET_Tau'), *hopt)
                hNoBase_Tau[reg][v][cat]       = ROOT.TH2D(suff('hNoBase_Tau'), *hopt)
                hNoBase_MET_NoTau[reg][v][cat] = ROOT.TH2D(suff('hNoBase_MET_NoTau'), *hopt)

                hBase_MET_Tau[reg][v][cat] = ROOT.TH2D(suff('hBase_MET_Tau'), *hopt)

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
                    if sel.sel_category(cat) and entries.ditau_deltaR > 0.5:

                        if v == variables_2D[0]:
                            ntotal[cat] += 1

                        if eval(met_region):
                            reg = 'met'
                        elif eval(tau_region):
                            reg = 'tau'
                        elif eval(ditau_region):
                            reg = 'ditau'

                        met_turnon_expr = entries.metnomu_et > met_turnon
                        tau_turnon_expr = ((entries.dau1_pt > tau_turnon and args.channel=='tautau') or
                                           (entries.dau2_pt > tau_turnon and args.channel!='tautau'))
                        pass_trg = sel.pass_triggers(triggers[channel]) and entries.isLeptrigger
                        pass_met = sel.pass_triggers(('METNoMu120',)) and met_turnon_expr and not entries.isLeptrigger
                        pass_tau = sel.pass_triggers(('IsoTau180',)) and tau_turnon_expr and not entries.isLeptrigger
                        
                        if pass_trg:
                            jBase[reg][v][cat].Fill(entries[v], evt_weight)

                        if not pass_trg and pass_met:
                            jNoBase_MET[reg][v][cat].Fill(entries[v], evt_weight)

                        if not pass_trg and not pass_met and pass_tau:
                            jNoBase_NoMET_Tau[reg][v][cat].Fill(entries[v], evt_weight)

            for v in variables_2D:
                vx, vy = v
                for cat in categories:
                    if sel.sel_category(cat) and entries.ditau_deltaR > 0.5:

                        if v == variables_2D[0]:
                            ntotal[cat] += 1

                        if eval(met_region):
                            reg = 'met'
                        elif eval(tau_region):
                            reg = 'tau'
                        elif eval(ditau_region):
                            reg = 'ditau'
                        else:
                            if v == variables_2D[0]:
                                norphans[cat] += 1
                            continue

                        met_turnon_expr = entries.metnomu_et > met_turnon
                        tau_turnon_expr = ((entries.dau1_pt > tau_turnon and args.channel=='tautau') or
                                           (entries.dau2_pt > tau_turnon and args.channel!='tautau'))
                        pass_trg = sel.pass_triggers(triggers[channel]) and entries.isLeptrigger
                        pass_met = sel.pass_triggers(('METNoMu120',)) and met_turnon_expr and not entries.isLeptrigger
                        pass_tau = sel.pass_triggers(('IsoTau180',)) and tau_turnon_expr and not entries.isLeptrigger
                        
                        if pass_trg:
                            hBase[reg][v][cat].Fill(entries[vx], entries[vy], evt_weight)

                        if pass_met:
                            hMET[reg][v][cat].Fill(entries[vx], entries[vy], evt_weight)
                            
                        if pass_tau:
                            hTau[reg][v][cat].Fill(entries[vx], entries[vy], evt_weight)

                        # intersections with two triggers
                        if pass_trg and pass_met:
                            hBase_MET[reg][v][cat].Fill(entries[vx], entries[vy], evt_weight)

                        if pass_trg and pass_tau:
                            hBase_Tau[reg][v][cat].Fill(entries[vx], entries[vy], evt_weight)

                        if pass_met and pass_tau:
                            hMET_Tau[reg][v][cat].Fill(entries[vx], entries[vy], evt_weight)

                        # negations
                        if not pass_trg and pass_met:
                            hNoBase_MET[reg][v][cat].Fill(entries[vx], entries[vy], evt_weight)

                        if not pass_trg and not pass_met and pass_tau:
                            hNoBase_NoMET_Tau[reg][v][cat].Fill(entries[vx], entries[vy], evt_weight)

                        if not pass_trg and pass_tau:
                            hNoBase_Tau[reg][v][cat].Fill(entries[vx], entries[vy], evt_weight)

                        if not pass_trg and pass_met and not pass_tau:
                            hNoBase_MET_NoTau[reg][v][cat].Fill(entries[vx], entries[vy], evt_weight)

                        # intersection with three triggers
                        if pass_trg and pass_met and pass_tau:
                            hBase_MET_Tau[reg][v][cat].Fill(entries[vx], entries[vy], evt_weight)

    f_out = ROOT.TFile(outname, 'RECREATE')
    f_out.cd()
    for reg in regions:
        for cat in categories:
            for v in variables:
                suff = lambda x : x + '_' + reg + '_' + v + '_' + cat
                jBase[reg][v][cat].Write(suff('jBase'))
                jNoBase_MET[reg][v][cat].Write(suff('jNoBase_MET'))
                jNoBase_NoMET_Tau[reg][v][cat].Write(suff('jNoBase_NoMET_Tau'))
                
            for v in variables_2D:
                suff = lambda x: x + '_' + reg + '_' + v[0] + '_VS_' + v[1] + '_' + cat
                hBase[reg][v][cat].Write(suff('hBase'))
                hMET[reg][v][cat].Write(suff('hMET'))
                hTau[reg][v][cat].Write(suff('hTau'))

                hBase_MET[reg][v][cat].Write(suff('hBase_MET'))
                hBase_Tau[reg][v][cat].Write(suff('hBase_Tau'))
                hMET_Tau[reg][v][cat].Write(suff('hMET_Tau'))

                hNoBase_MET[reg][v][cat].Write(suff('hNoBase_MET'))
                hNoBase_NoMET_Tau[reg][v][cat].Write(suff('hNoBase_NoMET_Tau'))
                hNoBase_Tau[reg][v][cat].Write(suff('hNoBase_Tau'))
                hNoBase_MET_NoTau[reg][v][cat].Write(suff('hNoBase_MET_NoTau'))

                hBase_MET_Tau[reg][v][cat].Write(suff('hBase_MET_Tau'))
                
    f_out.Close()
    for cat in categories:
        print('Category {} (m(X)={}GeV) had {} orphans ({}%)'.format(cat, sample, norphans[cat], float(norphans[cat])/ntotal[cat]))
    print('Raw histograms saved in {}.'.format(outname), flush=True)

if __name__ == '__main__':
    triggers = {#'etau': ('Ele32', 'EleIsoTauCustom'),
                #'mutau': ('IsoMu24', 'IsoMuIsoTauCustom'),
                'tautau': ('IsoDoubleTauCustom',)
                }
    binning = {'metnomu_et': (20, 0, 400),
               'dau1_pt': (30, 0, 350),
               'dau1_eta': (20, -2.5, 2.5),
               'dau2_iso': (20, 0.88, 1.01),
               'dau2_pt': (30, 0, 300),
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
    variables_2D = (('dau1_pt',  'dau2_pt'),
                    ('dau1_iso', 'dau2_iso'),
                    # ('dau1_eta', 'dau2_eta'),
                    # ('dau1_eta', 'dau1_iso'),
                    # ('dau1_pt',  'dau1_eta'),
                    # ('dau1_pt',  'dau1_iso'),
                    # ('dau2_eta', 'dau2_iso'),
                    # ('dau2_pt',  'dau2_eta'),
                    # ('dau1_pt',  'dau2_eta'),
                    # ('dau1_pt',  'dau2_iso'),
                    # ('dau2_pt',  'dau1_eta'),
                    # ('dau2_pt',  'dau1_iso'),
                    ('dau1_pt',  'metnomu_et'),
                    ('dau2_pt',  'metnomu_et'),
                    ('dau1_iso',  'metnomu_et'),
                    ('dau2_iso',  'metnomu_et'),
                    )

    #categories = ('baseline', 's1b1jresolvedMcut', 's2b0jresolvedMcut', 'sboostedLLMcut')
    categories = ('baseline',)
    met_turnon = 200
    tau_turnon = 190
    
    # Parse input arguments
    parser = argparse.ArgumentParser(description='Producer trigger histograms.')

    parser.add_argument('--indir', required=True, type=str,
                        help='Full path of ROOT input file')
    parser.add_argument('--samples', required=True, nargs='+', type=str,
                        help='Full path of ROOT input file')
    parser.add_argument('--channel', required=True, type=str,  
                        help='Select the channel over which the workflow will be run.' )
    parser.add_argument('--plot', action='store_true',
                        help='Reuse previously produced data for quick plot changes.')
    parser.add_argument('--copy', action='store_true',
                        help='Do not copy the outputs to EOS at the end.')
    parser.add_argument('--sequential', action='store_true',
                        help='Do not use the multiprocess package.')
    args = utils.parse_args(parser)

    region_cuts = '200', '200'
    pt_cuts = '40', '40'
    main_dir = 'Region_CUT1_' + region_cuts[0] + '_' + region_cuts[1] + '_CUT2_' + pt_cuts[0] + '_' + pt_cuts[1]

    regions = ('ditau', 'met', 'tau')
    met_region = ('(entries.dau2_pt < 40 and entries.dau1_pt < {}) or '.format(region_cuts[0]) +
                  '(entries.dau1_pt < 40 and entries.dau2_pt < {})'.format(region_cuts[1])) 
    tau_region = 'entries.dau2_pt < 40 and entries.dau1_pt >= {}'.format(region_cuts[0])
    ditau_region = 'entries.dau1_pt > {} and entries.dau2_pt > {}'.format(pt_cuts[0], pt_cuts[1])

    #### run major function ###
    if args.sequential:
        if not args.plot:
            for sample in args.samples:
                test_trigger_regions(args.indir, sample, args.channel)
    else:
        if not args.plot:
            pool = multiprocessing.Pool(processes=6)
            pool.starmap(test_trigger_regions,
                         zip(it.repeat(args.indir), args.samples, it.repeat(args.channel)))
    ###########################

    from_directory = os.path.join(main_dir, args.channel)
    for sample in args.samples:
        outname = get_outname(suffix=sample+'_'+args.channel, mode='met', ext='root')
        f_in = ROOT.TFile(outname, 'READ')
        f_in.cd()

        # write csv header, one per category
        out_counts = []
        for cat in categories:
            out_counts.append( os.path.join(from_directory, sample, cat, 'counts') )
            utils.create_single_dir(out_counts[-1])
            with open(os.path.join(out_counts[-1], 'table.csv'), 'w') as f:
                reader = csv.writer(f, delimiter=',', quotechar='|')
                header_row = ['Region', 'Base', 'MET', 'Tau',
                              'Base && MET', 'Base && Tau', 'MET && Tau',
                              '!Base && MET', '!Base && !MET && Tau',
                              '!Base && Tau', '!Base && MET && !Tau',
                              'Base && MET && Tau']
                reader.writerow(header_row)

        # plot histograms and fill CSV with histogram integrals
        c_ditau_trg, c_met_trg, c_tau_trg = ({} for _ in range(3))
        
        for reg in regions:
            for cat in categories:
                for v in variables:                    
                    suff = lambda x : x + '_' + reg + '_' + v + '_' + cat
                    jBase = f_in.Get(suff('jBase'))
                    jNoBase_MET = f_in.Get(suff('jNoBase_MET'))
                    jNoBase_NoMET_Tau = f_in.Get(suff('jNoBase_NoMET_Tau'))

                    jopt = (v, args.channel, sample, reg, cat, from_directory)
                    plot(jBase, jNoBase_MET, jNoBase_NoMET_Tau, *jopt)
                    
                for v in variables_2D:
                    suff = lambda x : x+'_'+reg+'_'+v[0]+'_VS_'+v[1]+'_'+cat
                    
                    hBase = f_in.Get(suff('hBase'))
                    hMET  = f_in.Get(suff('hMET'))
                    hTau  = f_in.Get(suff('hTau'))
 
                    hBase_MET = f_in.Get(suff('hBase_MET'))
                    hBase_Tau = f_in.Get(suff('hBase_Tau'))
                    hMET_Tau  = f_in.Get(suff('hMET_Tau'))

                    hNoBase_MET       = f_in.Get(suff('hNoBase_MET'))
                    hNoBase_NoMET_Tau = f_in.Get(suff('hNoBase_NoMET_Tau'))
                    hNoBase_Tau       = f_in.Get(suff('hNoBase_Tau'))
                    hNoBase_MET_NoTau = f_in.Get(suff('hNoBase_MET_NoTau'))
 
                    hBase_MET_Tau = f_in.Get(suff('hBase_MET_Tau'))

                    plot2D(hBase, v, args.channel, sample, 'Base', cat, from_directory, reg)
                    plot2D(hMET, v, args.channel, sample, 'MET', cat, from_directory, reg)
                    plot2D(hTau, v, args.channel, sample, 'Tau', cat, from_directory, reg)

                # count only for one of the variables, it does not matter which
                nbinsx, nbinsy = hBase.GetNbinsX(), hBase.GetNbinsY()
                cBase = round(hBase.Integral(0, nbinsx+1, 0, nbinsy+1), 2)
                cMET  = round(hMET.Integral(0, nbinsx+1, 0, nbinsy+1), 2)
                cTau  = round(hTau.Integral(0, nbinsx+1, 0, nbinsy+1), 2)
                
                cBase_MET = round(hBase_MET.Integral(0, nbinsx+1, 0, nbinsy+1), 2)
                cBase_Tau = round(hBase_Tau.Integral(0, nbinsx+1, 0, nbinsy+1), 2)
                cMET_Tau  = round(hMET_Tau.Integral(0, nbinsx+1, 0, nbinsy+1), 2)
                
                cNoBase_MET       = round(hNoBase_MET.Integral(0, nbinsx+1, 0, nbinsy+1), 2)
                cNoBase_NoMET_Tau = round(hNoBase_NoMET_Tau.Integral(0, nbinsx+1, 0, nbinsy+1), 2)
                cNoBase_Tau       = round(hNoBase_Tau.Integral(0, nbinsx+1, 0, nbinsy+1), 2)
                cNoBase_MET_NoTau = round(hNoBase_MET_NoTau.Integral(0, nbinsx+1, 0, nbinsy+1), 2)
                
                cBase_MET_Tau = round(hBase_MET_Tau.Integral(0, nbinsx+1, 0, nbinsy+1), 2)

                if reg=='ditau':
                    subsets = (cBase, cMET, cBase_MET, cTau, cBase_Tau, cMET_Tau, cBase_MET_Tau)
                    plt.figure(figsize=(10,8))
                    v = venn3(subsets=subsets, set_labels = (r'$\tau\tau$', 'MET', 'Tau'))
                    v.get_patch_by_id('100').set_alpha(1.0)
                    c = venn3_circles(subsets=subsets, linestyle='dashed')
                    plt.title('Baseline selection')
                    for ext in ('png', 'pdf'):
                        plt.savefig(os.path.join(out_counts[categories.index(cat)], 'venn.' + ext))

                with open(os.path.join(out_counts[categories.index(cat)], 'table.csv'), 'a') as f:
                    reader = csv.writer(f, delimiter=',', quotechar='|')
                    row = [reg, cBase, cMET, cTau,
                           cBase_MET, cBase_Tau, cMET_Tau,
                           cNoBase_MET, cNoBase_NoMET_Tau, cNoBase_Tau, cNoBase_MET_NoTau,
                           cBase_MET_Tau]
                    reader.writerow(row)

                if reg=='ditau':
                    c_ditau_trg[reg] = cBase
                    c_met_trg[reg] = cNoBase_MET
                    c_tau_trg[reg] = cNoBase_NoMET_Tau
                elif reg=='met':
                    c_ditau_trg[reg] = 0.
                    c_met_trg[reg] = cMET
                    c_tau_trg[reg] = cNoBase_NoMET_Tau
                elif reg=='tau':
                    c_ditau_trg[reg] = 0.
                    c_met_trg[reg] = cNoBase_MET_NoTau
                    c_tau_trg[reg] = cTau
                        
        f_in.Close()

        square_diagram(sample, c_ditau_trg, c_met_trg, c_tau_trg,
                       region_cuts, pt_cuts,
                       output_name=os.path.join(out_counts[categories.index(cat)], 'diagram.html'))
    
    if args.copy:
        import subprocess
        to_directory = os.path.join('/eos/user/b/bfontana/www/TriggerScaleFactors', main_dir)
        to_directory = os.path.join(to_directory, args.channel)
     
        for sample in args.samples:
            sample_from = os.path.join(from_directory, sample)
            print('Copying: {}\t\t--->\t{}'.format(sample_from, to_directory), flush=True)
            subprocess.run(['rsync', '-ah', sample_from, to_directory])
        
    print('Done.')
