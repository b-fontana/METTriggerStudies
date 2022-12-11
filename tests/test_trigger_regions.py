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
import numpy as np
import h5py
from collections import defaultdict

import inclusion
from inclusion import selection
from inclusion.config import main
from inclusion.utils import utils

import ROOT

from bokeh.plotting import figure, output_file, save
from bokeh.models import Range1d, Label
import matplotlib.pyplot as plt
#from matplotlib_venn import venn3, venn3_circles

tau = '\u03C4'
mu  = '\u03BC'
pm  = '\u00B1'
ditau = tau+tau

def contamination_save(savepath, label, c1, c2, e1, e2, m):
    with h5py.File(os.path.join(savepath, label + '.hdf5'), 'w') as f:
        dset = f.create_dataset(label, (5,len(m)), dtype='f')
        dset[0, :] = m
        dset[1, :] = c1
        dset[2, :] = c2
        dset[3, :] = e1
        dset[4, :] = e2
        dset.cols = ['mass [GeV]',
                     'contamination ' + ditau + '(%)',
                     'contamination ' + ditau + ' MET (%)',
                     'uncertainty' + ditau,
                     'uncertainty' + ditau + ' MET']
        
def square_diagram(c_ditau_trg, c_met_trg, c_tau_trg, channel,
                   region_cuts, pt_cuts, text, bigtau=False, notau=False, nomet=False):
    base = {'etau': 'e+e'+tau, 'mutau': mu+'+'+mu+tau, 'tautau': ditau}
    output_file(text['out'])

    topr = 9.9
    shft = 0.1
    start, b1, b2, b3 = 0.0, 2, 2.5, 7.5
    xgap = b1+shft if pt_cuts[0] == '40' and pt_cuts[1] == '40' else b2+sft
    
    p = figure(title='m(X)={}GeV'.format(text['mass']), width=600, height=400,
               tools='save')
    p.x_range = Range1d(0, 11.5)
    p.y_range = Range1d(0, 10)
    p.outline_line_color = None
    p.toolbar.logo = None
    p.xgrid.grid_line_color = None
    p.ygrid.grid_line_color = None
    p.xaxis.ticker = [b1,b3] if pt_cuts[0] == '40' else [b1, b2, b3]
    p.xaxis.major_label_overrides = ({b1: '40',
                                      b3: region_cuts[0]} if pt_cuts[0] == '40' and pt_cuts[1] == '40'
                                     else {b1: '40',
                                           b2: pt_cuts[0],
                                           b3: region_cuts[0]})
    p.xaxis.axis_label = 'dau1_pT [GeV]'
     
    p.yaxis.ticker = [b1, b3] if pt_cuts[1] == '40' else [b1, b2, b3]
    p.yaxis.major_label_overrides = ({b1: '40',
                                      b3: region_cuts[1]} if pt_cuts[0] == '40' and pt_cuts[1] == '40'
                                     else {b1: '40',
                                           b2: pt_cuts[1], b3:
                                           region_cuts[1]})
    p.yaxis.axis_label = 'dau2_pT [GeV]'

    # add a square renderer with a size, color, and alpha
    polyg_opt = dict(alpha=0.3)
    p.multi_polygons(color='green',
                     xs=[[[[xgap,b3-shft,b3-shft,xgap]]]] if bigtau else [[[[xgap,topr,topr,xgap]]]],
                     ys=[[[[b3-shft,b3-shft,xgap,xgap]]]] if bigtau else [[[[topr,topr,xgap,xgap]]]],
                     legend_label=base[channel], **polyg_opt)
    if not nomet:
        p.multi_polygons(color='blue',
                         xs=[[[[start+shft,b1-shft,b1-shft,b3-shft,b3-shft,shft]]]],
                         ys=[[[[b3-shft,b3-shft,b1-shft,b1-shft,shft,shft]]]],
                         legend_label='MET', **polyg_opt)
    if not notau:
        p.multi_polygons(color='red',
                         xs=([[[[shft,shft,topr,topr,b3+shft,b3+shft,shft]]]] if bigtau else
                             [[[[shft,shft,b1-shft,b1-shft], [b3+shft,b3+shft,topr,topr]]]]),
                         ys=([[[[b3+shft,topr,topr,shft,shft,b3+shft,b3+shft]]]] if bigtau else
                             [[[[b3+shft,topr,topr,b3+shft], [shft,b1-shft,b1-shft,shft]]]]),
                         legend_label=tau, **polyg_opt)
    p.legend.title = 'Regions'
    p.legend.title_text_font_style = 'bold'
    p.legend.border_line_color = None
    p.legend.background_fill_color = 'white'
    p.legend.click_policy = 'hide'
     
    label_opt = dict(x_units='data', y_units='data', text_font_size='10pt')

    gain = (100*float(c_met_trg['ditau']+c_tau_trg['ditau']) /
            (c_ditau_trg['ditau']+c_met_trg['ditau']+c_tau_trg['ditau']))
    gain = str(round(gain,2))
    stats_ditau = {'ditau': Label(x=b1+0.3, y=b1+1.5, text=ditau+': '+str(c_ditau_trg['ditau']),
                                  text_color='black', **label_opt),
                   'met':   Label(x=b1+0.3, y=b1+1.1, text='met && !'+ditau+': '+str(c_met_trg['ditau']),
                                  text_color='black', **label_opt),
                   'tau':   Label(x=b1+0.3, y=b1+0.7, text=tau+' && !met && !'+ditau+': '+str(c_tau_trg['ditau']),
                                  text_color='black', **label_opt),
                   'gain':  Label(x=b1+0.3, y=b1+0.3, text='Gain: '+gain+'%',
                                  text_color='blue', **label_opt),}
    for key,elem in stats_ditau.items():
        if nomet and notau and key=='gain':
            continue
        p.add_layout(elem)

    if not nomet:
        try:
            contam_by_tau = (100*(float(c_tau_trg['met'])+c_ditau_trg['met']) /
                             (c_met_trg['met']+c_tau_trg['met']+c_ditau_trg['met']))
            contam_by_tau = str(round(contam_by_tau,2))
        except ZeroDivisionError:
            contam_by_tau = '0'
        stats_met = {'met':   Label(x=b1+0.3, y=1.3, text='met: '+str(c_met_trg['met']),
                                    text_color='black', **label_opt),
                     'tau':   Label(x=b1+0.3, y=0.5, text=tau+' && !'+ditau+' && !met: '+str(c_tau_trg['met']),
                                    text_color='black', **label_opt),
                     'ditau': Label(x=b1+0.3, y=0.9, text=ditau+' && !met: '+str(c_ditau_trg['met']),
                                    text_color='black', **label_opt),
                     'contamination': Label(x=b1+0.3, y=0.1, text='Contam.: '+contam_by_tau+'%',
                                            text_color='blue', **label_opt),}
        for elem in stats_met.values():
            p.add_layout(elem)

    if not notau:
        num_ditau = float(c_ditau_trg['tau'])
        num_both = num_ditau + float(c_met_trg['tau'])
        den = float(c_met_trg['tau']+c_tau_trg['tau']+c_ditau_trg['tau'])
        if num_ditau == 0 or num_both == 0:
            contam_ditau = '0'
            contam_both = '0'
            err_ditau = '0'
            err_both = '0'
        else:
            contam_ditau = 100*num_ditau/den
            contam_both = 100*num_both/den

            enum_ditau = np.sqrt(c_ditau_trg['tau'])
            eden_ditau = enum_ditau + np.sqrt(c_met_trg['tau']) + np.sqrt(c_tau_trg['tau'])
            enum_both  = np.sqrt(c_met_trg['tau']) + np.sqrt(c_ditau_trg['tau'])
            eden_both  = enum_both + np.sqrt(c_tau_trg['tau'])
            err_ditau = contam_ditau * np.sqrt(enum_ditau**2/num_ditau**2 + eden_ditau**2/den**2)
            err_both  = contam_both * np.sqrt(enum_both**2/num_both**2 + eden_both**2/den**2)
        
            contam_ditau = str(round(contam_ditau,2))
            contam_both  = str(round(contam_both,2))
            err_ditau    = str(round(err_ditau,2))
            err_both     = str(round(err_both,2))

        stats_tau = {'tau':   Label(x=b3+0.2, y=1.3, text=tau+': '+str(c_tau_trg['tau']),
                                    text_color='black', **label_opt),
                     'met':   Label(x=b3+0.2, y=0.5, text='met && !'+ditau+' && !'+tau+': '+str(c_met_trg['tau']),
                                    text_color='black', **label_opt),
                     'ditau': Label(x=b3+0.2, y=0.9, text=ditau+' && !'+tau+': '+str(c_ditau_trg['tau']),
                                    text_color='black', **label_opt),
                     'contamination_both': Label(x=b3+0.2, y=0.1,
                                                  text='Contam.: ('+str(contam_both)+pm+str(err_both)+')%',
                                                  text_color='blue', **label_opt),}
        for elem in stats_tau.values():
            p.add_layout(elem)
     
    line_opt = dict(color='black', line_dash='dashed', line_width=2)
    p.line(x=[b1,b1], y=[start, topr+shft], **line_opt)
    p.line(x=[start,topr+shft], y=[b1,b1], **line_opt)
    if pt_cuts[0] != '40':
        p.line(x=[b2,b2], y=[start,topr+shft], **line_opt)
    if pt_cuts[1] != '40':
        p.line(x=[start,topr+shft], y=[b2,b2], **line_opt)
    p.line(x=[b3,b3], y=[start,topr+shft], **line_opt)
    p.line(x=[start,topr+shft], y=[b3,b3], **line_opt)
     
    p.output_backend = 'svg'
    save(p)
    return contam_ditau, contam_both, err_ditau, err_both
    
def get_outname(sample, channel, region_cuts, pt_cuts, met_turnon, tau_turnon,
                bigtau, notau, nomet):
    utils.create_single_dir('data')

    name = sample + '_' + channel + '_'
    name += '_'.join((*region_cuts, 'ptcuts', *pt_cuts, 'turnon', met_turnon, tau_turnon))
    if bigtau:
        name += '_BIGTAU'
    if notau:
        name += '_NOTAU'
    if nomet:
        name += '_NOMET'
    name += '.root'

    s = 'data/regions_{}'.format(name)
    return s

def set_plot_definitions():
    ROOT.gROOT.SetBatch(ROOT.kTRUE)
    ROOT.gStyle.SetOptStat(ROOT.kFALSE)
    ret = {'XTitleSize'     : 0.045,
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
        hbase2.Scale(0.)
    try:
        hmet2.Scale(1/hmet2.Integral())
    except ZeroDivisionError:
        hmet2.Scale(0.)
    try:
        htau2.Scale(1/htau2.Integral())
    except ZeroDivisionError:
        htau2.Scale(0.)

    shift_scale = 5.
    amax = hbase2.GetMaximum() + (hbase2.GetMaximum()-hbase2.GetMinimum()) /  shift_scale
    amax = max(amax, hmet2.GetMaximum() + (hmet2.GetMaximum()-hmet2.GetMinimum()) / shift_scale)
    amax = max(amax, htau2.GetMaximum() + (htau2.GetMaximum()-htau2.GetMinimum()) / shift_scale)
    hbase2.SetMaximum(amax)
    
    hbase2.GetXaxis().SetTitleSize(defs['XTitleSize']);
    hbase2.GetXaxis().SetTitle(var + ' [GeV]');
    hbase2.GetYaxis().SetTitleSize(defs['YTitleSize']);
    hbase2.GetYaxis().SetTitle('a. u.');
    hbase2.SetLineWidth(defs['LineWidth']);
    hbase2.SetLineColor(4);

    hmet2.SetLineWidth(defs['LineWidth']);
    hmet2.SetLineColor(8);

    htau2.SetLineWidth(defs['LineWidth']);
    htau2.SetLineColor(2);

    hbase2.Draw('hist')
    hmet2.Draw('histsame')
    htau2.Draw('histsame')

    leg2 = ROOT.TLegend(0.69, 0.77, 0.90, 0.9)
    leg2.SetNColumns(1)
    leg2.SetFillStyle(0)
    leg2.SetBorderSize(0)
    leg2.SetTextFont(43)
    leg2.SetTextSize(10)
    leg2.AddEntry(hbase2, 'ditau')
    leg2.AddEntry(hmet2, '!ditau && MET')
    leg2.AddEntry(htau2, '!ditau && !MET && tau')

    leg2.Draw('same')
    
    c2.Update();
    cat_dir = os.path.join(directory, sample, category)
    utils.create_single_dir(cat_dir)
    for ext in extensions:
        c2.SaveAs( os.path.join(cat_dir, 'norm_' + region + '_' + var + '.' + ext) )
    c2.Close()

def plot2D(histo, two_vars, channel, sample, trigger_str, category, directory, region, norm=False):
    histo = histo.Clone('histo_clone')
    defs = set_plot_definitions()    

    c = ROOT.TCanvas('c', '', 800, 800)
    c.cd()
    
    pad1 = ROOT.TPad('pad1', 'pad1', 0., 0., 1., 1.)
    pad1.SetFrameLineWidth(defs['FrameLineWidth'])
    pad1.SetLeftMargin(0.12);
    pad1.SetRightMargin(0.14);
    pad1.SetBottomMargin(0.12);
    pad1.SetTopMargin(0.1);
    pad1.Draw()
    pad1.cd()

    ROOT.gStyle.SetTitleX(0.44) #title X location
    ROOT.gStyle.SetTitleY(1.015) #title Y location
    ROOT.gStyle.SetTitleW(0.7) #title width
    ROOT.gStyle.SetTitleH(0.15) #title height
    ROOT.gStyle.SetTitleColor(ROOT.kBlue, 't')
    histo.SetTitle('Region: ' + region +' | Triggers: ' + trigger_str);
    histo.SetTitleSize(0.02, 't')
    histo.GetXaxis().SetTitle(two_vars[0])
    histo.GetXaxis().SetTitleSize(0.04)
    histo.GetYaxis().SetTitle(two_vars[1])
    histo.GetYaxis().SetTitleSize(0.04)

    if norm:
        try:
            histo.Scale(1/histo.Integral())
        except ZeroDivisionError:
            histo.Scale(0.)
    histo.Draw('colz')

    c.cd()
    cat_folder = os.path.join(directory, sample, category)
    utils.create_single_dir(cat_folder)
    c.Update();
    tstr = trigger_str.replace(' ', '_').replace('+', '_PLUS_').replace('!', 'NOT')
    for ext in extensions:
        c.SaveAs(os.path.join(cat_folder,
                              tstr + '_' + 'reg' + region + '_' + '_VS_'.join(two_vars) + '.' + ext))
    c.Close()

def test_trigger_regions(indir, sample, channel):
    outname = get_outname(sample, channel, region_cuts, pt_cuts[args.channel], met_turnon, tau_turnon,
                          args.bigtau, args.notau, args.nomet)

    if channel == 'etau' or channel == 'mutau':
        iso1 = (24, 0, 8)
    elif channel == 'tautau':
        iso1 = (20, 0.965, 1.005)
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
    hNoBase_MET, hBase_NoMET, hNoBase_NoMET_Tau = (defaultdict(lambda: defaultdict(dict)) for _ in range(3)) # negations
    hNoBase_Tau, hBase_NoTau, hNoBase_MET_NoTau = (defaultdict(lambda: defaultdict(dict)) for _ in range(3)) # negations
    hBase_MET_Tau = defaultdict(lambda: defaultdict(dict)) # intersection of three triggers
    hBase_NoMET_NoTau = defaultdict(lambda: defaultdict(dict)) # passes only the base
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
                hBase_NoMET[reg][v][cat]       = ROOT.TH2D(suff('hBase_NoMET'), *hopt)
                hNoBase_NoMET_Tau[reg][v][cat] = ROOT.TH2D(suff('hNoBase_NoMET_Tau'), *hopt)
                hNoBase_Tau[reg][v][cat]       = ROOT.TH2D(suff('hNoBase_Tau'), *hopt)
                hBase_NoTau[reg][v][cat]       = ROOT.TH2D(suff('hBase_NoTau'), *hopt)
                hNoBase_MET_NoTau[reg][v][cat] = ROOT.TH2D(suff('hNoBase_MET_NoTau'), *hopt)

                hBase_MET_Tau[reg][v][cat] = ROOT.TH2D(suff('hBase_MET_Tau'), *hopt)
                hBase_NoMET_NoTau[reg][v][cat] = ROOT.TH2D(suff('hBase_NoMET_NoTau'), *hopt)

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

        #triggers and turnon cuts
        met_turnon_expr = entries.metnomu_et > float(met_turnon)
        tau_turnon_expr = entries.dau1_pt > float(tau_turnon) or entries.dau2_pt > float(tau_turnon)
        pass_trg = sel.pass_triggers(triggers[channel]) and entries.isLeptrigger
        pass_met = sel.pass_triggers(('METNoMu120',)) and met_turnon_expr
        pass_tau = sel.pass_triggers(('IsoTau180',)) and tau_turnon_expr
        
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
                        else:
                            continue
                        assert reg in regions
                                                    
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
                        assert reg in regions
                                                
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

                        if pass_trg and not pass_met:
                            hBase_NoMET[reg][v][cat].Fill(entries[vx], entries[vy], evt_weight)

                        if not pass_trg and not pass_met and pass_tau:
                            hNoBase_NoMET_Tau[reg][v][cat].Fill(entries[vx], entries[vy], evt_weight)

                        if not pass_trg and pass_tau:
                            hNoBase_Tau[reg][v][cat].Fill(entries[vx], entries[vy], evt_weight)

                        if pass_trg and not pass_tau:
                            hBase_NoTau[reg][v][cat].Fill(entries[vx], entries[vy], evt_weight)

                        if not pass_trg and pass_met and not pass_tau:
                            hNoBase_MET_NoTau[reg][v][cat].Fill(entries[vx], entries[vy], evt_weight)

                        # intersection with three triggers
                        if pass_trg and pass_met and pass_tau:
                            hBase_MET_Tau[reg][v][cat].Fill(entries[vx], entries[vy], evt_weight)

                        # only passes the base
                        if pass_trg and not pass_met and not pass_tau:
                            hBase_NoMET_NoTau[reg][v][cat].Fill(entries[vx], entries[vy], evt_weight)

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
                hBase_NoMET[reg][v][cat].Write(suff('hBase_NoMET'))
                hNoBase_NoMET_Tau[reg][v][cat].Write(suff('hNoBase_NoMET_Tau'))
                hNoBase_Tau[reg][v][cat].Write(suff('hNoBase_Tau'))
                hBase_NoTau[reg][v][cat].Write(suff('hBase_NoTau'))
                hNoBase_MET_NoTau[reg][v][cat].Write(suff('hNoBase_MET_NoTau'))

                hBase_MET_Tau[reg][v][cat].Write(suff('hBase_MET_Tau'))
                hBase_NoMET_NoTau[reg][v][cat].Write(suff('hBase_NoMET_NoTau'))
                
    f_out.Close()
    for cat in categories:
        orph_frac = float(norphans[cat])/ntotal[cat]
        if orph_frac > 0.1:
            print('{}% orphans ({}/{}). This is unusual.'.format(orph_frac, norphans[cat], ntotal[cat]))
            print(met_region)
            print(tau_region)
            print(ditau_region)
        print('Category {} (m(X)={}GeV) had {} orphans ({}%)'.format(cat, sample, norphans[cat], orph_frac))
    print('Raw histograms saved in {}.'.format(outname), flush=True)

if __name__ == '__main__':
    extensions = ('png',) #('png', 'pdf')
    triggers = {'etau': ('Ele32', 'EleIsoTauCustom'),
                'mutau': ('IsoMu24', 'IsoMuIsoTauCustom'),
                'tautau': ('IsoDoubleTauCustom',)
                }
    binning = {'metnomu_et': (20, 0, 450),
               'dau1_pt': (30, 0, 450),
               'dau1_eta': (20, -2.5, 2.5),
               'dau2_iso': (20, 0.88, 1.005),
               'dau2_pt': (30, 0, 400),
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
                    # ('dau1_iso', 'dau2_iso'),
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
                    # ('dau1_pt',  'metnomu_et'),
                    # ('dau2_pt',  'metnomu_et'),
                    # ('dau1_iso',  'metnomu_et'),
                    # ('dau2_iso',  'metnomu_et'),
                    )

    categories = ('baseline',) #('baseline', 's1b1jresolvedMcut', 's2b0jresolvedMcut', 'sboostedLLMcut')
    
    # Parse input arguments
    parser = argparse.ArgumentParser(description='Producer trigger histograms.')

    parser.add_argument('--indir', required=True, type=str,
                        help='Full path of ROOT input file')
    parser.add_argument('--masses', required=True, nargs='+', type=str,
                        help='Resonance mass')
    parser.add_argument('--channel', required=True, type=str,  
                        help='Select the channel over which the workflow will be run.' )
    parser.add_argument('--plot', action='store_true',
                        help='Reuse previously produced data for quick plot changes.')
    parser.add_argument('--copy', action='store_true',
                        help='Copy the outputs to EOS at the end.')
    parser.add_argument('--sequential', action='store_true',
                        help='Do not use the multiprocess package.')
    parser.add_argument('--bigtau', action='store_true',
                        help='Consider a larger single tau region, reducing the ditau one.')
    parser.add_argument('--notau', action='store_true',
                        help='Remove the single tau region (default analysis).')
    parser.add_argument('--nomet', action='store_true',
                        help='Remove the MET region (default analysis).')
    parser.add_argument('--met_turnon', required=False, type=str,  default='200',
                        help='MET trigger turnon cut [GeV].' )
    args = utils.parse_args(parser)

    met_turnon = args.met_turnon
    tau_turnon = '190'
    region_cuts = ('190', '190')
    pt_cuts = {'etau':   ('20', '20'),
               'mutau':  ('20', '20'),
               'tautau': ('40', '40')}
    main_dir = os.path.join('/eos/user/b/bfontana/www/TriggerScaleFactors/',
                            '_'.join(('Region', *region_cuts, 'PT', *pt_cuts[args.channel], 'TURNON',
                                      met_turnon, tau_turnon)))
    if args.bigtau:
        main_dir += '_BIGTAU'
    if args.notau:
        main_dir += '_NOTAU'
    if args.nomet:
        main_dir += '_NOMET'

    regions = ('ditau', 'met', 'tau')
    met_region = ('(entries.dau2_pt < 40 and entries.dau1_pt < {}) or '.format(region_cuts[0]) +
                  '(entries.dau1_pt < 40 and entries.dau2_pt < {})'.format(region_cuts[1]))

    if args.bigtau:
        tau_region = 'entries.dau1_pt >= {} or entries.dau2_pt >= {}'.format(*region_cuts)
        ditau_region = ('entries.dau1_pt > {} and entries.dau2_pt > {} and '.format(*pt_cuts[args.channel]) +
                        'entries.dau1_pt < {} and entries.dau2_pt < {}'.format(*region_cuts))

    else:
        tau_region = ('(entries.dau2_pt < 40 and entries.dau1_pt >= {}) or '.format(region_cuts[0]) +
                      '(entries.dau1_pt < 40 and entries.dau2_pt >= {})'.format(region_cuts[1])) #this one is never realized due to the 190GeV trigger cut
        ditau_region = 'entries.dau1_pt > {} and entries.dau2_pt > {}'.format(*pt_cuts[args.channel])

    if args.notau:
        tau_region = 'False'
    if args.nomet:
        met_region = 'False'
        
    #### run main function ###
    if args.sequential:
        if not args.plot:
            for sample in args.masses:
                test_trigger_regions(args.indir, sample, args.channel)
    else:
        if not args.plot:
            pool = multiprocessing.Pool(processes=6)
            pool.starmap(test_trigger_regions,
                         zip(it.repeat(args.indir), args.masses, it.repeat(args.channel)))
    ###########################

    contam1, contam2, contam1_errors, contam2_errors = ([] for _ in range(4))
    from_directory = os.path.join(main_dir, args.channel)
    for sample in args.masses:
        outname = get_outname(sample, args.channel, region_cuts, pt_cuts[args.channel], met_turnon, tau_turnon,
                              args.bigtau, args.notau, args.nomet)
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
                              'Base && !MET', 'Base && !Tau',
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
                    # plot(jBase, jNoBase_MET, jNoBase_NoMET_Tau, *jopt)
                    
                for v in variables_2D:
                    suff = lambda x : x+'_'+reg+'_'+v[0]+'_VS_'+v[1]+'_'+cat
                    
                    hBase = f_in.Get(suff('hBase'))
                    hMET  = f_in.Get(suff('hMET'))
                    hTau  = f_in.Get(suff('hTau'))
 
                    hBase_MET = f_in.Get(suff('hBase_MET'))
                    hBase_Tau = f_in.Get(suff('hBase_Tau'))
                    hMET_Tau  = f_in.Get(suff('hMET_Tau'))

                    hNoBase_MET       = f_in.Get(suff('hNoBase_MET'))
                    hBase_NoMET       = f_in.Get(suff('hBase_NoMET'))
                    hNoBase_NoMET_Tau = f_in.Get(suff('hNoBase_NoMET_Tau'))
                    hNoBase_Tau       = f_in.Get(suff('hNoBase_Tau'))
                    hBase_NoTau       = f_in.Get(suff('hBase_NoTau'))
                    hNoBase_MET_NoTau = f_in.Get(suff('hNoBase_MET_NoTau'))
 
                    hBase_MET_Tau = f_in.Get(suff('hBase_MET_Tau'))
                    hBase_NoMET_NoTau = f_in.Get(suff('hBase_NoMET_NoTau'))

                    # plot2D(hBase, v, args.channel, sample, 'Base', cat, from_directory, reg)
                    # plot2D(hBase_NoTau, v, args.channel, sample, 'Base + !Tau', cat, from_directory, reg)
                    # #plot2D(hMET, v, args.channel, sample, 'MET', cat, from_directory, reg)
                    # plot2D(hTau, v, args.channel, sample, 'Tau', cat, from_directory, reg)

                # count only for one of the variables, it does not matter which
                nbinsx, nbinsy = hBase.GetNbinsX(), hBase.GetNbinsY()
                cBase = round(hBase.Integral(0, nbinsx+1, 0, nbinsy+1), 2)
                cMET  = round(hMET.Integral(0, nbinsx+1, 0, nbinsy+1), 2)
                cTau  = round(hTau.Integral(0, nbinsx+1, 0, nbinsy+1), 2)
                
                cBase_MET = round(hBase_MET.Integral(0, nbinsx+1, 0, nbinsy+1), 2)
                cBase_Tau = round(hBase_Tau.Integral(0, nbinsx+1, 0, nbinsy+1), 2)
                cMET_Tau  = round(hMET_Tau.Integral(0, nbinsx+1, 0, nbinsy+1), 2)
                
                cNoBase_MET       = round(hNoBase_MET.Integral(0, nbinsx+1, 0, nbinsy+1), 2)
                cBase_NoMET       = round(hBase_NoMET.Integral(0, nbinsx+1, 0, nbinsy+1), 2)
                cNoBase_NoMET_Tau = round(hNoBase_NoMET_Tau.Integral(0, nbinsx+1, 0, nbinsy+1), 2)
                cNoBase_Tau       = round(hNoBase_Tau.Integral(0, nbinsx+1, 0, nbinsy+1), 2)
                cBase_NoTau       = round(hBase_NoTau.Integral(0, nbinsx+1, 0, nbinsy+1), 2)
                cNoBase_MET_NoTau = round(hNoBase_MET_NoTau.Integral(0, nbinsx+1, 0, nbinsy+1), 2)
                
                cBase_MET_Tau = round(hBase_MET_Tau.Integral(0, nbinsx+1, 0, nbinsy+1), 2)
                cBase_NoMET_NoTau = round(hBase_NoMET_NoTau.Integral(0, nbinsx+1, 0, nbinsy+1), 2)

                # if reg=='ditau':
                #     subsets = (cBase, cMET, cBase_MET, cTau, cBase_Tau, cMET_Tau, cBase_MET_Tau)
                #     plt.figure(figsize=(10,8))
                #     v = venn3(subsets=subsets, set_labels = (r'$\tau\tau$', 'MET', 'Tau'))
                #     v.get_patch_by_id('100').set_alpha(1.0)
                #     c = venn3_circles(subsets=subsets, linestyle='dashed')
                #     plt.title('Baseline selection')
                #     for ext in extensions:
                #         plt.savefig(os.path.join(out_counts[categories.index(cat)], 'venn_' + reg + '.' + ext))

                with open(os.path.join(out_counts[categories.index(cat)], 'table.csv'), 'a') as f:
                    reader = csv.writer(f, delimiter=',', quotechar='|')
                    row = [reg, cBase, cMET, cTau,
                           cBase_MET, cBase_Tau, cMET_Tau,
                           cNoBase_MET, cNoBase_NoMET_Tau,
                           cNoBase_Tau, cNoBase_MET_NoTau,
                           cBase_NoMET, cBase_NoTau,
                           cBase_MET_Tau]
                    reader.writerow(row)

                if reg=='ditau':
                    c_ditau_trg[reg] = cBase
                    c_met_trg[reg] = cNoBase_MET
                    c_tau_trg[reg] = cNoBase_NoMET_Tau
                elif reg=='met':
                    c_ditau_trg[reg] = cBase_NoMET
                    c_met_trg[reg] = cMET
                    c_tau_trg[reg] = cNoBase_NoMET_Tau
                elif reg=='tau':
                    c_ditau_trg[reg] = cBase_NoTau
                    c_met_trg[reg] = cNoBase_MET_NoTau
                    c_tau_trg[reg] = cTau
                        
        f_in.Close()

        text = {'mass': sample,
                'out': os.path.join(out_counts[categories.index(cat)],
                                    'diagram.html')}
        sq_res = square_diagram(c_ditau_trg, c_met_trg, c_tau_trg, args.channel,
                                region_cuts, pt_cuts[args.channel], text=text,
                                bigtau=args.bigtau, notau=args.notau, nomet=args.nomet)
        c1, c2, e1, e2 = sq_res
        contam1.append(c1)
        contam2.append(c2)
        contam1_errors.append(e1)
        contam2_errors.append(e2)

    contam1 = [float(x) for x in contam1]
    contam2 = [float(x) for x in contam2]
    contam1_errors  = [float(x) for x in contam1_errors]
    contam2_errors  = [float(x) for x in contam2_errors]
    masses = [float(x) for x in args.masses]
    contamination_save('data', '_'.join(region_cuts) + '_' + args.channel,
                       contam1, contam2, contam1_errors, contam2_errors, masses)

    if args.copy:
        import subprocess
        to_directory = os.path.join('/eos/user/b/bfontana/www/TriggerScaleFactors', main_dir)
        to_directory = os.path.join(to_directory, args.channel)
     
        for sample in args.masses:
            sample_from = os.path.join(from_directory, sample)
            print('Copying: {}\t\t--->\t{}'.format(sample_from, to_directory), flush=True)
            subprocess.run(['rsync', '-ah', sample_from, to_directory])
        
    print('Done.')
