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
from collections import defaultdict as dd
import importlib

import inclusion
from inclusion import selection
from inclusion.config import main
from inclusion.utils import utils

import ROOT
import hist
import pickle

from bokeh.plotting import figure, output_file, save
from bokeh.models import Range1d, Label

tau = '\u03C4'
mu  = '\u03BC'
pm  = '\u00B1'
ditau = tau+tau

def contamination_save(savepath, label, c1, c2, e1, e2, m, mode='w'):
    with h5py.File(os.path.join(savepath, label + '.hdf5'), mode) as f:
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

def stats_save(savepath, label, c1, e1, m, mode='w'):
    with h5py.File(os.path.join(savepath, label + '.hdf5'), mode) as f:
        dset = f.create_dataset(label+'_stats', (3,len(m)), dtype='f')
        dset[0, :] = m
        dset[1, :] = c1
        dset[2, :] = e1
        dset.cols = ['mass [GeV]',
                     'stats',
                     'uncertainty']
        
def rec_dd():
    return dd(rec_dd)

def square_diagram(c_legacy_trg, c_met_trg, c_tau_trg, channel,
                   ptcuts, text, notext=False, bigtau=False, notau=False, nomet=False):
    base = {'etau': 'e+e'+tau, 'mutau': mu+'+'+mu+tau, 'tautau': ditau}
    output_file(text['out'])
    print('Saving file {}'.format(text['out']))

    topr = 9.9
    shft = 0.1
    start, b1, b2, b3 = 0.0, 2, 2.5, 7.5
    xgap = b1+shft# else b2+shft
    
    p = figure(title='m(X)={}GeV'.format(text['mass']), width=600, height=400,
               tools='save')
    p.x_range = Range1d(0, 11.5)
    p.y_range = Range1d(0, 10)
    p.outline_line_color = None
    p.toolbar.logo = None
    p.xgrid.grid_line_color = None
    p.ygrid.grid_line_color = None
    p.xaxis.ticker = [b1,b3]#else [b1, b2, b3]
    p.xaxis.major_label_overrides = ({b1: str(ptcuts[0]), b3: str(regcuts[0])})

    for aaxis in (p.xaxis, p.yaxis):
        aaxis.axis_label_text_font_size = "12pt"
        aaxis.axis_label_text_font_style = "normal"
        aaxis.major_label_text_font_size = "11pt"

    p.yaxis.axis_label_standoff = 0
    p.yaxis.ticker = [b1, b3]# else [b1, b2, b3]
    if len(ptcuts) > 1:
        p.yaxis.major_label_overrides = ({b1: str(ptcuts[0]), b3: str(regcuts[1])})
    else:
        p.yaxis.major_label_overrides = ({b1: str(ptcuts[0])})

    p.xaxis.axis_label = r'\(p_T(\tau_1)  [GeV]\)'
    p.yaxis.axis_label = r'\(p_T(\tau_2)  [GeV]\)'

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

    gain = (100*float(c_met_trg['legacy']+c_tau_trg['legacy']) /
            (c_legacy_trg['legacy']+c_met_trg['legacy']+c_tau_trg['legacy']))
    gain = str(round(gain,2))
    if not notext:
        stats_ditau = {'legacy': Label(x=b1+0.3, y=b1+1.5, text=ditau+': '+str(c_legacy_trg['legacy']),
                                       text_color='black', **label_opt),
                       'met':   Label(x=b1+0.3, y=b1+1.1, text='met && !'+ditau+': '+str(c_met_trg['legacy']),
                                      text_color='black', **label_opt),
                       'tau':   Label(x=b1+0.3, y=b1+0.7, text=tau+' && !met && !'+ditau+': '+str(c_tau_trg['legacy']),
                                      text_color='black', **label_opt),
                       'gain':  Label(x=b1+0.3, y=b1+0.3, text='Gain: '+gain+'%',
                                      text_color='blue', **label_opt),}
        for key,elem in stats_ditau.items():
            if (nomet and notau and key=='gain') or key=='gain':
                continue
            p.add_layout(elem)
 
    if not nomet:
        try:
            contam_by_tau = (100*(float(c_tau_trg['met'])+c_legacy_trg['met']) /
                             (c_met_trg['met']+c_tau_trg['met']+c_legacy_trg['met']))
            contam_by_tau = str(round(contam_by_tau,2))
        except ZeroDivisionError:
            contam_by_tau = '0'
        if not notext:
            stats_met = {'met':   Label(x=b1+0.3, y=1.3, text='met: '+str(c_met_trg['met']),
                                        text_color='black', **label_opt),
                         'tau':   Label(x=b1+0.3, y=0.5, text=tau+' && !'+ditau+' && !met: '+str(c_tau_trg['met']),
                                        text_color='black', **label_opt),
                         'legacy': Label(x=b1+0.3, y=0.9, text=ditau+' && !met: '+str(c_legacy_trg['met']),
                                         text_color='black', **label_opt),
                         'contamination': Label(x=b1+0.3, y=0.1, text='Contam.: '+contam_by_tau+'%',
                                                text_color='blue', **label_opt),}
            for key,elem in stats_met.items():
                if key != 'contamination':
                    p.add_layout(elem)
 
    if not notau:
        num_ditau = float(c_legacy_trg['tau'])
        num_both = num_ditau + float(c_met_trg['tau'])
        den = float(c_met_trg['tau']+c_tau_trg['tau']+c_legacy_trg['tau'])
        if num_ditau == 0 or num_both == 0:
            contam_ditau = '0'
            contam_both = '0'
            err_ditau = '0'
            err_both = '0'
        else:
            contam_ditau = 100*num_ditau/den
            contam_both = 100*num_both/den
 
            enum_ditau = np.sqrt(c_legacy_trg['tau'])
            eden_ditau = enum_ditau + np.sqrt(c_met_trg['tau']) + np.sqrt(c_tau_trg['tau'])
            enum_both  = np.sqrt(c_met_trg['tau']) + np.sqrt(c_legacy_trg['tau'])
            eden_both  = enum_both + np.sqrt(c_tau_trg['tau'])
            err_ditau = contam_ditau * np.sqrt(enum_ditau**2/num_ditau**2 + eden_ditau**2/den**2)
            err_both  = contam_both * np.sqrt(enum_both**2/num_both**2 + eden_both**2/den**2)
        
            contam_ditau = str(round(contam_ditau,2))
            contam_both  = str(round(contam_both,2))
            err_ditau    = str(round(err_ditau,2))
            err_both     = str(round(err_both,2))

        if not notext:
            stats_tau = {'tau':   Label(x=b3+0.2, y=1.3, text=tau+': '+str(c_tau_trg['tau']),
                                        text_color='black', **label_opt),
                         'met':   Label(x=b3+0.2, y=0.5, text='met && !'+ditau+' && !'+tau+': '+str(c_met_trg['tau']),
                                        text_color='black', **label_opt),
                         'legacy': Label(x=b3+0.2, y=0.9, text=ditau+' && !'+tau+': '+str(c_legacy_trg['tau']),
                                         text_color='black', **label_opt),
                         'contamination_both': Label(x=b3+0.2, y=0.1,
                                                     text='Contam.: ('+str(contam_both)+pm+str(err_both)+')%',
                                                     text_color='blue', **label_opt),}
            for key,elem in stats_tau.items():
                if key != 'contamination_both':
                    p.add_layout(elem)
 
    line_opt = dict(color='black', line_dash='dashed', line_width=2)
    p.line(x=[b1,b1], y=[start, topr+shft], **line_opt)
    p.line(x=[start,topr+shft], y=[b1,b1], **line_opt)
    p.line(x=[b3,b3], y=[start,topr+shft], **line_opt)
    p.line(x=[start,topr+shft], y=[b3,b3], **line_opt)
     
    p.output_backend = 'svg'
    save(p)
    return contam_ditau, contam_both, err_ditau, err_both
    
def get_outname(sample, channel, regcuts, ptcuts, met_turnon, bigtau, notau, nomet):
    utils.create_single_dir('data')

    name = sample + '_' + channel + '_'
    name += '_'.join((*[str(x) for x in regcuts], 'ptcuts',
                      *[str(x) for x in ptcuts], 'turnon', str(met_turnon)))
    if bigtau:
        name += '_BIGTAU'
    if notau:
        name += '_NOTAU'
    if nomet:
        name += '_NOMET'
    name += '.pkl'

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
    leg2.AddEntry(hbase2, 'legacy')
    leg2.AddEntry(hmet2, '!legacy && MET')
    leg2.AddEntry(htau2, '!legacy && !MET && tau')

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

def which_region(ent, year, ptcuts, regcuts, channel, met_turnon, bigtau=False, notau=False, nomet=False):
    """
    Select one of the threee non-overlapping regions: leg(acy), met and tau
    - ptcuts: pT cuts, delimiting the three regions based on the HLT pT thresholds
    - regcuts: region cuts, delimiting the three regions in the pT space
    Example in tautau:
    - ptcuts = (40, 40)
    - regcuts = (190, 190)
    """
    dau1_eta, dau2_eta = ent.dau1_eta <= 2.1, ent.dau2_eta <= 2.1
    
    if bigtau:
        if channel == "tautau":
            leg = (ent.dau1_pt >= ptcuts[0] and ent.dau2_pt >= ptcuts[1] and
                   ent.dau1_pt < regcuts[0] and ent.dau2_pt < regcuts[1] and
                   dau1_eta and dau2_eta)
            tau = ((ent.dau1_pt >= regcuts[0] and dau1_eta) or
                   (ent.dau2_pt >= regcuts[1] and dau2_eta))
 
        elif channel == "etau" and year == "2016":
            leg = ent.dau1_pt >= ptcuts[0] and ent.dau2_pt < regcuts[1]
            tau = ent.dau2_pt >= regcuts[1] and dau2_eta
 
        else: #mutau or etau non-2016
            leg = (ent.dau1_pt >= ptcuts[0] or (ent.dau1_pt >= ptcuts[1] and ent.dau2_pt >= ptcuts[2]) and
                   ent.dau2_pt < regcuts[1] and dau2_eta)
            tau = ent.dau2_pt >= regcuts[1] and dau2_eta
            
    else:
        if channel == "tautau":
            leg = (ent.dau1_pt >= ptcuts[0] and ent.dau2_pt >= ptcuts[1] and
                   dau1_eta and dau2_eta)
            tau = ((ent.dau2_pt < ptcuts[1] and ent.dau1_pt >= regcuts[0] and dau1_eta) or 
                   (ent.dau1_pt < ptcuts[0] and ent.dau2_pt >= regcuts[1] and dau2_eta))
 
        elif channel == "etau" and year == "2016":
            leg = ent.dau1_pt >= ptcuts[0]
            tau = ent.dau1_pt < ptcuts[0] and ent.dau2_pt >= regcuts[1] and dau2_eta
 
        else: #mutau or etau non-2016
            leg = ent.dau1_pt >= ptcuts[0] or (ent.dau1_pt >= ptcuts[1] and ent.dau2_pt >= ptcuts[2]) and dau2_eta
            tau = ent.dau1_pt < ptcuts[1] and ent.dau2_pt >= regcuts[1] and dau2_eta

    met = ent.metnomu_et > met_turnon and not leg and not tau
        
    # only one True: non-overlapping regions
    assert int(leg)+int(met)+int(tau)<=1

    if notau:
        tau = False
    if nomet:
        met = False

    return leg, met, tau


def test_trigger_regions(indir, sample, channel, spin, year, deltaR):
    outname = get_outname(sample, channel, regcuts, ptcuts, met_turnon,
                          args.bigtau, args.notau, args.nomet)
    config_module = importlib.import_module(args.configuration)
    if channel == 'etau' or channel == 'mutau':
        iso1 = (24, 0, 8)
    elif channel == 'tautau':
        iso1 = (20, 0.965, 1.005)
    binning.update({'HHKin_mass': (20, float(sample)-300, float(sample)+300),
                    'dau1_iso': iso1})
    
    full_sample = {0: 'GluGluToRadionToHHTo2B2Tau_M-' + sample + '_',
                   2: 'GluGluToBulkGravitonToHHTo2B2Tau_M-' + sample + '_'}[spin]
    
    norphans, ntotal = ({k:0 for k in categories} for _ in range(2))
    
    ahistos = rec_dd()
    for reg in regions:                
        for cat in categories:
            for htype in htypes:
                ahistos[htype][reg][cat] = hist.Hist(
                    hist.axis.Regular(*binning["dau1_pt"], name="dau1pt"),
                    hist.axis.Regular(*binning["dau2_pt"], name="dau2pt"),
                )

    t_in = ROOT.TChain('HTauTauTree')
    glob_files = glob.glob( os.path.join(indir, full_sample, 'output_*.root') )
    if len(glob_files) < 1:
        raise RuntimeError("No files!")
    for f in glob_files:
        t_in.Add(f)
    t_in.SetBranchStatus('*', 0)
    _entries = utils.define_used_tree_variables(cut=config_module.custom_cut)
    _entries += tuple(variables)
    for ientry in _entries:
        t_in.SetBranchStatus(ientry, 1)
  
    for entry in t_in:
        # this is slow: do it once only
        entries = utils.dot_dict({x: getattr(entry, x) for x in _entries})

        in_legacy, in_met, in_tau = which_region(entries, year, ptcuts, regcuts, channel, met_turnon,
                                                 bigtau=args.bigtau, notau=args.notau, nomet=args.nomet)
        
        sel = selection.EventSelection(entries, isdata=False, configuration=config_module)

        pass_trg = sel.pass_triggers(triggers[channel]) and entries.isLeptrigger
        pass_met = sel.pass_triggers(('METNoMu120',))
        pass_tau = sel.pass_triggers(('IsoTau180',))

        cuts = {'Base'           : pass_trg,
                'MET'            : pass_met,
                'Tau'            : pass_tau,
                'VBF'            : False,
                'BaseMET'        : pass_trg and pass_met,
                'BaseTau'        : pass_trg and pass_tau,
                'METTau'         : pass_met and pass_tau,
                'NoBaseMET'      : not pass_trg and pass_met,
                'BaseNoMET'      : pass_trg and not pass_met,
                'NoBaseNoMETTau' : not pass_trg and not pass_met,
                'NoBaseTau'      : not pass_trg and pass_tau,
                'BaseNoTau'      : pass_trg and not pass_tau,
                'NoBaseMETNoTau' : not pass_trg and pass_met and not pass_tau,
                'BaseMETTau'     : pass_trg and pass_met and pass_tau,
                'BaseNoMETNoTau' : pass_trg and not pass_met and not pass_tau,
                'LegacyKin'      : pass_trg and in_legacy,
                'METKin'         : pass_met and in_met,
                'TauKin'         : pass_tau and in_tau,
                'VBFKin'         : False
                }
        assert htypes == list(cuts.keys())
        
        w_mc     = entries.MC_weight
        w_pure   = entries.PUReweight
        w_l1pref = entries.L1pref_weight
        w_trig   = entries.trigSF
        w_idiso  = entries.IdSF_deep_2d
        w_jetpu  = entries.PUjetID_SF
        w_btag   = entries.bTagweightReshape
        
        if utils.is_nan(w_mc)     : w_mc=1
        if utils.is_nan(w_pure)   : w_pure=1
        if utils.is_nan(w_l1pref) : w_l1pref=1
        if utils.is_nan(w_trig)   : w_trig=1
        if utils.is_nan(w_idiso)  : w_idiso=1
        if utils.is_nan(w_jetpu)  : w_jetpu=1
        if utils.is_nan(w_btag)   : w_btag=1
  
        evt_weight = w_mc * w_pure * w_l1pref * w_trig * w_idiso * w_jetpu * w_btag

        if utils.is_channel_consistent(channel, entries.pairType):
            if not sel.selection_cuts(lepton_veto=True, bjets_cut=True,
                                      mass_cut=config_module.mass_cut):
                continue

            for cat in categories:
                if sel.sel_category(cat) and entries.ditau_deltaR > deltaR:
                    ntotal[cat] += 1

                    if in_met:
                        reg = 'met'
                    elif in_tau:
                        reg = 'tau'
                    elif in_legacy:
                        reg = 'legacy'
                    else:
                        norphans[cat] += 1
                        continue
                    assert reg in regions

                    for key,cut in cuts.items():
                        if cut:
                            ahistos[key][reg][cat].fill(dau1pt=entries["dau1_pt"],
                                                        dau2pt=entries["dau2_pt"], weight=evt_weight)

    with open(outname, "wb") as f:
        pickle.dump(ahistos, f)
                
    for cat in categories:
        orph_frac = float(norphans[cat])/ntotal[cat]
        if orph_frac > 0.1:
            print('{}% orphans ({}/{}). This is unusual.'.format(orph_frac, norphans[cat], ntotal[cat]))
            print(met_region)
            print(tau_region)
            print(legacy_region)
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

    categories = ('baseline',) #('baseline', 's1b1jresolvedMcut', 's2b0jresolvedMcut', 'sboostedLLMcut')

    htypes = ['Base', 'MET', 'Tau', 'VBF', 'BaseMET', 'BaseTau', 'METTau', 'NoBaseMET', 'BaseNoMET',
              'NoBaseNoMETTau', 'NoBaseTau', 'BaseNoTau', 'NoBaseMETNoTau', 'BaseMETTau', 'BaseNoMETNoTau',
              'LegacyKin', 'METKin', 'TauKin', 'VBFKin']
    
    # Parse input arguments
    desc = 'Producer trigger histograms.\n'
    desc += "Run example: python tests/test_trigger_regions.py --indir /data_CMS/cms/alves/HHresonant_SKIMS/SKIMS_UL18_EOSv4_Signal/ --masses 400 500 600 700 800 900 1000 1250 1500 --channels ETau --met_turnon 180 --region_cuts 40 40 --copy"
    parser = argparse.ArgumentParser(description=desc, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('--indir', required=True, type=str,
                        help='Full path of ROOT input file')
    parser.add_argument('--masses', required=True, nargs='+', type=str,
                        help='Resonance mass')
    parser.add_argument('--channel', required=True, type=str,  
                        help='Select the channel over which the workflow will be run.' )
    parser.add_argument('--year', required=True, type=str, choices=('2016', '2017', '2018'),
                        help='Select the year over which the workflow will be run.' )
    parser.add_argument('--spin', required=True, type=int, choices=(0, 2), 
                        help='Select the spin hypothesis over which the workflow will be run.' )
    parser.add_argument('--deltaR', type=float, default=0.5, help='DeltaR between the two leptons.' )
    parser.add_argument('--plot', action='store_true',
                        help='Reuse previously produced data for quick plot changes.')
    parser.add_argument('--copy', action='store_true',
                        help='Copy the outputs to EOS at the end.')
    parser.add_argument('--notext', action='store_true', help='Square diagram without text.')
    parser.add_argument('--sequential', action='store_true',
                        help='Do not use the multiprocess package.')
    parser.add_argument('--bigtau', action='store_true',
                        help='Consider a larger single tau region, reducing the ditau one.')
    parser.add_argument('--notau', action='store_true',
                        help='Remove the single tau region (default analysis).')
    parser.add_argument('--nomet', action='store_true',
                        help='Remove the MET region (default analysis).')
    parser.add_argument('--met_turnon', type=float, default=180,
                        help='MET trigger turnon cut [GeV].' )
    parser.add_argument('--region_cuts', required=False, type=float, nargs=2, default=(190, 190),
                        help='High/low regions pT1 and pT2 selection cuts [GeV].' )
    parser.add_argument('--configuration', dest='configuration', required=True,
                        help='Name of the configuration module to use.')
    args = utils.parse_args(parser)

    met_turnon = args.met_turnon
    regcuts = args.region_cuts
    ptcuts = utils.get_ptcuts(args.channel, args.year)
        
    main_dir = os.path.join('/eos/home-b/bfontana/www/TriggerScaleFactors/',
                            '_'.join(('Spin' + str(args.spin), args.channel, *[str(x) for x in regcuts],
                                      'DR', str(args.deltaR), 'PT', *[str(x) for x in ptcuts], 'TURNON',
                                      str(met_turnon))))
    if args.bigtau:
        main_dir += '_BIGTAU'
    if args.notau:
        main_dir += '_NOTAU'
    if args.nomet:
        main_dir += '_NOMET'

    regions = ('legacy', 'met', 'tau')
        
    #### run main function ###
    if not args.plot:
        if args.sequential:
            for sample in args.masses:
                test_trigger_regions(args.indir, sample, args.channel, args.spin, args.year, args.deltaR)
        else:
            pool = multiprocessing.Pool(processes=6)
            pool.starmap(test_trigger_regions,
                         zip(it.repeat(args.indir), args.masses, it.repeat(args.channel),
                             it.repeat(args.spin), it.repeat(args.year), it.repeat(args.deltaR)))

    ###########################

    sum_stats, err_sum_stats = ([] for _ in range(2))
    contam1, contam2, contam1_errors, contam2_errors = ([] for _ in range(4))
    from_directory = os.path.join(main_dir, args.channel)
    for sample in args.masses:
        outname = get_outname(sample, args.channel,
                              (str(x) for x in regcuts), (str(x) for x in ptcuts),
                              str(met_turnon), args.bigtau, args.notau, args.nomet)
        with open(outname, "rb") as f:
            ahistos = pickle.load(f)

        # write csv header, one per category
        out_counts = []
        for cat in categories:
            out_counts.append( os.path.join(from_directory, sample, cat, 'counts') )
            utils.create_single_dir(out_counts[-1])
            with open(os.path.join(out_counts[-1], 'table.csv'), 'w') as f:
                reader = csv.writer(f, delimiter=',', quotechar='|')
                header_row = ['Region']
                header_row.extend(htypes)
                reader.writerow(header_row)

        # plot histograms and fill CSV with histogram integrals
        c_legacy_trg, c_met_trg, c_tau_trg = ({} for _ in range(3))

        acounts = rec_dd()
        for reg in regions:
            for key,cut in ahistos.items():
                acounts[key][reg] = round(ahistos[key][reg]["baseline"].values().sum(), 2)
            
            # append to table, one line per region
            with open(os.path.join(out_counts[categories.index(cat)], 'table.csv'), 'a') as f:
                reader = csv.writer(f, delimiter=',', quotechar='|')
                row = [reg]
                row.extend([acounts[k][reg] for k in acounts.keys()])
                reader.writerow(row)

            if reg=='legacy':
                c_legacy_trg[reg] = acounts["Base"][reg]
                c_met_trg[reg]    = acounts["NoBaseMET"][reg]
                c_tau_trg[reg]    = acounts["NoBaseNoMETTau"][reg]
            elif reg=='met':
                c_legacy_trg[reg] = acounts["BaseNoMET"][reg]
                c_met_trg[reg]    = acounts["MET"][reg]
                c_tau_trg[reg]    = acounts["NoBaseNoMETTau"][reg]
            elif reg=='tau':
                c_legacy_trg[reg] = acounts["BaseNoTau"][reg]
                c_met_trg[reg]    = acounts["NoBaseMETNoTau"][reg]
                c_tau_trg[reg]    = acounts["Tau"][reg]

        text = {'mass': sample,
                'out': os.path.join(out_counts[categories.index(cat)], 'diagram.html')}
        sq_res = square_diagram(c_legacy_trg, c_met_trg, c_tau_trg, args.channel,
                                [str(x) for x in ptcuts], text=text, notext=args.notext,
                                bigtau=args.bigtau, notau=args.notau, nomet=args.nomet)
        c1, c2, e1, e2 = sq_res
        contam1.append(c1)
        contam2.append(c2)
        contam1_errors.append(e1)
        contam2_errors.append(e2)

        stats_l = [c_legacy_trg['legacy'],c_met_trg['legacy'],c_tau_trg['tau'],c_legacy_trg['tau']]
        stats = sum(stats_l)
        estats = sum(np.sqrt(stats_l))
        sum_stats.append(stats)
        err_sum_stats.append(estats)

    contam1 = [float(x) for x in contam1]
    contam2 = [float(x) for x in contam2]
    contam1_errors  = [float(x) for x in contam1_errors]
    contam2_errors  = [float(x) for x in contam2_errors]
    masses = [float(x) for x in args.masses]
    contamination_save('data', '_'.join([str(x) for x in regcuts]) + '_' + args.channel,
                       contam1, contam2, contam1_errors, contam2_errors, masses)
    stats_save('data', '_'.join([str(x) for x in regcuts]) + '_' + args.channel,
               sum_stats, err_sum_stats, masses, mode='a')

    if args.copy:
        import subprocess
        to_directory = os.path.join('/eos/home-b/bfontana/www/TriggerScaleFactors', main_dir)
        to_directory = os.path.join(to_directory, args.channel)
     
        for sample in args.masses:
            sample_from = os.path.join(from_directory, sample)
            print('Copying: {}\t\t--->\t{}'.format(sample_from, to_directory), flush=True)
            subprocess.run(['rsync', '-ah', sample_from, to_directory])
