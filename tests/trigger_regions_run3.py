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
import uproot

from bokeh.plotting import figure, output_file, save
from bokeh.models import Range1d, Label

tau = '\u03C4'
mu  = '\u03BC'
pm  = '\u00B1'
ditau = tau+tau

def contamination_save(savepath, label, c1, c2, e1, e2, mode='w'):
    with h5py.File(os.path.join(savepath, label + '.hdf5'), mode) as f:
        dset = f.create_dataset(label, (4, len(c1)), dtype='f')
        dset[0, :] = c1
        dset[1, :] = c2
        dset[2, :] = e1
        dset[3, :] = e2
        dset.cols = ['contamination ' + ditau + '(%)',
                     'contamination ' + ditau + ' MET (%)',
                     'uncertainty' + ditau,
                     'uncertainty' + ditau + ' MET']

def stats_save(savepath, label, c1, e1, mode='w'):
    with h5py.File(os.path.join(savepath, label + '.hdf5'), mode) as f:
        dset = f.create_dataset(label+'_stats', (2, len(c1)), dtype='f')
        dset[0, :] = c1
        dset[1, :] = e1
        dset.cols = ['stats',
                     'uncertainty']
        
def rec_dd():
    return dd(rec_dd)

def square_diagram(c_legacy_trg, c_met_trg, c_tau_trg, channel,
                   ptcuts, text, notext=False, bigtau=False):
    base = {'etau': 'e+e'+tau, 'mutau': mu+'+'+mu+tau, 'tautau': ditau}
    output_file(text['out'])
    print('Saving file {}'.format(text['out']))

    topr = 9.9
    shft = 0.1
    start, b1, b2, b3 = 0.0, 2, 2.5, 7.5
    xgap = b1+shft# else b2+shft
    
    p = figure(title='m(X)=0GeV', width=600, height=400,
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
            p.add_layout(elem)
 

    try:
        if c_met_trg['met']+c_tau_trg['met']+c_legacy_trg['met'] == 0.0:
            raise ZeroDivisionError
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
    
def get_outname(channel, bigtau):
    utils.create_single_dir('data')
    
    name = ""
    if bigtau:
        name += '_BIGTAU'
    name += '_all.root'

    s = 'data/regions_preEE_12p11-old-way-pf75{}'.format(name)
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

def match_trigger_object(off_eta, off_phi, obj_id, TrigObj_id, TrigObj_filterBits, TrigObj_eta, TrigObj_phi, bits):
    for iobj in range(len(TrigObj_id)):
        if TrigObj_id[iobj] != obj_id:
            continue
        dPhi = off_phi - TrigObj_phi[iobj]
        dEta = off_eta - TrigObj_eta[iobj]
        delR2 = dPhi * dPhi + dEta * dEta
        if delR2 > 0.25: #0.5 * 0.5
            continue
        # matched_bits = True
        # for bit in bits:
        #     if (TrigObj_filterBits[iobj] & (1<<bit)) == 0: 
        #         matched_bits = False
        #         break
        #     if not matched_bits:
        #         continue
        return True
    return False

def match_trigger_objects_QuadJet(entry):
    TrigObj_id, TrigObj_filterBits = entry['TrigObj_id'], entry['TrigObj_filterBits'], 
    TrigObj_eta, TrigObj_phi =  entry['TrigObj_eta'], entry['TrigObj_phi']
    dau1_eta, dau1_phi = entry['dau1_eta'], entry['dau1_phi']
    dau2_eta, dau2_phi = entry['dau2_eta'], entry['dau2_phi']
    bjet1_eta, bjet1_phi = entry['Jet_eta']['bjet1_JetIdx'], entry['Jet_phi']['bjet1_JetIdx']
    bjet2_eta, bjet2_phi = entry['Jet_eta']['bjet2_JetIdx'], entry['Jet_phi']['bjet2_JetIdx']
    pass_dau1 = match_trigger_object(dau1_eta, dau1_phi, 1, TrigObj_id, TrigObj_filterBits, TrigObj_eta, TrigObj_phi, 0)
    pass_dau2 = match_trigger_object(dau2_eta, dau2_phi, 1, TrigObj_id, TrigObj_filterBits, TrigObj_eta, TrigObj_phi, 0)
    pass_bjet1 = match_trigger_object(bjet1_eta, bjet1_phi, 1, TrigObj_id, TrigObj_filterBits, TrigObj_eta, TrigObj_phi, 0)
    pass_bjet2 = match_trigger_object(bjet2_eta, bjet2_phi, 1, TrigObj_id, TrigObj_filterBits, TrigObj_eta, TrigObj_phi, 0)
    return pass_dau1 and pass_dau2 and pass_bjet1 and pass_bjet2

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

def which_region(ent, year, ptcuts, regcuts, channel, sel, bigtau=False):
    """
    Select one of the threee non-overlapping regions: leg(acy), met and tau
    - ptcuts: pT cuts, delimiting the three regions based on the HLT pT thresholds
    - regcuts: region cuts, delimiting the three regions in the pT space
    Example in tautau:
    - ptcuts = (40, 40)
    - regcuts = (190, 190)

    Eta cuts are defined with different strategies:
    - tautau: eta trigger cut at 2.1, so that the MET trigger can take advantage of high-eta events for the legacy region
    - e/mutau: eta cuts correspond to the analysis selection
    Since different triggers in the etau and mutau channels can require different eta thresholds
    (eg: single-muon has no threshold and cross-muon-tau has one at 2.1), a fully correct approach
    would have to apply different cuts depending on the fired trigger bit. To avoid
    this extra complication, which would very likely bring a negligible signal acceptance improvement,
    the eta cuts applied correspond to the object selection in the analysis.
    """
    # eta1_sel = {"etau": abs(ent.dau1_eta) <= 2.5, "mutau": abs(ent.dau1_eta) <= 2.4}
    # eta1_trg = abs(ent.dau1_eta) <= 2.1
    # eta2_trg = abs(ent.dau2_eta) <= 2.1
    
    # if bigtau:
    #     if channel == "tautau":
    #         tau = (ent.dau1_pt >= regcuts[0] and eta1_trg) or (ent.dau2_pt >= regcuts[1] and eta2_trg)
    #         leg = ent.dau1_pt >= ptcuts[0] and ent.dau2_pt >= ptcuts[1] and eta1_trg and eta2_trg and not tau
            
    #     elif channel == "etau" and year == "2016":
    #         tau = ent.dau2_pt >= regcuts[1] and eta2_trg
    #         leg = ent.dau1_pt >= ptcuts[0] and eta1_trg and not tau
 
    #     else: #mutau or etau non-2016
    #         single_lepton_validity = ent.dau1_pt >= ptcuts[0] and eta1_sel[channel]
    #         cross_lepton_validity  = ent.dau1_pt >= ptcuts[1] and eta1_trg and ent.dau2_pt >= ptcuts[2] and eta2_trg

    #         tau = ent.dau2_pt >= regcuts[1] and eta2_trg
    #         leg = (single_lepton_validity or cross_lepton_validity) and not tau
            
    if channel == "tautau":
        leg = True 
        tau = False
        # leg = ent.dau1_pt >= ptcuts[0] and ent.dau2_pt >= ptcuts[1] and eta1_trg and eta2_trg
        # tau = ((ent.dau1_pt >= regcuts[0] and eta1_trg) or (ent.dau2_pt >= regcuts[1] and eta2_trg)) and not leg

    elif channel == "etau" and year == "2016":
        leg = ent.dau1_pt >= ptcuts[0] and eta1_trg
        tau = ent.dau2_pt >= regcuts[1] and eta2_trg and not leg

    else: #mutau or etau non-2016
        leg = True
        tau = False
        # leg = sel.check_bit(main.trig_map[year]['IsoMu24']['mc']) or sel.check_bit(main.trig_map[year]['IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1']['mc'])
        # tau = not leg and sel.check_bit(main.trig_map[year]['LooseDeepTauPFTauHPS180_L2NN_eta2p1']['mc'])
        # single_lepton_validity = ent.dau1_pt >= ptcuts[0] and eta1_sel[channel]
        # cross_lepton_validity  = ent.dau1_pt >= ptcuts[1] and eta1_trg and ent.dau2_pt >= ptcuts[2] and eta2_trg
        # single_tau_validity    = ent.dau2_pt >= regcuts[1] and eta2_trg

        # leg = single_lepton_validity or cross_lepton_validity
        # tau = single_tau_validity and not leg

    met = not leg and not tau and sel.check_bit(main.trig_map[year]['PFMETNoMu120_PFMHTNoMu120_IDTight']['mc'])
        
    # only one True: non-overlapping regions
    assert int(leg)+int(met)+int(tau)<=1

    return leg, met, tau

def pass_offline_selection(triggers, entry):
    res = False
    for trigger in triggers:
        trigger_pass = entry[trigger]
        if not trigger_pass:
            continue
        if trigger == "HLT_IsoMu24" and entry['pairType'] == 0:
            res = res or (entry['dau1_pt'] > 25 and abs(entry['dau1_eta']) < 2.4)
        elif trigger == "HLT_IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1" and entry['pairType'] == 0:
            res = res or (entry['dau1_pt'] > 21 and entry['dau2_pt'] > 32)
        elif trigger == "HLT_Ele30_WPTight_Gsf" and entry['pairType'] == 1:
            res = res or (entry['dau1_pt'] > 32 and abs(entry['dau1_eta']) < 2.4)
        elif trigger == "HLT_Ele24_eta2p1_WPTight_Gsf_LooseDeepTauPFTauHPS30_eta2p1_CrossL1" and entry['pairType'] == 1:
            res = res or (entry['dau1_pt'] > 26 and entry['dau2_pt'] > 35)
        elif trigger == "HLT_DoubleMediumDeepTauPFTauHPS35_L2NN_eta2p1" and entry['pairType'] == 2:
            res = res or (entry['dau1_pt'] > 40 and entry['dau2_pt'] > 40)
        elif trigger == "HLT_DoubleTightChargedIsoPFTauHPS35_Trk1_eta2p1" and entry['pairType'] == 2:
            res = res or (entry['dau1_pt'] > 40 and entry['dau2_pt'] > 40)
        elif trigger == "HLT_LooseDeepTauPFTauHPS180_L2NN_eta2p1":
            if entry['pairType'] == 2:
                res = res or (entry['dau1_pt'] > 190 or entry['dau2_pt'] > 190)
            else:
                res = res or entry['dau2_pt'] > 190
        elif trigger == "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight":
            res = res or entry["MET_pt"] > 150
        elif trigger == "HLT_Mu50":
            res = res or (entry['dau1_pt'] > 51 and abs(entry['dau1_eta']) < 2.4)
        elif trigger == "HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet60" and entry['pairType'] == 2:
            res = res or (entry['dau1_pt'] > 35 and entry['dau2_pt'] > 35)
        elif trigger == "HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet75" and entry['pairType'] == 2:
            res = res or (entry['dau1_pt'] > 35 and entry['dau2_pt'] > 35 and (entry['bjet1_pt'] > 75 or entry['bjet2_pt'] > 75))
        elif trigger == "HLT_QuadPFJet70_50_40_30_PFBTagParticleNet_2BTagSum0p65":
            if entry['pairType'] == 2:
                res = res or (abs(entry['dau1_eta']) < 2.1 and abs(entry['dau2_eta']) < 2.1)
            else:
                res = res or (abs(entry['dau1_eta']) < 2.4 and abs(entry['dau2_eta']) < 2.1)
    return res


def pass_triggers(triggers, entry):
    res = False
    for trigger in triggers:
        res = res or entry[trigger]
    return res


def trigger_regions(indir, channel, year, deltaR):
    outname = get_outname(channel,
                          args.bigtau)
    config_module = importlib.import_module(args.configuration)


    if channel == 'etau' or channel == 'mutau':
        iso1 = range(0, 8.1, 8/24)
    elif channel == 'tautau':
        iso1 = range(0, 401, 8)
    pNet_dist = [x/12. for x in range(13)]
    binning.update({
                    'genHH_mass': ([250, 300, 350, 400, 450, 500, 550, 600, 675, 800, 1000, 1600],),
                    'dau1_iso': (iso1,),
                    'dau1_pt': ([0, 20, 30, 40, 60, 80, 100, 125, 150, 200, 250],),
                    'dau2_iso': (iso1,),
                    'dau2_pt': ([0, 20, 30, 40, 60, 80, 100, 125, 150, 200, 250],),
                    'dau1_eta': ([-3, -2.5, -2.1, -1.8, -1.5, -1.2, -0.8, -0.4, 0, 0.4, 0.8, 1.2, 1.5, 1.8, 2.1, 2.5, 3],),
                    'dau2_eta': ([-3, -2.5, -2.1, -1.8, -1.5, -1.2, -0.8, -0.4, 0, 0.4, 0.8, 1.2, 1.5, 1.8, 2.1, 2.5, 3],),
                    'dau1_tauIdVSjet': ([0, 1, 2, 3, 4, 5, 6, 7], ),
                    'dau2_tauIdVSjet': ([0, 1, 2, 3, 4, 5, 6, 7], ),
                    'bjet1_pNet': (pNet_dist,),
                    'bjet2_pNet': (pNet_dist,),
                    'bjet1_pt': ([0, 20, 30, 40, 60, 80, 100, 125, 150, 200, 250],),
                    'bjet2_pt': ([0, 20, 30, 40, 60, 80, 100, 125, 150, 200, 250],),
                    'bjet1_eta': ([-3, -2.5, -2.1, -1.8, -1.5, -1.2, -0.8, -0.4, 0, 0.4, 0.8, 1.2, 1.5, 1.8, 2.1, 2.5, 3],),
                    'bjet2_eta': ([-3, -2.5, -2.1, -1.8, -1.5, -1.2, -0.8, -0.4, 0, 0.4, 0.8, 1.2, 1.5, 1.8, 2.1, 2.5, 3],),
                    })
    
    # full_sample = {0: 'GluGluToRadionToHHTo2B2Tau_M-' + sample + '_',
    #                2: 'GluGluToBulkGravitonToHHTo2B2Tau_M-' + sample + '_'}[spin]
    
    norphans, ntotal = ({k:0 for k in categories} for _ in range(2))
    
    ahistos = rec_dd()
    for cat in categories:                
        for chn in ['mutau', 'etau', 'tautau']:
            for htype in htypes:
                for i in binning.keys():
                    ahistos[htype][cat][chn][i] = (
                        hist.Hist.new.Variable(*binning[i], name=i)
                        .Weight()
                    )

    t_in = ROOT.TChain('Events')
    glob_files = glob.glob( os.path.join(indir, 'data_*.root') )
    if len(glob_files) < 1:
        raise RuntimeError("No files!")
    for f in glob_files:
        t_in.Add(f)
    t_in.SetBranchStatus('*', 0)
    _entries, mc_corrections = utils.define_used_tree_variables(cut=config_module.custom_cut)

    _entries += tuple([ x + mc_corrections[x] for x in mc_corrections])
    for ientry in _entries:
        t_in.SetBranchStatus(ientry, 1)
  
    for entry in t_in:
        # this is slow: do it once only
        entries = utils.dot_dict({x: getattr(entry, x) for x in _entries})

        if entries.pairType == -1:
            continue

        if entries.pairType == 0 or entries.pairType == 1:
            if entries.isOS != 1 or entries.dau2_tauIdVSjet < 5:
                continue
        if entries.pairType == 2:
            if entries.isOS != 1 or entries.dau2_tauIdVSjet < 5 or entries.dau1_tauIdVSjet < 5:
                continue
        
        sel = selection.EventSelection(entries, year=year, isdata=False, configuration=config_module)
        # in_legacy, in_met, in_tau = which_region(entries, year, ptcuts, regcuts, channel, sel,
        #                                          bigtau=args.bigtau)

        pass_base_mu = pass_triggers(triggers['mutau'], entries)
        pass_base_e = pass_triggers(triggers['etau'], entries)
        pass_base_tau = pass_offline_selection(triggers['tautau'], entries)
        # pass_alt_base_tau = pass_offline_selection(('HLT_DoubleTightChargedIsoPFTauHPS35_Trk1_eta2p1',), entries)
        pass_met = pass_offline_selection(('HLT_PFMETNoMu120_PFMHTNoMu120_IDTight',), entries)
        pass_tau = pass_offline_selection(('HLT_LooseDeepTauPFTauHPS180_L2NN_eta2p1',), entries)
        # pass_mu50 = pass_offline_selection(('HLT_Mu50',), entries)
        # pass_ele28 = pass_offline_selection(('HLT_Ele28_eta2p1_WPTight_Gsf_HT150', ), entries)
        pass_ttjet = pass_offline_selection(("HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet60",), entries)
        # pass_4jets_deepJet = pass_triggers(("HLT_QuadPFJet103_88_75_15_DoublePFBTagDeepJet_1p3_7p7_VBF1", ), entries)
        pass_4jets_pNet = pass_triggers(("HLT_QuadPFJet70_50_40_35_PFBTagParticleNet_2BTagSum0p65", ), entries) 
        # pass_4jets_match = match_trigger_objects_QuadJet(entries) if pass_4jets_pNet else False
        cuts = {
                'All': True,
                'BaseMu'           : pass_base_mu,
                'BaseE'          : pass_base_e,
                'BaseTau'       : pass_base_tau and entries.dau1_pt > 40 and entries.dau2_pt > 40,
                # 'AltBaseTau'    : pass_alt_base_tau,
                # 'BaseTauNoAltBaseTau': pass_base_tau and not pass_alt_base_tau,
                # 'NoBaseTauAltBaseTau': not pass_base_tau and pass_alt_base_tau,
                # 'BaseTauAltBaseTau': pass_base_tau and pass_alt_base_tau,
                # 'MET'            : pass_met,
                # 'Tau'            : pass_tau,
                # 'VBF'            : False,
                # 'BaseMuMET'        : pass_base_mu and pass_met,
                # 'BaseEMET'        : pass_base_e and pass_met,
                # 'BaseTauMET'        : pass_base_tau and pass_met,
                # 'BaseMuORMET'      : pass_base_mu or pass_met,
                # 'BaseEORMET'      : pass_base_e or pass_met,
                # 'BaseTauORMET'      : pass_base_tau or pass_met,
                # 'BaseMuTau'        : pass_base_mu and pass_tau,
                # 'BaseETau'        : pass_base_e and pass_tau,
                # 'BaseTauTau'        : pass_base_tau and pass_tau,
                # 'BaseMuORTau'      : pass_base_mu or pass_tau,
                # 'BaseEORTau'      : pass_base_e or pass_tau,
                # 'BaseTauORTau'      : pass_base_tau or pass_tau,
                # 'BaseTauORTTJet'    : pass_base_tau or pass_ttjet,
                # 'BaseTauORTauORTTJet': pass_base_tau or pass_tau or pass_ttjet,
                # 'METTau'         : pass_met and pass_tau,
                # 'METORTau'       : pass_met or pass_tau,
                # 'BaseMuORMETORTau' : pass_base_mu or pass_met or pass_tau,
                # 'BaseEORMETORTau' : pass_base_e or pass_met or pass_tau,
                # 'BaseTauORMETORTau' : pass_base_tau or pass_met or pass_tau,
                # 'BaseTauORMETORTauORTTJet': pass_base_tau or pass_tau or pass_met or pass_ttjet,
                # 'NoBaseMuMET'      : not pass_base_mu and pass_met,
                # 'BaseMuNoMET'      : pass_base_mu and not pass_met,
                # 'NoBaseEMET'      : not pass_base_e and pass_met,
                # 'BaseENoMET'      : pass_base_e and not pass_met,
                # 'NoBaseTauMET'      : not pass_base_tau and pass_met,
                # 'BaseTauNoMET'      : pass_base_tau and not pass_met,
                # 'NoBaseMuNoMETTau' : not pass_base_mu and not pass_met and pass_tau,
                # 'NoBaseENoMETTau' : not pass_base_e and not pass_met and pass_tau,
                # 'NoBaseTauNoMETTau' : not pass_base_tau and not pass_met and pass_tau,
                # 'NoBaseTauTau'      : not pass_base_tau and pass_tau,
                # 'NoBaseMuMu50'     : not pass_base_mu and pass_mu50,
                # 'NoBaseMuTau'      : not pass_base_mu and pass_tau,
                # 'NoBaseETau'      : not pass_base_e and pass_tau,
                # 'NoBaseMuNoMETMu50': not pass_base_mu and not pass_met and pass_mu50,
                # 'NoBaseMuNoTauMu50': not pass_base_mu and not pass_tau and pass_mu50,
                # 'NoBaseMuMETNoMu50': not pass_base_mu and pass_met and not pass_mu50,
                # 'NoBaseMuMETMu50': not pass_base_mu and pass_met and pass_mu50,
                # 'BaseTauNoTau'      : pass_base_tau and not pass_tau,
                # 'NoBaseTauTTJet'    : pass_ttjet and not pass_base_tau,
                # 'NoBaseMETNoTau' : not pass_trg and pass_met and not pass_tau,
                # 'NoBaseMETORTau': not pass_trg and (pass_met or pass_tau),
                # 'NoBaseMETTau': not pass_trg and pass_met and pass_tau,
                # 'NoBaseMETORTauORMu50': not pass_trg and (pass_met or pass_tau or pass_mu50),
                # 'NoBaseMETORTauOREle28': not pass_trg and (pass_met or pass_tau or pass_ele28),
                # 'NoBaseMETNoTTJet': not pass_trg and pass_met and not pass_ttjet,
                # 'NoBaseNoMETTTJet': not pass_trg and not pass_met and pass_ttjet,
                # 'NoBaseMETTTJet': not pass_trg and pass_met and pass_ttjet,
                # 'NoBaseTauTTJet': not pass_trg and pass_tau and pass_ttjet,
                # 'NoBaseTauNoTTJet': not pass_trg and pass_tau and not pass_ttjet, 
                # 'NoBaseNoTauTTJet': not pass_trg and not pass_tau and pass_ttjet,
                # 'NoBaseTauMu50': not pass_trg and pass_tau and pass_mu50,
                # 'NoBaseTauNoMu50': not pass_trg and pass_tau and not pass_mu50, 
                # 'NoBaseNoTauMu50': not pass_trg and not pass_tau and pass_mu50,
                # 'BaseMETTau'     : pass_trg and pass_met and pass_tau,
                # 'BaseNoMETNoTau' : pass_trg and not pass_met and not pass_tau,
                # 'VBFKin'         : False,
                # 'NoBaseMETORTauORTTJet' : not pass_trg and (pass_met or pass_tau or pass_ttjet)
                # 'NoBaseMu4JetsPNet': not pass_base_mu and pass_4jets_pNet,
                # 'NoBaseE4JetsPNet':  not pass_base_e and pass_4jets_pNet,
                # 'BaseMuOR4JetsPNet': pass_base_mu or pass_4jets_pNet,
                # 'BaseEOR4JetsPNet': pass_base_e or pass_4jets_pNet,
                # 'BaseTauOR4JetsPNet': pass_base_tau or pass_4jets_pNet,
                # 'NoBaseTau4JetsDeepJet': not pass_base_tau and pass_4jets_deepJet,
                # 'NoBaseTau4JetsPNet': not pass_base_tau and pass_4jets_pNet,
                # 'NoBaseTau4JetsDeepJetNo4JetsPNet': not pass_base_tau and pass_4jets_deepJet and not pass_4jets_pNet,
                # 'NoBaseTauNo4JetsDeepJet4JetsPNet': not pass_base_tau and not pass_4jets_deepJet and pass_4jets_pNet,
                # 'NoBaseTau4JetsDeepJet4JetsPNet': not pass_base_tau and pass_4jets_deepJet and pass_4jets_pNet,
                # 'NoBaseTauNoTTJet4JetsPNet': not pass_base_tau and not pass_ttjet and pass_4jets_pNet,
                # 'NoBaseTauTTJetNo4JetsPNet': not pass_base_tau and pass_ttjet and not pass_4jets_pNet,
                # 'NoBaseTauTTJet4JetsPNet': not pass_base_tau and pass_ttjet and pass_4jets_pNet,
                # 'NoBaseTauORMETORTauORTTJetOR4JetsPNet': not pass_base_tau and (pass_met or pass_tau or pass_ttjet or pass_4jets_pNet),
                # 'NoBaseMuNoMET4JetsPNet': not pass_base_mu and not pass_met and pass_4jets_pNet,
                # 'NoBaseMuMETNo4JetsPNet': not pass_base_mu and pass_met and not pass_4jets_pNet,
                # 'NoBaseMuMET4JetsPNet': not pass_base_mu and pass_met and pass_4jets_pNet,
                # 'NoBaseENoMET4JetsPNet': not pass_base_e and not pass_met and pass_4jets_pNet,
                # 'NoBaseEMETNo4JetsPNet': not pass_base_e and pass_met and not pass_4jets_pNet,
                # 'NoBaseEMET4JetsPNet': not pass_base_e and pass_met and pass_4jets_pNet,
                # 'NoBaseTauNoMET4JetsPNet': not pass_base_tau and not pass_met and pass_4jets_pNet,
                # 'NoBaseTauMETNo4JetsPNet': not pass_base_tau and pass_met and not pass_4jets_pNet,
                # 'NoBaseTauMET4JetsPNet': not pass_base_tau and pass_met and pass_4jets_pNet,
                # 'NoBaseMuTauNo4JetsPNet': not pass_base_mu and pass_tau and not pass_4jets_pNet,
                # 'NoBaseMuTau4JetsPNet': not pass_base_mu and pass_tau and pass_4jets_pNet,
                # 'NoBaseENoTau4JetsPNet': not pass_base_e and not pass_tau and pass_4jets_pNet,
                # 'NoBaseETauNo4JetsPNet': not pass_base_e and pass_tau and not pass_4jets_pNet,
                # 'NoBaseETau4JetsPNet': not pass_base_e and pass_tau and pass_4jets_pNet,
                'NoBaseTauNoTau4JetsPNet': not pass_base_tau and not pass_tau and pass_4jets_pNet,
                'NoBaseTauTauNo4JetsPNet': not pass_base_tau and pass_tau and not pass_4jets_pNet,
                'NoBaseTauTau4JetsPNet': not pass_base_tau and pass_tau and pass_4jets_pNet,
                # 'BaseEORMETORTauOR4JetsPNet': pass_base_e or pass_met or pass_tau or pass_4jets_pNet,
                # 'BaseMuORMETORTauORMu50OR4JetsPNet': pass_base_mu or pass_met or pass_tau or pass_mu50 or pass_4jets_pNet,
                # 'BaseTauORMETORTauORTTJetOR4JetsPNet': pass_base_tau or pass_met or pass_tau or pass_ttjet or pass_4jets_pNet,
                'BaseTauORTauTauJet': pass_base_tau or pass_ttjet,
                'BaseTauOR4JetsPNet': pass_base_tau or pass_4jets_pNet,
                # 'BaseTauOR4JetsPNetMatch': pass_base_tau or pass_4jets_match,
                'BaseTauORTauTauJetOR4JetsPNet': pass_base_tau or pass_ttjet or pass_4jets_pNet,
                'BaseTauORTauTauJetNo4JetsPNet': (pass_base_tau or pass_ttjet) and not pass_4jets_pNet,
                'BaseTauORTauTauJet4JetsPNet' : (pass_base_tau or pass_ttjet) and pass_4jets_pNet,
                'NoBaseTauTauTauJetNo4JetsPNet': not pass_base_tau and pass_ttjet and not pass_4jets_pNet,
                'NoBaseTauNoTauTauJet4JetsPNet': not pass_base_tau and not pass_ttjet and pass_4jets_pNet,
                'NoBaseTauORNoTauTauJet4JetsPNet': not (pass_base_tau or pass_ttjet) and pass_4jets_pNet,
                # 'BaseTauORMETORTauTauJetOR4JetsPNet': pass_base_tau or pass_met or pass_ttjet or pass_4jets_pNet,
                # 'BaseTauORTauORTauTauJetOR4JetsPNet': pass_base_tau or pass_tau or pass_ttjet or pass_4jets_pNet
        }
        # assert htypes == list(cuts.keys())
        
        w_mc     = entries.genWeight
        # w_pure   = entries.puWweight
        # # w_l1pref = entries.L1pref_weight
        # w_trig   = 1 # entries.trigSF
        # w_idiso  = entries.IdSF_deep_2d
        # w_jetpu  = entries.PUjetID_SF
        # w_btag   = entries.bTagweightReshape
        
        if utils.is_nan(w_mc)     : w_mc=1
        # if utils.is_nan(w_pure)   : w_pure=1
        # # if utils.is_nan(w_l1pref) : w_l1pref=1
        # if utils.is_nan(w_trig)   : w_trig=1
        # if utils.is_nan(w_idiso)  : w_idiso=1
        # if utils.is_nan(w_jetpu)  : w_jetpu=1
        # if utils.is_nan(w_btag)   : w_btag=1
        evt_weight = 1 # * w_btag
        # if evt_weight < 0.:
        #     print(w_mc, w_pure, w_trig, w_idiso, w_jetpu)
        tau_gen_cut = {"etau": None, "mutau": None,
                       "tautau": 'self.entries["isTau1real"] == 1 and self.entries["isTau2real"] == 1'}
        # if not pass_trg and entries.pairType == 0:
        #     print(channel, entries.pairType, utils.is_channel_consistent(channel, entries.pairType))
        # if utils.is_channel_consistent(channel, entries.pairType):
        #     if not sel.selection_cuts(lepton_veto=True, bjets_cut=True,
        #                               mass_cut=config_module.mass_cut,
        #                               custom_cut=tau_gen_cut[channel]):
        #         print("Not pass selection cut")
        #         continue

        for cat in categories:
            # if sel.sel_category(cat): # and entries.ditau_deltaR > deltaR:
            ntotal[cat] += 1

            if entries.pairType == 0:
                chn = 'mutau'
            elif entries.pairType == 1:
                chn = 'etau'
            elif entries.pairType == 2:
                chn = 'tautau'

            # if in_met:
            #     reg = 'met'
            # elif in_tau:
            #     reg = 'tau'
            # elif in_legacy:
            #     reg = 'legacy'
            # else:
            #     norphans[cat] += 1
            #     continue
            # assert reg in regions

            for key,cut in cuts.items():
                if cut:
                    # print(key, cat, chn)
                    for i in binning.keys():
                        z = 0
                        if 'bjet1' in i:
                            if 'eta' in i:
                                z = entries['Jet_eta'][entries['bjet1_JetIdx']]
                            if 'pt' in i:
                                z = entries['Jet_pt'][entries['bjet1_JetIdx']]
                            if 'pNet' in i:
                                z = entries['Jet_btagPNetB'][entries['bjet1_JetIdx']]
                        elif 'bjet2' in i:
                            if 'eta' in i:
                                z = entries['Jet_eta'][entries['bjet2_JetIdx']]
                            if 'pt' in i:
                                z = entries['Jet_pt'][entries['bjet2_JetIdx']]
                            if 'pNet' in i:
                                z = entries['Jet_btagPNetB'][entries['bjet2_JetIdx']]
                        else:
                            z = entries[i]
                        ahistos[key][cat][chn][i].fill(z, weight=evt_weight)

    # all MC and signal must be rescaled to get the correct number of events
    # for key,_ in cuts.items():
    #     for reg in regions:
    #         for cat in categories:
    #             for cc in ['mutau', 'etau', 'tautau']:
    #                 print("lumi", utils.get_lumi(args.year), utils.total_sum_weights(glob_files[0].replace("PreprocessRDF", "PreCounter").replace("/cat_base_selection", "").replace(".root", ".json"), isdata=False))
    #                 ahistos[key][reg][cat][cc] *= (utils.get_lumi(args.year) /
    #                                     utils.total_sum_weights(glob_files[0].replace("PreprocessRDF", "PreCounter").replace("/cat_base_selection", "").replace(".root", ".json"), isdata=False))
    
    # with open(outname, "wb") as f:
    #     pickle.dump(ahistos, f)
    file = uproot.recreate(outname)
    for key,_ in cuts.items():
        for cat in categories:
            for chn in ['mutau', 'etau', 'tautau']:
                for i in binning.keys():
                    file[str(key) + '_' + str(cat) + '_' + str(chn) + "_" + str(i)] = ahistos[key][cat][chn][i]

    for cat in categories:
        orph_frac = float(norphans[cat])/ntotal[cat]
        if orph_frac > 0.1:
            print('{}% orphans ({}/{}). This is unusual.'.format(orph_frac, norphans[cat], ntotal[cat]))
            print(met_region)
            print(tau_region)
            print(legacy_region)
        print('Category {} (m(X)=0GeV) had {} orphans ({}%)'.format(cat, norphans[cat], orph_frac))
    print('Raw histograms saved in {}.'.format(outname), flush=True)

if __name__ == '__main__':
    extensions = ('png',) #('png', 'pdf')
    triggers = {'etau': ("HLT_Ele30_WPTight_Gsf",
                    "HLT_Ele24_eta2p1_WPTight_Gsf_LooseDeepTauPFTauHPS30_eta2p1_CrossL1"),
                'mutau': ('HLT_IsoMu24', 'HLT_IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1'),
                'tautau': ('HLT_DoubleMediumDeepTauPFTauHPS35_L2NN_eta2p1',)
                }
    binning = {
        # 'metnomu_et': (20, 0, 450),
        #        'dau1_pt': (30, 0, 450),
        #        'dau1_eta': (20, -2.5, 2.5),
        #        'dau2_iso': (20, 0.88, 1.005),
        #        'dau2_pt': (30, 0, 400),
        #        'dau2_eta': (20, -2.5, 2.5),
        #        'ditau_deltaR': (30, 0.3, 1.3),
        #        'dib_deltaR': (25, 0, 2.5),
        #        'bH_pt': (20, 70, 600),
        #        'bH_mass': (30, 0, 280),
        #        'tauH_mass': (30, 0, 170),
        #        'tauH_pt': (30, 0, 500),
        #        'tauH_SVFIT_mass': (30, 0, 250),
        #        'tauH_SVFIT_pt': (20, 200, 650),
        #        'bjet1_pt': (25, 10, 600),
        #        'bjet2_pt': (25, 10, 550),
        #        'bjet1_eta': (20, -2.5, 2.5),
        #        'bjet2_eta': (20, -2.5, 2.5),
               }
    variables = tuple(binning.keys()) + ('HHKin_mass', 'dau1_iso')

    categories = ('baseline',) #('baseline', 's1b1jresolvedMcut', 's2b0jresolvedMcut', 'sboostedLLMcut')

    htypes = [
                'All', 'BaseMu', 'BaseE', 'BaseTau', 'MET', 'Tau', 'VBF', 'BaseMuMET', 'BaseEMET',
                'BaseTauMET', 'BaseMuORMET', 'BaseEORMET', 'BaseTauORMET', 'BaseMuTau', 'BaseETau', 'BaseTauTau', 'BaseMuORTau',
                'BaseEORTau', 'BaseTauORTau', 'BaseTauORTTJet', 'BaseTauORTauORTTJet',
                'METTau', 'METORTau', 'BaseMuORMETORTau', 'BaseEORMETORTau' , 'BaseTauORMETORTau', 'BaseTauORMETORTauORTTJet', 'NoBaseMuMET',
                'BaseMuNoMET', 'NoBaseEMET', 'BaseENoMET', 'NoBaseTauMET', 'BaseTauNoMET', 'NoBaseMuNoMETTau', 'NoBaseENoMETTau', 'NoBaseTauNoMETTau', 'NoBaseTauTau', 'NoBaseMuTau', 'NoBaseETau',
                # 'NoBaseMuNoMETMu50', 'NoBaseMuNoTauMu50', 'NoBaseMuMETNoMu50', 'NoBaseMuMETMu50', 
                'BaseTauNoTau', 'NoBaseTauTTJet', 
                'NoBaseMu4JetsPNet', 'NoBaseE4JetsPNet', 'BaseMuOR4JetsPNet', 'BaseEOR4JetsPNet', 'BaseTauOR4JetsPNet', 'NoBaseTau4JetsDeepJet', 'NoBaseTau4JetsPNet', 'NoBaseTau4JetsDeepJetNo4JetsPNet',
                'NoBaseTauNo4JetsDeepJet4JetsPNet', 'NoBaseTau4JetsDeepJet4JetsPNet',
                'NoBaseTauNoTTJet4JetsPNet', 'NoBaseTauTTJetNo4JetsPNet', 'NoBaseTauTTJet4JetsPNet', 'NoBaseTauORMETORTauORTTJetOR4JetsPNet',
                'NoBaseMuNoMET4JetsPNet', 'NoBaseMuMETNo4JetsPNet', 'NoBaseMuMET4JetsPNet', 'NoBaseENoMET4JetsPNet', 'NoBaseEMETNo4JetsPNet', 'NoBaseEMET4JetsPNet', 'NoBaseTauNoMET4JetsPNet' , 'NoBaseTauMETNo4JetsPNet' ,
                'NoBaseTauMET4JetsPNet', 'NoBaseMuTauNo4JetsPNet', 'NoBaseMuTau4JetsPNet', 'NoBaseENoTau4JetsPNet', 'NoBaseETauNo4JetsPNet', 'NoBaseETau4JetsPNet', 'NoBaseTauNoTau4JetsPNet', 'NoBaseTauTauNo4JetsPNet', 'NoBaseTauTau4JetsPNet',
                'BaseEORMETORTauOR4JetsPNet', 'BaseTauORMETORTauORTTJetOR4JetsPNet',
                'BaseTauORTauTauJet', 'BaseTauOR4JetsPNet', 'BaseTauORTauTauJetOR4JetsPNet', 'BaseTauORTauTauJetNo4JetsPNet', 'BaseTauORTauTauJet4JetsPNet', 'NoBaseTauTauTauJetNo4JetsPNet',
                'NoBaseTauNoTauTauJet4JetsPNet', 'BaseTauORMETORTauTauJetOR4JetsPNet', 'BaseTauORTauORTauTauJetOR4JetsPNet', 'NoBaseTauORNoTauTauJet4JetsPNet']    
    
    # Parse input arguments
    desc = 'Producer trigger histograms.\n'
    desc += "Run example: python tests/test_trigger_regions.py --indir /data_CMS/cms/alves/HHresonant_SKIMS/SKIMS_UL18_EOSv4_Signal/ --masses 400 500 600 700 800 900 1000 1250 1500 --channels ETau --met_turnon 180 --region_cuts 40 40 --copy"
    parser = argparse.ArgumentParser(description=desc, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('--indir', required=True, type=str,
                        help='Full path of ROOT input file')
    # parser.add_argument('--masses', required=True, nargs='+', type=str,
    #                     help='Resonance mass')
    parser.add_argument('--channel', required=True, type=str,  
                        help='Select the channel over which the workflow will be run.' )
    parser.add_argument('--year', required=True, type=str, choices=('2016', '2016APV', '2017', '2018', '2022'),
                        help='Select the year over which the workflow will be run.' )
    # parser.add_argument('--spin', required=True, type=int, choices=(0, 2), 
    #                     help='Select the spin hypothesis over which the workflow will be run.' )
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
        
    main_dir = os.path.join(os.path.join('/t3home/', os.environ['USER'], 'TriggerScaleFactors'),
                            '_'.join((args.channel, *[str(x) for x in regcuts],
                                      'DR', str(args.deltaR), 'PT', *[str(x) for x in ptcuts], 'TURNON',
                                      str(met_turnon))))
    if args.bigtau:
        main_dir += '_BIGTAU'

    regions = ('legacy', 'met', 'tau')
        
    #### run main function ###
    if not args.plot:
        # for sample in args.masses:
        trigger_regions(args.indir, args.channel, args.year, args.deltaR)
        # if args.sequential:
        #     pass
        # else:
        #     pool = multiprocessing.Pool(processes=6)
        #     pool.starmap(trigger_regions,
        #                  zip(it.repeat(args.indir), it.repeat(args.channel), it.repeat(args.year), it.repeat(args.deltaR)))

    ###########################

    sum_stats, err_sum_stats = ([] for _ in range(2))
    contam1, contam2, contam1_errors, contam2_errors = ([] for _ in range(4))
    from_directory = os.path.join(main_dir, args.channel)
    outname = get_outname(args.channel, args.bigtau)
    with open(outname, "rb") as f:
        ahistos = pickle.load(f)

    # write csv header, one per category
    out_counts = []
    for cat in categories:
        out_counts.append( os.path.join(from_directory, cat, 'counts') )
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

    text = {'out': os.path.join(out_counts[categories.index(cat)], 'diagram.html')}
    sq_res = square_diagram(c_legacy_trg, c_met_trg, c_tau_trg, args.channel,
                            [str(x) for x in ptcuts], text=text, notext=args.notext,
                            bigtau=args.bigtau)
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
    # masses = [float(x) for x in args.masses]
    contamination_save('data', '_'.join([str(x) for x in regcuts]) + '_' + args.channel,
                       contam1, contam2, contam1_errors, contam2_errors)
    stats_save('data', '_'.join([str(x) for x in regcuts]) + '_' + args.channel,
               sum_stats, err_sum_stats, mode='a')

    if args.copy:
        import subprocess
        to_directory = os.path.join('/eos/home-b/bfontana/www/TriggerScaleFactors', main_dir)
        to_directory = os.path.join(to_directory, args.channel)
     
        for sample in args.masses:
            sample_from = os.path.join(from_directory, sample)
            print('Copying: {}\t\t--->\t{}'.format(sample_from, to_directory), flush=True)
            subprocess.run(['rsync', '-ah', sample_from, to_directory])
