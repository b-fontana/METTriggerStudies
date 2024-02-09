# coding: utf-8

_all_ = [ 'test_trigger_gains' ]

import os
import sys
parent_dir = os.path.abspath(__file__ + 2 * '/..')
sys.path.insert(0, parent_dir)

import csv
import argparse
from inclusion.utils import utils
import numpy as np
from collections import defaultdict as dd
import hist
from hist.intervals import clopper_pearson_interval as clop
import pickle

import bokeh
from bokeh.plotting import figure, output_file, save
from bokeh.models import Whisker
from bokeh.layouts import gridplot
#from bokeh.io import export_svg

tau = '\u03C4'
mu  = '\u03BC'
pm  = '\u00B1'
ditau = tau+tau

def get_outname(sample, channel, regcuts, ptcuts, met_turnon, tau_turnon,
                bigtau, notau, nomet):
    utils.create_single_dir('data')

    name = sample + '_' + channel + '_'
    name += '_'.join((*regcuts, 'ptcuts', *ptcuts, 'turnon', met_turnon, tau_turnon))
    if bigtau:
        name += '_BIGTAU'
    if notau:
        name += '_NOTAU'
    if nomet:
        name += '_NOMET'
    name += '.pkl'

    s = 'data/regions_{}'.format(name)
    return s

def pp(chn):
    if chn == "tautau":
        return ditau
    elif chn == "etau":
        return "e" + tau
    elif chn == "mutau":
        return mu + tau

def rec_dd():
    return dd(rec_dd)

def set_fig(fig, legend=True):
    fig.output_backend = 'svg'
    fig.toolbar.logo = None
    # if legend:
    #     fig.legend.click_policy='hide'
    #     fig.legend.location = 'top_left'
    #     fig.legend.label_text_font_size = '8pt'
    fig.min_border_bottom = 5
    fig.xaxis.visible = True
    fig.title.align = "left"
    fig.title.text_font_size = "15px"
    fig.xaxis.axis_label_text_font_style = "bold"
    fig.yaxis.axis_label_text_font_style = "bold"
    fig.xaxis.axis_label_text_font_size = "13px"
    fig.yaxis.axis_label_text_font_size = "13px"

def main(args):
    channels = args.channels
    linear_x = [k for k in range(1,len(args.masses)+1)]
    edges_x  = [k-0.5 for k in range(1,len(args.masses)+1)] + [len(args.masses)+0.5]
    ptcuts = {chn: utils.get_ptcuts(chn, args.year) for chn in args.channels}

    nevents = dd(lambda: dd(dict))
    yone, yboth, ykin = (dd(lambda: dd(dict)) for _ in range(3))
    eone, eboth, ekin = (dd(lambda: dd(dict)) for _ in range(3))
    for adir in main_dir:
        dRstr = str(args.deltaR).replace('.', 'p')
        if len(args.channels) == 1:
            output_name = os.path.join(base_dir, 'trigger_gains_{}_{}_DR{}'.format(args.channels[0],
                                                                                        args.year, dRstr))
        elif len(args.channels) == 2:
            output_name = os.path.join(base_dir, 'trigger_gains_{}_{}_{}_DR{}'.format(*args.channels[:2],
                                                                                           args.year, dRstr))
        elif len(args.channels) == 3:
            output_name = os.path.join(base_dir, 'trigger_gains_all_{}_DR{}'.format(args.year, dRstr))
        if args.bigtau:
            output_name += "_BIGTAU"
        output_name += ".html"
        output_file(output_name)
        print('Saving file {}.'.format(output_name))

        for chn in channels:
            md = adir[chn]
            in_base = os.path.join(base_dir, md)

            nevents[md][chn]['base'], nevents[md][chn]['met'],  nevents[md][chn]['tau'] = [], [], []
            
            yone[md][chn]['met'],  yone[md][chn]['tau'], yone[md][chn]['vbf'] = [], [], []
            yboth[md][chn]['met'], yboth[md][chn]['two'] = [], []
            ykin[md][chn]['met'],  ykin[md][chn]['tau']  = [], []
            ykin[md][chn]['two'],  ykin[md][chn]['vbf']  = [], []
            
            eone[md][chn]['met'],  eone[md][chn]['tau'], eone[md][chn]['vbf']  = [], [], []
            eboth[md][chn]['met'], eboth[md][chn]['two'] = [], []
            ekin[md][chn]['met'],  ekin[md][chn]['tau']  = [], []
            ekin[md][chn]['two'],  ekin[md][chn]['vbf']  = [], []
  
            for mass in args.masses:
                outname = get_outname(mass, chn, args.region_cuts, ptcuts[chn],
                                      args.met_turnon, args.tau_turnon,
                                      args.bigtau, args.notau, args.nomet)

                with open(outname, "rb") as f:
                    ahistos = pickle.load(f)

                    # legacy region
                    l1 = lambda x : round(x["legacy"]["baseline"].values().sum(), 2)
                    sum_base     = l1(ahistos["Base"])
                    sum_vbf      = l1(ahistos["VBF"])
                    sum_met      = l1(ahistos["NoBaseMET"])
                    sum_only_tau = l1(ahistos["NoBaseNoMETTau"])
                    sum_tau      = l1(ahistos["NoBaseTau"])
                    sum_basekin  = l1(ahistos["LegacyKin"])

                    # MET region
                    l2 = lambda x : round(x["met"]["baseline"].values().sum(), 2)
                    sum_metkin   = l2(ahistos["METKin"])

                    # Single Tau region
                    l3 = lambda x : round(x["tau"]["baseline"].values().sum(), 2)
                    sum_taukin   = l3(ahistos["TauKin"])

                    # hypothetical VBF region
                    sum_vbfkin   = l2(ahistos["VBFKin"]) + l3(ahistos["VBFKin"])

                    valerr = lambda num, den : ((num/den)*100., clop(num,den)*100.)
                    frac_one_met, err_one_met = valerr(sum_met, sum_base)
                    frac_one_tau, err_one_tau = valerr(sum_tau, sum_base)
                    frac_both, err_both       = valerr(sum_met + sum_only_tau, sum_base)
                    frac_one_vbf, err_one_vbf = valerr(sum_vbf, sum_base)
                    frac_metkin,  err_metkin  = valerr(sum_metkin, sum_base)
                    frac_taukin,  err_taukin  = valerr(sum_taukin, sum_base)
                    frac_bothkin, err_bothkin = valerr(sum_metkin + sum_taukin, sum_base)
                    frac_vbfkin,  err_vbfkin  = valerr(sum_vbfkin, sum_base)

                    nevents[md][chn]['base'].append(sum_basekin)
                    nevents[md][chn]['met'].append(sum_basekin + sum_metkin)
                    nevents[md][chn]['tau'].append(sum_basekin + sum_metkin + sum_taukin)
                    
                    yone[md][chn]['met'].append(frac_one_met)
                    yone[md][chn]['tau'].append(frac_one_tau)
                    yone[md][chn]['vbf'].append(frac_one_vbf)
                    yboth[md][chn]['met'].append(frac_one_met)
                    yboth[md][chn]['two'].append(frac_both)
                    ykin[md][chn]['met'].append(frac_metkin)
                    ykin[md][chn]['tau'].append(frac_taukin)
                    ykin[md][chn]['two'].append(frac_bothkin)
                    ykin[md][chn]['vbf'].append(frac_vbfkin)

                    eone[md][chn]['met'].append(err_one_met)
                    eone[md][chn]['tau'].append(err_one_tau)
                    eone[md][chn]['vbf'].append(err_one_vbf)
                    eboth[md][chn]['met'].append(err_one_met)
                    eboth[md][chn]['two'].append(err_both)
                    ekin[md][chn]['met'].append(err_metkin)
                    ekin[md][chn]['tau'].append(err_taukin)
                    ekin[md][chn]['two'].append(err_bothkin)
                    ekin[md][chn]['vbf'].append(err_vbfkin)

    opt_points = dict(size=8)
    opt_line = dict(width=1.5)
    colors = ('green', 'blue', 'red', 'brown')
    styles = ('solid', 'dashed', 'dotdash')
    legends = {'base': 'Legacy',
               'met': 'MET', 'tau': 'Single Tau',
               'two': 'MET + Single Tau', 'vbf': 'VBF'}
     
    x_str = [str(k) for k in args.masses]
    xticks = linear_x[:]
    yticks = [x for x in range(0,110,5)]
    shift_one = {'met': [-0.15, 0., 0.15],  'tau': [-0.20, -0.05, 0.1],
                 'vbf': [-0.10, 0.05, 0.20]}
    shift_both = {'met': [-0.15, 0., 0.15], 'two': [-0.20, -0.05, 0.1]}
    shift_kin = {'met': [-0.09, 0., 0.15],  'tau': [0.03, -0.05, 0.1],
                 'two': [-0.03, 0.05, 0.20], 'vbf': [0.09, 0.1, 0.25]}
     
    for adir in main_dir:
        p_opt = dict(width=800, height=400, x_axis_label='x', y_axis_label='y')
        p1 = figure(title='Event number (' + pp(channels[0]) + ')', y_axis_type="linear", **p_opt)
        p2 = figure(title='Acceptance Gain (' + pp(channels[0]) + ')', **p_opt) if len(channels)==1 else figure(**p_opt)
        p1.yaxis.axis_label = 'Weighted number of events'
        p2.yaxis.axis_label = 'Trigger acceptance gain (w.r.t. trigger baseline) [%]'
        pics = (p1, p2)
        for p in pics:
            set_fig(p)

        for ichn,chn in enumerate(channels):
            md = adir[chn]
            
            p1.quad(top=nevents[md][chn]["base"], bottom=0,
                    left=edges_x[:-1], right=edges_x[1:],
                    legend_label=legends["base"]+(' ('+pp(chn)+')' if len(channels)>1 else ''),
                    fill_color="dodgerblue", line_color="black")
            p1.quad(top=nevents[md][chn]["met"], bottom=nevents[md][chn]["base"],
                    left=edges_x[:-1], right=edges_x[1:],
                    legend_label=legends["met"]+(' ('+pp(chn)+')' if len(channels)>1 else ''),
                    fill_color="green", line_color="black")
            p1.quad(top=nevents[md][chn]["tau"], bottom=nevents[md][chn]["met"],
                    left=edges_x[:-1], right=edges_x[1:],
                    legend_label=legends["tau"]+(' ('+pp(chn)+')' if len(channels)>1 else ''),
                    fill_color="red", line_color="black")
            # p1.multi_line([(x,x) for x in linear_x],
            #               [(-0.1,0.1) for x in linear_x],
            #               color=colors[itd], **opt_line)

            for itd,td in enumerate(('met', 'tau', 'two')):
                p2.circle([x+shift_kin[td][ichn] for x in linear_x],
                          ykin[md][chn][td], color=colors[itd], fill_alpha=1.,
                          **opt_points)
                p2.line([x+shift_kin[td][ichn] for x in linear_x],
                        ykin[md][chn][td], color=colors[itd], line_dash=styles[ichn],
                        legend_label=legends[td]+(' ('+pp(chn)+')' if len(channels)>1 else ''), **opt_line)
                p2.multi_line(
                    [(x+shift_kin[td][ichn],x+shift_kin[td][ichn]) for x in linear_x], 
                    [(y[0],y[1])
                     for x,y in zip(ykin[md][chn][td],ekin[md][chn][td])],
                    color=colors[itd], **opt_line)

        p1.legend.location = 'top_right'
        p2.legend.location = 'top_left'                
        for p in pics:
            p.xaxis[0].ticker = xticks
            p.xgrid[0].ticker = xticks
            p.xgrid.grid_line_alpha = 0.2
            p.xgrid.grid_line_color = 'black'
            # p.yaxis[0].ticker = yticks
            # p.ygrid[0].ticker = yticks
            p.ygrid.grid_line_alpha = 0.2
            p.ygrid.grid_line_color = 'black'
             
            p.xaxis.axis_label = "m(X) [GeV]"
         
            p.xaxis.major_label_overrides = dict(zip(linear_x,x_str))
     
            p.legend.click_policy='hide'        
     
            p.output_backend = 'svg'
            #export_svg(p, filename='line_graph.svg')
         
        g = gridplot([[p] for p in pics])
        save(g, title=md)

if __name__ == '__main__':
    desc = "Produce plots of trigger gain VS resonance mass.\n"
    desc += "Uses the output of test_trigger_regions.py."
    desc += "When running on many channels, one should keep in mind each channel has different pT cuts."
    desc += "This might imply moving sub-folders (produced by the previous script) around."
    parser = argparse.ArgumentParser(description=desc, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('--masses', required=True, nargs='+', type=str,
                        help='Resonance mass')
    parser.add_argument('--channels', required=True, nargs='+', type=str, 
                        choices=('etau', 'mutau', 'tautau'),
                        help='Select the channel over which the workflow will be run.' )
    parser.add_argument('--year', required=True, type=str, choices=('2016', '2017', '2018'),
                        help='Select the year over which the workflow will be run.' )
    parser.add_argument('--deltaR', type=float, default=0.5, help='DeltaR between the two leptons.')
    parser.add_argument('--bigtau', action='store_true',
                        help='Consider a larger single tau region, reducing the ditau one.')
    parser.add_argument('--notau', action='store_true',
                        help='Remove the single tau region (default analysis).')
    parser.add_argument('--nomet', action='store_true',
                        help='Remove the MET region (default analysis).')
    parser.add_argument('--met_turnon', required=False, type=str,  default='180',
                        help='MET trigger turnon cut [GeV].' )
    parser.add_argument('--tau_turnon', required=False, type=str,  default='190',
                        help='MET trigger turnon cut [GeV].' )
    parser.add_argument('--region_cuts', required=False, type=str, nargs=2, default=('200', '200'),
                        help='High/low regions pT1 and pT2 selection cuts [GeV].' )

    args = utils.parse_args(parser)

    base_dir = '/eos/home-b/bfontana/www/TriggerScaleFactors/'
    main_dir = [{"etau":   "Region_Spin0_190_190_PT_33_25_35_DR_{}_TURNON_200_190".format(args.deltaR),
                 "mutau":  "Region_Spin0_190_190_PT_25_21_32_DR_{}_TURNON_200_190".format(args.deltaR),
                 "tautau": "Region_Spin0_190_190_PT_40_40_DR_{}_TURNON_200_190".format(args.deltaR)},
                ]
    
    main(args)
