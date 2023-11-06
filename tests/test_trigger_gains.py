# coding: utf-8

_all_ = [ 'test_trigger_gains' ]

import os
import sys
parent_dir = os.path.abspath(__file__ + 2 * '/..')
sys.path.insert(0, parent_dir)

from pathlib import Path
import csv
import argparse
from inclusion.utils import utils
import numpy as np
from collections import defaultdict as dd
 
import bokeh
from bokeh.plotting import figure, output_file, save
output_file('plots/trigger_gains.html')
from bokeh.models import Whisker
from bokeh.layouts import gridplot
#from bokeh.io import export_svg

tau = '\u03C4'
mu  = '\u03BC'
pm  = '\u00B1'
ditau = tau+tau

def pp(chn):
    if chn == "tautau":
        return ditau
    elif chn == "etau":
        return "e" + tau
    elif chn == "mutau":
        return mu + tau
    
def ZeroDivError(func):
    try:
        res = func()
    except ZeroDivisionError:
        res = 0.
    return res

def error(v1, v2):
    return ZeroDivError(lambda : np.sqrt(1/v1 + 1/v2))

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


def main(args):
    channels = args.channels
    linear_x = [k for k in range(1,len(args.masses)+1)]
     
    yone, yboth, ykin = (dd(lambda: dd(dict)) for _ in range(3))
    eone, eboth, ekin = (dd(lambda: dd(dict)) for _ in range(3))
    for md in main_dir:
        d_base = Path(base_dir) / md
        output_file(d_base / 'trigger_gains.html')
        print('Saving file {}.'.format(d_base / 'trigger_gains.html'))
        for chn in channels:
            yone[md][chn]['met'],  yone[md][chn]['tau'], yone[md][chn]['vbf'] = [], [], []
            yboth[md][chn]['met'], yboth[md][chn]['two'] = [], []
            ykin[md][chn]['met'],  ykin[md][chn]['tau'] = [], []
            ykin[md][chn]['two'], ykin[md][chn]['vbf'] = [], []
            
            eone[md][chn]['met'],  eone[md][chn]['tau'], eone[md][chn]['vbf']  = [], [], []
            eboth[md][chn]['met'], eboth[md][chn]['two'] = [], []
            ekin[md][chn]['met'], ekin[md][chn]['tau'] = [], []
            ekin[md][chn]['two'], ekin[md][chn]['vbf'] = [], []
  
            for mass in args.masses:
                d = d_base / chn / str(mass)
                fullpath = Path(d) / 'baseline/counts/table.csv'
  
                with open(fullpath) as f:
                    reader = csv.reader(f, delimiter=',', quotechar='|')
                    next(reader, None) #ignore header line

                    line = next(reader, None)
                    assert line[0] == "ditau"
                    sum_base     = float(line[1])
                    sum_vbf      = float(line[4])
                    sum_met      = float(line[8])
                    sum_only_tau = float(line[9])
                    sum_tau      = float(line[10])

                    line = next(reader, None)
                    assert line[0] == "met"
                    sum_metkin   = float(line[15])
                    sum_vbfkin   = float(line[17])

                    line = next(reader, None)
                    assert line[0] == "tau"
                    sum_taukin   = float(line[16])
                        
                    frac_one_met = sum_met / sum_base
                    frac_one_tau = sum_tau / sum_base
                    frac_both    = (sum_met + sum_only_tau) / sum_base
                    frac_one_vbf     = sum_vbf / sum_base

                    err_one_met = error(sum_base, sum_met)
                    err_one_tau = error(sum_base, sum_tau)
                    err_one_vbf = error(sum_base, sum_vbf)
                    err_both    = error(sum_base, sum_met+sum_only_tau)
                    err_vbf     = error(sum_base, sum_vbf)
                    
                    frac_metkin  = sum_metkin / sum_base
                    frac_taukin  = sum_taukin / sum_base
                    frac_bothkin = frac_metkin + frac_taukin
                    frac_vbfkin  = sum_vbfkin / sum_base

                    err_metkin  = error(sum_base, sum_metkin)
                    err_taukin  = error(sum_base, sum_taukin)
                    err_bothkin = error(sum_base, sum_metkin+sum_taukin)
                    err_vbfkin  = error(sum_base, sum_vbfkin)

                    yone[md][chn]['met'].append(frac_one_met*100)
                    yone[md][chn]['tau'].append(frac_one_tau*100)
                    yone[md][chn]['vbf'].append(frac_one_vbf*100)
                    yboth[md][chn]['met'].append(frac_one_met*100)
                    yboth[md][chn]['two'].append(frac_both*100)
                    ykin[md][chn]['met'].append(frac_metkin*100)
                    ykin[md][chn]['tau'].append(frac_taukin*100)
                    ykin[md][chn]['two'].append(frac_bothkin*100)
                    ykin[md][chn]['vbf'].append(frac_vbfkin*100)

                    eone[md][chn]['met'].append(err_one_met*100)
                    eone[md][chn]['tau'].append(err_one_tau*100)
                    eone[md][chn]['vbf'].append(err_one_vbf*100)
                    eboth[md][chn]['met'].append(err_one_met*100)
                    eboth[md][chn]['two'].append(err_both*100)
                    ekin[md][chn]['met'].append(err_metkin*100)
                    ekin[md][chn]['tau'].append(err_taukin*100)
                    ekin[md][chn]['two'].append(err_bothkin*100)
                    ekin[md][chn]['vbf'].append(err_vbfkin*100)

    opt_points = dict(size=8)
    opt_line = dict(width=1.5)
    colors = ('green', 'blue', 'red', 'brown')
    styles = ('solid', 'dashed', 'dotdash')
    legends = {'met': ' (MET)', 'tau': ' (Tau)',
               'two': ' (MET + Tau)', 'vbf': ' (VBF)'}
     
    x_str = [str(k) for k in args.masses]
    xticks = linear_x[:]
    yticks = [x for x in range(0,110,10)]
    shift_one = {'met': [-0.15, 0., 0.15],  'tau': [-0.20, -0.05, 0.1],
                 'vbf': [-0.10, 0.05, 0.20]}
    shift_both = {'met': [-0.15, 0., 0.15], 'two': [-0.20, -0.05, 0.1]}
    shift_kin = {'met': [-0.15, 0., 0.15],  'tau': [-0.20, -0.05, 0.1],
                 'two': [-0.10, 0.05, 0.20], 'vbf': [-0.05, 0.1, 0.25]}
     
    for md in main_dir:
        p_opt = dict(width=800, height=400, x_axis_label='x', y_axis_label='y')
        p1 = figure(title='Inclusion', **p_opt)
        p2 = figure(title='Acceptance gain', **p_opt)
        p3 = figure(title='Acceptance gain in bb' + tau + tau + ' kin region', **p_opt)
        for p in (p1, p2, p3):
            set_fig(p)
        for ichn,chn in enumerate(channels):

            for itd,td in enumerate(('met', 'tau', 'vbf')):
                p1.circle([x+shift_one[td][ichn] for x in linear_x],
                          yone[md][chn][td], color=colors[itd], fill_alpha=1.,
                          **opt_points)
                p1.line([x+shift_one[td][ichn] for x in linear_x],
                        yone[md][chn][td], color=colors[itd], line_dash=styles[ichn],
                        legend_label=pp(chn)+legends[td], **opt_line)
                p1.multi_line(
                    [(x+shift_one[td][ichn],x+shift_one[td][ichn]) for x in linear_x],
                    [(max(0,x-y/2),min(100,x+y/2))
                     for x,y in zip(yone[md][chn][td],eone[md][chn][td])],
                    color=colors[itd], **opt_line)

            for itd,td in enumerate(('met', 'two')):
                p2.circle([x+shift_both[td][ichn] for x in linear_x],
                          yboth[md][chn][td], color=colors[itd], fill_alpha=1.,
                          **opt_points)
                p2.line([x+shift_both[td][ichn] for x in linear_x],
                        yboth[md][chn][td], color=colors[itd], line_dash=styles[ichn],
                        legend_label=pp(chn)+legends[td], **opt_line)
                p2.multi_line(
                    [(x+shift_both[td][ichn],x+shift_both[td][ichn]) for x in linear_x], 
                    [(max(0,x-y/2),min(100,x+y/2))
                     for x,y in zip(yboth[md][chn][td],eboth[md][chn][td])],
                    color=colors[itd], **opt_line)

            for itd,td in enumerate(('met', 'tau', 'two', 'vbf')):
                p3.circle([x+shift_kin[td][ichn] for x in linear_x],
                          ykin[md][chn][td], color=colors[itd], fill_alpha=1.,
                          **opt_points)
                p3.line([x+shift_kin[td][ichn] for x in linear_x],
                        ykin[md][chn][td], color=colors[itd], line_dash=styles[ichn],
                        legend_label=pp(chn)+legends[td], **opt_line)
                p3.multi_line(
                    [(x+shift_kin[td][ichn],x+shift_kin[td][ichn]) for x in linear_x], 
                    [(max(0,x-y/2),min(100,x+y/2))
                     for x,y in zip(ykin[md][chn][td],ekin[md][chn][td])],
                    color=colors[itd], **opt_line)

        for p in (p1, p2, p3):
            p.xaxis[0].ticker = xticks
            p.xgrid[0].ticker = xticks
            p.xgrid.grid_line_alpha = 0.2
            p.xgrid.grid_line_color = 'black'
            p.yaxis[0].ticker = yticks
            p.ygrid[0].ticker = yticks
            p.ygrid.grid_line_alpha = 0.2
            p.ygrid.grid_line_color = 'black'
             
            p.legend.location = 'top_left'
            p.xaxis.axis_label = "m(X) [GeV]"
            p.yaxis.axis_label = 'Trigger acceptance gain (w.r.t. trigger baseline) [%]'
         
            p.xaxis.major_label_overrides = dict(zip(linear_x,x_str))
     
            p.legend.click_policy='hide'        
     
            p.output_backend = 'svg'
            #export_svg(p, filename='line_graph.svg')
         
        g = gridplot([[p1], [p2], [p3]])
        save(g, title=md)

if __name__ == '__main__':
    base_dir = '/eos/home-b/bfontana/www/TriggerScaleFactors/'
    main_dir = ['Region_Spin0_190_190_PT_40_40_TURNON_200_190',]
        #'Region_1000_1000_PT_40_40_TURNON_200_190',]
    #'TriggerStudy_MET200_SingleTau190_CUT_entries_ditau_deltaR_GT_0_5',
    #'TriggerStudy_MET200_SingleTau190_CUT_entries_ditau_deltaR_ST_0_5']
    
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
    args = utils.parse_args(parser)
    main(args)
