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

import bokeh
from bokeh.plotting import figure, output_file, save
output_file('plots/trigger_gains.html')
from bokeh.models import Whisker
from bokeh.layouts import gridplot
#from bokeh.io import export_svg

def ZeroDivError(func):
    try:
        res = func()
    except ZeroDivisionError:
        res = 0.
    return res


def main(args):
    channels = args.channels
    linear_x = [k for k in range(1,len(args.masses)+1)]
     
    yindep, yboth = ({} for _ in range(2))
    errindep, errboth = ({} for _ in range(2))
    for md in main_dir:
        d_base = Path(base_dir) / md
        output_file(d_base / 'trigger_gains.html')
        print('Saving file {}.'.format(d_base / 'trigger_gains.html'))
        yindep[md], yboth[md] = ({} for _ in range(2))
        errindep[md], errboth[md] = ({} for _ in range(2))
        for chn in channels:
            yindep[md][chn], yboth[md][chn] = ({} for _ in range(2))
            errindep[md][chn], errboth[md][chn] = ({} for _ in range(2))
            yindep[md][chn]['met'], yindep[md][chn]['tau'] = ([] for _ in range(2))
            errindep[md][chn]['met'], errindep[md][chn]['tau'] = ([] for _ in range(2))
            yboth[md][chn]['met'],  yboth[md][chn]['add_met_tau']  = ([] for _ in range(2))
            errboth[md][chn]['met'],  errboth[md][chn]['add_met_tau']  = ([] for _ in range(2))
  
            for mass in args.masses:
                d = d_base / chn / str(mass)
                fullpath = Path(d) / 'baseline/counts/table.csv'
  
                with open(fullpath) as f:
                    reader = csv.reader(f, delimiter=',', quotechar='|')
                    next(reader, None) #ignore header line
                    ditau_line = next(reader, None)
  
                    sum_base     = float(ditau_line[1])
                    sum_met      = float(ditau_line[7])
                    sum_tau      = float(ditau_line[9])
                    sum_only_tau = float(ditau_line[8])
                        
                    frac_indep_met = sum_met / sum_base
                    err_indep_met = ZeroDivError(lambda : np.sqrt(1/sum_base + 1/sum_met))
                    frac_indep_tau = sum_tau / sum_base
                    err_indep_tau = ZeroDivError(lambda : np.sqrt(1/sum_base + 1/sum_tau))
                    frac_both  = (sum_met + sum_only_tau) / sum_base
                    err_rel_both = sum_met + sum_only_tau
                    err_both = ZeroDivError(lambda : np.sqrt(1/sum_base + 1/err_rel_both))

                    yindep[md][chn]['met'].append(frac_indep_met*100)
                    yindep[md][chn]['tau'].append(frac_indep_tau*100)
                    yboth[md][chn]['met'].append(frac_indep_met*100)
                    yboth[md][chn]['add_met_tau'].append(frac_both*100)
                    errindep[md][chn]['met'].append(err_indep_met*100)
                    errindep[md][chn]['tau'].append(err_indep_tau*100)
                    errboth[md][chn]['met'].append(err_indep_met*100)
                    errboth[md][chn]['add_met_tau'].append(err_both*100)

    opt_points = dict(size=8)
    opt_line = dict(width=1.5)
    colors = ('green', 'blue', 'red')
    styles = ('solid', 'dashed')
    legends = {'met': ' (MET)', 'tau': ' (Tau)', 'add_met_tau': ' (MET + Tau)'}
     
    x_str = [str(k) for k in args.masses]
    xticks = linear_x[:]
    yticks = [x for x in range(0,110,10)]
    shift_indep = {'met': [-0.15, 0., 0.15], 'tau': [-0.20, -0.05, 0.1]}
    shift_both = {'met': [-0.15, 0., 0.15], 'add_met_tau': [-0.20, -0.05, 0.1]}
     
    for md in main_dir:
        p_opt = dict(width=800, height=400, x_axis_label='x', y_axis_label='y')
        p1 = figure(title='Inclusion of MET or Single Tau triggers', **p_opt)
        p2 = figure(title='Acceptance gain of MET + SingleTau triggers', **p_opt)
        p1.toolbar.logo = None
        p2.toolbar.logo = None
        for ichn,chn in enumerate(channels):
            for itd,td in enumerate(('met', 'tau')):
                p1.circle([x+shift_indep[td][ichn] for x in linear_x],
                          yindep[md][chn][td], color=colors[ichn], fill_alpha=1.,
                          **opt_points)
                p1.line([x+shift_indep[td][ichn] for x in linear_x],
                        yindep[md][chn][td], color=colors[ichn], line_dash=styles[itd],
                        legend_label=chn+legends[td], **opt_line)
                p1.multi_line(
                    [(x+shift_indep[td][ichn],x+shift_indep[td][ichn]) for x in linear_x],
                    [(max(0,x-y/2),min(100,x+y/2))
                     for x,y in zip(yindep[md][chn][td],errindep[md][chn][td])],
                    color=colors[ichn], **opt_line)
            for itd,td in enumerate(('met', 'add_met_tau')):
                p2.circle([x+shift_both[td][ichn] for x in linear_x],
                          yboth[md][chn][td], color=colors[ichn], fill_alpha=1.,
                          **opt_points)
                p2.line([x+shift_both[td][ichn] for x in linear_x],
                        yboth[md][chn][td], color=colors[ichn], line_dash=styles[itd],
                        legend_label=chn+legends[td], **opt_line)
                p2.multi_line(
                    [(x+shift_both[td][ichn],x+shift_both[td][ichn]) for x in linear_x], 
                    [(max(0,x-y/2),min(100,x+y/2))
                     for x,y in zip(yboth[md][chn][td],errboth[md][chn][td])],
                    color=colors[ichn], **opt_line)
                
        for p in (p1, p2):
            p.xaxis[0].ticker = xticks
            p.xgrid[0].ticker = xticks
            p.xgrid.grid_line_alpha = 0.2
            p.xgrid.grid_line_color = 'black'
            p.yaxis[0].ticker = yticks
            p.ygrid[0].ticker = yticks
            p.ygrid.grid_line_alpha = 0.2
            p.ygrid.grid_line_color = 'black'
             
            p.legend.location = 'top_left'
            p.xaxis.axis_label = 'mHH [GeV]'
            p.yaxis.axis_label = 'Trigger acceptance gain (w.r.t. trigger baseline) [%]'
         
            p.xaxis.major_label_overrides = dict(zip(linear_x,x_str))
     
            p.legend.click_policy='hide'        
     
            p.output_backend = 'svg'
            #export_svg(p, filename='line_graph.svg')
         
        g = gridplot([[p1], [p2]])
        save(g, title=md)

if __name__ == '__main__':
    base_dir = '/eos/user/b/bfontana/www/TriggerScaleFactors/'
    main_dir = ['Region_190_190_PT_40_40_TURNON_200_190_NOTAU',]
    #'TriggerStudy_MET200_SingleTau190_CUT_entries_ditau_deltaR_GT_0_5',
    #'TriggerStudy_MET200_SingleTau190_CUT_entries_ditau_deltaR_ST_0_5']
    

    parser = argparse.ArgumentParser(description='Produce plots of trigger gain VS resonance mass.')

    parser.add_argument('--masses', required=True, nargs='+', type=str,
                        help='Resonance mass')
    parser.add_argument('--channels', required=True, nargs='+', type=str,  
                        help='Select the channel over which the workflow will be run.' )
    args = utils.parse_args(parser)
    main(args)
