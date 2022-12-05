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

import bokeh
from bokeh.plotting import figure, output_file, save
output_file('plots/trigger_gains.html')
from bokeh.layouts import gridplot
#from bokeh.io import export_svg

def main(args):
    channels = args.channels
    linear_x = [k for k in range(1,len(args.masses)+1)]
     
    yindep, yboth = ({} for _ in range(2))
    for md in main_dir:
        d_base = Path(base_dir) / md
        output_file(d_base / 'trigger_gains.html')
        yindep[md], yboth[md] = ({} for _ in range(2))
        for chn in channels:
            yindep[md][chn], yboth[md][chn] = ({} for _ in range(2))
            yindep[md][chn]['met'], yindep[md][chn]['tau'] = ([] for _ in range(2))
            yboth[md][chn]['met'],  yboth[md][chn]['add_met_tau']  = ([] for _ in range(2))
  
            for mass in args.masses:
                d = d_base / chn / str(mass)
                fullpath = Path(d) / 'baseline/counts/table.csv'
  
                with open(fullpath) as f:
                    reader = csv.reader(f, delimiter=',', quotechar='|')
                    next(reader, None) #ignore header line
                    ditau_line = next(reader, None)
                    next(reader, None) #ignore met region line
                    tau_line = next(reader, None)
  
                    sum_base     = float(ditau_line[1])
                    sum_met      = float(ditau_line[7]) + float(tau_line[7])
                    sum_tau      = float(ditau_line[9]) + float(tau_line[9])
                    sum_only_tau = float(ditau_line[8]) + float(tau_line[8])
                        
                    frac_indep_met = sum_met / sum_base
                    frac_indep_tau = sum_tau / sum_base
                    print(sum_tau, sum_base)
                    frac_both  = (sum_met + sum_only_tau) / sum_base
                    yindep[md][chn]['met'].append(frac_indep_met*100)
                    yindep[md][chn]['tau'].append(frac_indep_tau*100)
                    yboth[md][chn]['met'].append(frac_indep_met*100)
                    yboth[md][chn]['add_met_tau'].append(frac_both*100)

    print(yindep[md][chn]['tau'])
    quit()

    opt_points = dict(size=8)
    opt_line = dict(width=1.5)
    colors = ('green', 'blue', 'red')
    styles = ('solid', 'dashed')
    legends = {'met': ' (MET)', 'tau': ' (Tau)', 'add_met_tau': ' (MET + Tau)'}
     
    x_str = [str(k) for k in args.masses]
    xticks = linear_x[:]
    yticks = [x for x in range(0,100,10)]
     
    for md in main_dir:
        p_opt = dict(width=800, height=400, x_axis_label='x', y_axis_label='y')
        p1 = figure(title='Inclusion of MET and Single Tau triggers', **p_opt)
        p2 = figure(title='Acceptance gain of MET + SingleTau triggers', **p_opt)
        p1.toolbar.logo = None
        p2.toolbar.logo = None
        for ichn,chn in enumerate(channels):
            for itd,td in enumerate(('met', 'tau')):
                p1.circle(linear_x, yindep[md][chn][td], color=colors[ichn], fill_alpha=1.,
                          **opt_points)
                p1.line(linear_x, yindep[md][chn][td], color=colors[ichn], line_dash=styles[itd],
                        legend_label=chn+legends[td], **opt_line)
            for itd,td in enumerate(('met', 'add_met_tau')):
                p2.circle(linear_x, yboth[md][chn][td], color=colors[ichn], fill_alpha=1.,
                          **opt_points)
                p2.line(linear_x, yboth[md][chn][td], color=colors[ichn], line_dash=styles[itd],
                        legend_label=chn+legends[td], **opt_line)
                
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
