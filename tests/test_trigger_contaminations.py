# coding: utf-8

_all_ = [ 'test_trigger_contaminations' ]

import os
import sys
parent_dir = os.path.abspath(__file__ + 2 * '/..')
sys.path.insert(0, parent_dir)

import inclusion
from inclusion.config import main
from inclusion.utils import utils

import argparse
import h5py
import numpy as np
from bokeh.plotting import figure, output_file, save
from bokeh.palettes import Set1 as ColorSet
from bokeh.models import Range1d, ColumnDataSource
from bokeh.layouts import layout

tau = '\u03C4'
mu  = '\u03BC'
pm  = '\u00B1'
ditau = tau+tau

def main(args):
    output_file( os.path.join(basedir, 'contaminations_' + args.region_vary + '_' + args.channel + '.html') )
    p_opt = dict(width=800, height=400, x_axis_label='x', y_axis_label='y')
    title_d = {'dau1_pt': '(varying first lepton pT cut, second set to 190 GeV)',
               'dau2_pt': '(varying second lepton pT cut, first set to 190 GeV)',
               'both': '(varying both lepton pT cuts)'}
    p1 = figure(title='Contaminations in the SingleTau regions '+title_d[args.region_vary],
                tools='save,box_zoom,reset', **p_opt)
    p1.xaxis.axis_label = 'm(HH) [GeV]'
    p1.yaxis.axis_label = 'Contamination [%]'
    p1.x_range = Range1d(0, 14.5)
    p1.y_range = Range1d(-0.5, 20.)

    p2 = figure(title='Statistics in the Single and DiTau regions '+title_d[args.region_vary],
                tools='save,box_zoom,reset', **p_opt)
    p2.xaxis.axis_label = 'm(HH) [GeV]'
    p2.yaxis.axis_label = '#Events'
    # p2.x_range = Range1d(0, 14.5)
    # p2.y_range = Range1d(-0.5, 20.)

    #linearize x axis
    linear_x = [k for k in range(1,len(args.masses)+1)]
    xticks = linear_x[:]
    for p in (p1,p2):
        p.toolbar.logo = None
        p.xaxis[0].ticker = xticks
        p.xgrid[0].ticker = xticks
        p.xgrid.grid_line_alpha = 0.2
        p.xgrid.grid_line_color = 'black'

    nshifts = len(args.region_cuts)
    ns2 = int(nshifts/2)
    diff = 0.10 if ns2%2==0 else 0.05
    shifts = [round((-ns2+x)*diff,2) for x in range(nshifts)]
    colors = ColorSet[9]
    for icut, cut in enumerate(args.region_cuts):
        if args.region_vary == 'dau1_pt':
            cut1, cut2 = cut, '190'
        elif args.region_vary == 'dau2_pt':
            cut1, cut2 = '190', cut
        elif args.region_vary == 'both':
            cut1, cut2 = cut, cut
        label = '_'.join((str(cut1),str(cut2),args.channel))
        label_stats = label + '_stats'
        
        with h5py.File(os.path.join('data', label + '.hdf5'), 'r') as f:
            masses  = f[label][0]
            x_str = [str(int(k)) for k in masses]
            contam1 = f[label][1]
            contam2 = f[label][2]
            econtam1 = f[label][3]
            econtam2 = f[label][4]
            np.all(f[label_stats][1]==masses)

            stats1  = f[label_stats][1]
            estats1 = f[label_stats][2]

            source = ColumnDataSource({'x': [x+shifts[icut] for x in linear_x],
                                       'contam1': contam1,
                                       'contam2': contam2,
                                       'stats1': stats1,})
            opt = dict(color=colors[icut], source=source)

            leg0 = (label[:3] if args.region_vary=='dau1_pt' else label[4:7]) + 'GeV '
            leg1 = leg0 + ditau
            leg2 = leg1 + '+MET'
            p1.line('x', 'contam1', line_width=1, legend_label=leg1, **opt)
            p1.line('x', 'contam2', line_width=1, legend_label=leg2, line_dash='4 4', **opt)
            p1.multi_line([(x+shifts[icut],x+shifts[icut]) for x in linear_x], 
                          [(max(0,x-y/2),min(100,x+y/2)) for x,y in zip(contam1,econtam1)],
                          color=colors[icut], line_width=2, legend_label=leg1)
            p1.multi_line([(x+shifts[icut],x+shifts[icut]) for x in linear_x], 
                          [(max(0,x-y/2),min(100,x+y/2)) for x,y in zip(contam2,econtam2)],
                          color=colors[icut], line_width=2, legend_label=leg2)
            ga = p1.circle('x', 'contam1', size=6, legend_label=leg1, **opt)
            gb = p1.triangle('x', 'contam2', size=8, legend_label=leg2, **opt)
            ga = ga.glyph
            gb = gb.glyph
            ga.line_color = 'black'
            gb.line_color = 'black'
            p1.xaxis.major_label_overrides = dict(zip(linear_x,x_str))

            # p2.line('x', 'stats1', line_width=1, legend_label=leg0, **opt)
            p2.multi_line([(x+shifts[icut],x+shifts[icut]) for x in linear_x], 
                          [(max(0,x-y/2),x+y/2) for x,y in zip(stats1,estats1)],
                          color=colors[icut], line_width=2, legend_label=leg0)
            ga = p2.circle('x', 'stats1', size=6, legend_label=leg0, **opt)
            ga = ga.glyph
            ga.line_color = 'black'
            p2.xaxis.major_label_overrides = dict(zip(linear_x,x_str))

    for p in (p1,p2):
        p.legend.glyph_height = 15
        p.legend.glyph_width = 15
        p.legend.label_height = 9
        p.legend.label_width = 15
        p.legend.label_text_font_size = '7pt'
        p.legend.click_policy = 'hide'
        p.output_backend = 'svg'

    save(layout([[p1],[p2]]))


if __name__ == '__main__':
    basedir = '/eos/user/b/bfontana/www/TriggerScaleFactors/'

    parser = argparse.ArgumentParser(description='Produce plots of trigger gain VS resonance mass.')
    parser.add_argument('--masses', required=True, nargs='+', type=str,
                        help='Resonance mass')
    parser.add_argument('--region_cuts', required=True, nargs='+', type=str,
                        help='Region SingleTau pT cuts')
    parser.add_argument('--region_vary', required=True, type=str, choices=('dau1_pt', 'dau2_pt', 'both'),
                        help='SingleTau region to vary')
    parser.add_argument('--channel', required=True, type=str,  
                        help='Select the channel over which the workflow will be run.' )
    args = utils.parse_args(parser)
    
    main(args)
