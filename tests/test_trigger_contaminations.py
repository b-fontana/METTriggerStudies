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
from bokeh.plotting import figure, output_file, save

tau = '\u03C4'
mu  = '\u03BC'
pm  = '\u00B1'
ditau = tau+tau

def main(args):
    output_file( os.path.join(basedir, 'contaminations.html') )
    p_opt = dict(width=800, height=400, x_axis_label='x', y_axis_label='y')
    p = figure(title='Contaminations in the Single Tau regions', tools='save', **p_opt)
    p.xaxis.axis_label = 'm(HH) [GeV]'
    p.yaxis.axis_label = 'Contamination [%]'
    p.toolbar.logo = None
    # p.x_range = Range1d(0, 11.5)
    # p.y_range = Range1d(0, 10)

    #linearize x axis
    linear_x = [k for k in range(1,len(args.masses)+1)]
    xticks = linear_x[:]
    p.xaxis[0].ticker = xticks
    p.xgrid[0].ticker = xticks
    p.xgrid.grid_line_alpha = 0.2
    p.xgrid.grid_line_color = 'black'

    shifts = (-0.15,-0.05,0.05,0.15)
    colors = ('red', 'green', 'blue', 'brown')
    for icut, cut in enumerate(args.region_cuts):
        label = '_'.join((str(cut),str(cut),args.channel))
        with h5py.File(os.path.join('data', label + '.hdf5'), 'r') as f:
            masses  = f[label][0]
            contam1 = f[label][1]
            contam2 = f[label][2]
            errors1 = f[label][3]
            errors2 = f[label][4]

            p.square([x+shifts[icut] for x in linear_x], contam1, fill_alpha=1., size=6,
                     color=colors[icut],
                     legend_label=label[:3] + ' ' + ditau)
            p.triangle([x+shifts[icut] for x in linear_x], contam2, fill_alpha=1., size=8,
                       color=colors[icut], alpha=0.5,
                       legend_label=label[:3] + ' ' + ditau + '+MET')
            p.line([x+shifts[icut] for x in linear_x], contam1,
                   color=colors[icut], line_width=1)
            p.line([x+shifts[icut] for x in linear_x], contam2, alpha=0.5,
                   color=colors[icut], line_width=1, line_dash='dashed',)
            p.multi_line([(x+shifts[icut],x+shifts[icut]) for x in linear_x], 
                         [(max(0,x-y/2),min(100,x+y/2)) for x,y in zip(contam1,errors1)],
                         color=colors[icut], line_width=2)
            p.multi_line([(x+shifts[icut],x+shifts[icut]) for x in linear_x], 
                         [(max(0,x-y/2),min(100,x+y/2)) for x,y in zip(contam2,errors2)],
                         color=colors[icut], line_width=2, alpha=0.5)

            x_str = [str(int(k)) for k in masses]
            p.xaxis.major_label_overrides = dict(zip(linear_x,x_str))

    p.legend.click_policy = 'hide'
    p.output_backend = 'svg'
    save(p)

if __name__ == '__main__':
    basedir = '/eos/user/b/bfontana/www/TriggerScaleFactors/'

    parser = argparse.ArgumentParser(description='Produce plots of trigger gain VS resonance mass.')
    parser.add_argument('--masses', required=True, nargs='+', type=str,
                        help='Resonance mass')
    parser.add_argument('--region_cuts', required=True, nargs='+', type=str,
                        help='Region SingleTau pT cuts')
    parser.add_argument('--channel', required=True, type=str,  
                        help='Select the channel over which the workflow will be run.' )
    args = utils.parse_args(parser)
    
    main(args)
