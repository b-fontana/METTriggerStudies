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

def main(args):
    output_file( os.path.join(basedir, 'contaminations.html') )
    p_opt = dict(width=800, height=400, x_axis_label='x', y_axis_label='y')
    p = figure(title='Contaminations in the Single Tau regions', tools='save', **p_opt)
    p.toolbar.logo = None
    p.xaxis.axis_label = 'm(HH) [GeV]'
    p.yaxis.axis_label = 'Contamination [%]'
    p.output_backend = 'svg'

    colors = ['orange', 'green', 'brown', 'purple']
    for icut, cut in enumerate(args.region_cuts):
        label = str(cut) + '_' + str(cut)
        with h5py.File(os.path.join('data', label + '.hdf5'), 'r') as f:
            masses = f[label][0]
            contam = f[label][1]
            errors = f[label][2]
            p.square(masses, contam, fill_alpha=1., size=6,
                     color=colors[icut],
                     legend_label='tau region contamin. by met')
            p.line(masses, contam, color=colors[icut], line_width=1)
            p.multi_line([(x,x) for x in masses], 
                         [(max(0,x-y/2),min(100,x+y/2)) for x,y in zip(contam,errors)],
                         color=colors[icut], line_width=2)

    save(p)

if __name__ == '__main__':
    basedir = '/eos/user/b/bfontana/www/TriggerScaleFactors/'

    parser = argparse.ArgumentParser(description='Produce plots of trigger gain VS resonance mass.')

    parser.add_argument('--masses', required=True, nargs='+', type=str,
                        help='Resonance mass')
    parser.add_argument('--region_cuts', required=True, nargs='+', type=str,
                        help='Region SingleTau pT cuts')
    parser.add_argument('--channels', required=True, nargs='+', type=str,  
                        help='Select the channel over which the workflow will be run.' )
    args = utils.parse_args(parser)
    main(args)
