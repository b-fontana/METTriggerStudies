# coding: utf-8

_all_ = [ 'discriminator' ]

import os
import sys
parent_dir = os.path.abspath(__file__ + 3 * '/..')
sys.path.insert(0, parent_dir)

import inclusion
from inclusion.utils import utils
from inclusion.utils.utils import join_name_trigger_intersection as joinNTC

import glob
import re
import json
import argparse

# This function is very much prone to change depending on the results obtained
def discriminator(args, chn):
    """
    Associates each trigger combination to a set of variables, ordered by importance.
    The variables will be used to retrieve the corresponding trigger efficiencies when evaluating the scale factors.
    """
    result = {}

    list_1D = {'etau': discr_vars_1D_etau,
    constant_list = ['dau1_pt', 'dau1_eta', 'dau2_pt', 'dau2_eta']
    for tcomb in utils.generate_trigger_combinations(chn, args.triggers):        
        result[joinNTC(tcomb)] = [ constant_list, #always the same 1D variables
                                   [utils.discr_vars_1D[chn][tcomb]], #1D changing variables
                                   [], #2D pairs of changing variables
                                  ]

    return result

def discriminator_exec_outputs(args, chn):
    out = os.path.join(args.outdir, '{}_{}.json'.format( os.path.basename(__file__).split('.')[0], chn))
    return out

def discriminator_exec(args, chn):
    match = re.compile('')
    
    if args.debug:
        print('[=debug=] Open file: {}'.format(name_data))

    out = discriminator_exec_outputs(args, chn)
    ordered_vars = discriminator(args, chn)

    with open(out, 'w') as f:
        json.dump(ordered_vars, f)

    print('File {} produced.'.format(out))


parser = argparse.ArgumentParser(description='Choose the most significant variables to draw the efficiencies.')

parser.add_argument('--indir',  help='Inputs directory',  required=True)
parser.add_argument('--outdir', help='Outputs directory', required=True)
parser.add_argument('--triggers', dest='triggers', nargs='+', type=str,
                    required=True, help='Triggers included in the workfow.')
parser.add_argument('--channel',   dest='channel', required=True,
                    help='Select the channel over which the discrimination will be run.' )
parser.add_argument('--variables',   dest='variables', required=True, nargs='+', type=str,
                    help='Select the variables over which the workflow will be run.' )
parser.add_argument('--tag', help='string to differentiate between different workflow runs', required=True)
parser.add_argument('--subtag', dest='subtag', required=True, help='subtag')
parser.add_argument('--debug', action='store_true', help='debug verbosity')
args = utils.parse_args(parser)

discriminator_exec(args, args.channel)
