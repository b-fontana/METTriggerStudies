# coding: utf-8

_all_ = [ 'define_binning', 'define_binning_outputs' ]

import os
import sys
parent_dir = os.path.abspath(__file__ + 3 * '/..')
sys.path.insert(0, parent_dir)

import glob
import h5py
import uproot as up
import pandas as pd
import argparse

import inclusion
from inclusion.utils import utils
from inclusion import config

def skip_data_loop(args):
    for var in args.variables:
        if var not in config.binedges.keys():
            return False
    return True

@utils.set_pure_input_namespace
def define_binning_outputs(args):
    assert os.path.splitext(args.binedges_filename)[1] == '.hdf5'
    return os.path.join(args.outdir, args.binedges_filename)

@utils.set_pure_input_namespace
def define_binning(args):
    """
    Determine histogram quantiles
    """
    quant_down, quant_up = 0.00, .95

    # Initialization
    nTotEntries = 0
    _maxedge, _minedge = ({} for _ in range(2))
    for var in args.variables:
        _maxedge[var], _minedge[var] = ({} for _ in range(2))
        for chn in args.channels:
            _maxedge[var][chn], _minedge[var][chn] = ({} for _ in range(2))

    ###############################################
    ############## Data Loop: Start ###############
    ###############################################
    if not skip_data_loop(args):
        
        # Loop over all data datasets, calculating the quantiles in each dataset per variable
        for sample in args.data_vals:

            # to avoid long waiting times, define the binning using the first dataset only
            # the MET skims are 2GB, the EGamma one 25GB!
            if sample != args.data_vals[0]:
                continue
            
            #### Parse input list
            filelist = utils.get_root_input_files(sample, args.indir, include_tree=True)

            treesize = 0
            quantiles = {k: [] for k in args.channels }
            branches = args.variables + ('pairType',)
 
            for ib,batch in enumerate(up.iterate(files=filelist, expressions=branches,
                                                 step_size='10 GB', library='pd')):
                print('{}/{} {} files\r'.format(ib+1,len(filelist),sample),
                      end='' if ib+1!=len(filelist) else '\n', flush=True)
                treesize += batch.shape[0]
                df_chn = {}
                for chn in args.channels:
                    if chn == 'all':
                        sel_chn = batch['pairType'] <  config.sel[chn]['pairType'][1]
                    else:
                        sel_chn = batch['pairType'] == config.sel[chn]['pairType'][1]
                        
                    df_chn = batch[ sel_chn ]
                    quantiles[chn].append( df_chn.quantile([quant_down, quant_up]) )
                
            for chn in args.channels:
                quantiles[chn] = pd.concat(quantiles[chn], axis=0).groupby(level=0).mean()
            nTotEntries += treesize

            for var in args.variables:
                for chn in args.channels:
                    if (var not in config.binedges or
                        chn not in config.binedges[var]):
                        _minedge[var][chn][sample] = treesize*quantiles[chn].loc[quant_down, var]
                        _maxedge[var][chn][sample] = treesize*quantiles[chn].loc[quant_up, var]

        # Do weighted average based on the number of events in each dataset
        maxedge, minedge = ({} for _ in range(2))
        for var in args.variables:
            maxedge[var], minedge[var] = ({} for _ in range(2))
            for chn in args.channels:
                if var not in config.binedges or chn not in config.binedges[var]:
                    maxedge[var].update({chn: sum(_maxedge[var][chn].values()) / nTotEntries})
                    minedge[var].update({chn: sum(_minedge[var][chn].values()) / nTotEntries})
                    if maxedge[var][chn] <= minedge[var][chn]:
                        print('MaxEdge: {}, MinEdge: {}'.format(maxedge[var][chn], minedge[var][chn]))
                        print('Channel: {}, Var: {}'.format(chn, var))
                        raise ValueError('Wrong binning!')
                    if args.debug:
                        print( '{} quantiles in channel {}:  Q({})={}, Q({})={}'
                               .format(var, chn,
                                       quant_down, minedge[var][chn],
                                       quant_up,   maxedge[var][chn]) )
                else:
                    if args.debug:
                        print('Quantiles were not calculated. '
                              'Custom bins for variable {} and channel {} were instead used.'.format(var, chn))


    ###############################################
    ############## Data Loop: End #################
    ###############################################
    for _ in (True,): #breakable scope (otherwise 'break' cannot be used)
        with h5py.File( define_binning_outputs(args), 'a') as f:
            try:
                group = f.create_group(args.subtag)
            except ValueError:
                print('The HD5 group already existed. Skipping binning definition...')
                break
            for v in args.variables:
                vargroup = group.create_group(v)
                for chn in args.channels:
                    try:
                        if chn in config.binedges[v]:
                            dset = vargroup.create_dataset(chn, dtype=float,
                                                           shape=(len(config.binedges[v][chn]),))
                            if args.debug:
                                print( '[' + os.path.basename(__file__) + '] ' +
                                       'Using custom binning for variable {}: {}'
                                       .format(v, config.binedges[v][chn]) )
                            dset[:] = config.binedges[v][chn]
                        else:
                            raise KeyError #a "go-to" to the except clause

                    except KeyError:
                        dset = vargroup.create_dataset(chn, dtype=float, shape=(args.nbins+1,))
                        _binwidth = (maxedge[v][chn]-minedge[v][chn])/args.nbins
                        _data = [minedge[v][chn]+k*_binwidth for k in range(args.nbins+1)]
                        if args.debug:
                            print( '[' + os.path.basename(__file__) + '] Using regular binning for variable {}: {}'
                                   .format(v, _data) )
                        dset[:] = _data
                        

# -- Parse options
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Command line parser')

    parser.add_argument('--binedges_filename', dest='binedges_filename',
                        required=True, help='in directory')
    parser.add_argument('--nbins', dest='nbins', required=True, help='number of X bins')
    parser.add_argument('-a', '--indir', dest='indir', required=True, help='in directory')
    parser.add_argument('-o', '--outdir', dest='outdir', required=True, help='out directory')
    parser.add_argument('-t', '--tag', dest='tag', required=True, help='tag')
    parser.add_argument('--subtag', dest='subtag', required=True, help='subtag')
    parser.add_argument('--data_vals', dest='data_vals', required=True, nargs='+', type=str,
                        help='list of datasets')                          
    parser.add_argument('--variables', dest='variables', required=True, nargs='+', type=str,
                        help='Select the variables to calculate the binning' )
    parser.add_argument('-c', '--channels',   dest='channels',         required=True, nargs='+', type=str,
                        help='Select the channels over which the workflow will be run.' )
    parser.add_argument('--debug', action='store_true', help='debug verbosity')
    args = utils.parse_args(parser)

    define_binning( args )
