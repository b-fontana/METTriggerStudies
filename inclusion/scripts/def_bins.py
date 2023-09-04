# coding: utf-8

_all_ = [ 'define_binning', 'define_binning_outputs' ]

import os
import sys
parent_dir = os.path.abspath(__file__ + 3 * '/..')
sys.path.insert(0, parent_dir)

import glob
import h5py
import uproot as up
import numpy as np
import argparse
import importlib

import inclusion
from inclusion.utils import utils
from inclusion.config import main

def skip_data_loop(args, cfg):
    for var in args.variables:
        if var not in cfg.binedges.keys():
            return False
    return True

@utils.set_pure_input_namespace
def define_binning_outputs(args):
    assert os.path.splitext(args.binedges_filename)[1] == '.hdf5'
    return os.path.join(args.outdir, args.binedges_filename)

def set_min_max(acc, obj):
    for k in acc:
        if obj[k][0] < acc[k][0]:
            acc[k][0] = obj[k][0]
        if obj[k][1] > acc[k][1]:
            acc[k][1] = obj[k][1]
            
@utils.set_pure_input_namespace
def define_binning(args):
    """
    Determine histogram quantiles
    """
    quant_down, quant_up = 0., 1.
    cfg = importlib.import_module(args.configuration)
    
    # Initialization
    # nTotEntries = 0
    # _maxedge, _minedge = ({} for _ in range(2))
    # for var in args.variables:
    #     _maxedge[var], _minedge[var] = ({} for _ in range(2))
    #     for chn in args.channels:
    #         _maxedge[var][chn], _minedge[var][chn] = ({} for _ in range(2))

    min_max = {}
    for k in args.channels:
        min_max[k] = {}
        for v in args.variables:
            min_max[k][v] = [1e10, -1e10]
        
    ###############################################
    ############## Data Loop: Start ###############
    ###############################################
    if not skip_data_loop(args, cfg):

        # Loop over all data datasets, calculating the quantiles in each dataset per variable
        for sample in args.data_vals:

            # to avoid long waiting times, define the binning using the first dataset only
            # the MET skims are 2GB, the EGamma one 25GB!
            if sample != args.data_vals[0]:
                continue
            
            #### Parse input list
            filelist, _ = utils.get_root_inputs(sample, args.indir, include_tree=True)

            treesize = 0
            percentage = 0.1
            qup, qlo = 0.99, 0.
            branches = args.variables + ('pairType',)
            nfiles = len(filelist)
            for ib,batch in enumerate(up.iterate(files=filelist, expressions=branches,
                                                 step_size='50 MB', library='ak')):
                if ib+1 > int(percentage*nfiles) and nfiles > 30: # assume X% of the files are enough to estimate quantiles
                    break

                print('{}/{} {} files (only {} will be processed)\r'.format(ib+1, nfiles, sample, int(percentage*nfiles)),
                      end='' if ib+1!=nfiles else '\n', flush=True)

                
                treesize += len(batch)
                for chn in args.channels:
                    if chn == 'all':
                        sel_chn = batch.pairType < main.sel[chn]['pairType'][1]
                    else:
                        sel_chn = batch.pairType == main.sel[chn]['pairType'][1]
                    batch_sel = batch[sel_chn]
                    if len(batch_sel) != 0:
                        # trim outliers
                        quant = {v:np.quantile(batch_sel[v], q=np.array([qlo, qup])) for v in args.variables}
                        set_min_max(min_max[chn], quant)
                    else:
                        mes = ("Channel {} is not present in batch {} of sample {} in folders {}."
                               .format(chn, ib, sample, args.indir))
                        raise RuntimeError(mes)

            # for var in args.variables:
            #     for chn in args.channels:
            #         if (var not in cfg.binedges or chn not in cfg.binedges[var]):
            #             _minedge[var][chn][sample] = treesize*quantiles[chn].loc[quant_down, var]
            #             _maxedge[var][chn][sample] = treesize*quantiles[chn].loc[quant_up, var]

            # nTotEntries += treesize
            
        # Do weighted average based on the number of events in each dataset
        # maxedge, minedge = ({} for _ in range(2))
        # for var in args.variables:
        #     maxedge[var], minedge[var] = ({} for _ in range(2))
        #     for chn in args.channels:
        #         if var not in cfg.binedges or chn not in cfg.binedges[var]:
        #             maxedge[var].update({chn: sum(_maxedge[var][chn].values()) / nTotEntries})
        #             minedge[var].update({chn: sum(_minedge[var][chn].values()) / nTotEntries})
        #             if maxedge[var][chn] <= minedge[var][chn]:
        #                 print('MaxEdge: {}, MinEdge: {}'.format(maxedge[var][chn], minedge[var][chn]))
        #                 print('Channel: {}, Var: {}'.format(chn, var))
        #                 raise ValueError('Wrong binning!')
        #             if args.debug:
        #                 print( '{} quantiles in channel {}:  Q({})={}, Q({})={}'
        #                        .format(var, chn,
        #                                quant_down, minedge[var][chn],
        #                                quant_up,   maxedge[var][chn]) )
        #         else:
        #             if args.debug:
        #                 print('Quantiles were not calculated. '
        #                       'Custom bins for variable {} and channel {} were instead used.'.format(var, chn))


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
                vgroup = group.create_group(v)
                for chn in args.channels:
                    try:
                        if chn in cfg.binedges[v]:                            
                            if args.debug:
                                mes = 'Using custom binning for variable {}: {}'.format(v, cfg.binedges[v][chn])
                                utils.debug(mes, args.debug, __file__)

                            if len(cfg.binedges[v][chn]) < 2:
                                mes = 'The binning of variable {} must have at least two values.'.format(v)
                                raise ValueError(mes)
                            elif len(cfg.binedges[v][chn]) == 2:
                                dset = vgroup.create_dataset(chn, dtype=float, shape=(args.nbins+1,))
                                dset[:] = np.linspace(cfg.binedges[v][chn][0], cfg.binedges[v][chn][1], args.nbins+1)
                            else:
                                nb_tmp = len(cfg.binedges[v][chn])-1
                                dset = vgroup.create_dataset(chn, dtype=float, shape=(nb_tmp+1,))
                                dset[:] = cfg.binedges[v][chn]
                        else:
                            raise KeyError #a "go-to" to the except clause

                    except KeyError:
                        dset = vgroup.create_dataset(chn, dtype=float, shape=(args.nbins+1,))
                        _binwidth = (min_max[chn][v][1]-min_max[chn][v][0])/args.nbins
                        _data = [min_max[chn][v][0]+k*_binwidth for k in range(args.nbins+1)]
                        if args.debug:
                            mes = 'Using regular binning for variable {}: {}'.format(v, _data)
                            utils.debug(mes, args.debug, __file__)
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
    parser.add_argument('--configuration', dest='configuration', required=True,
                        help='Name of the configuration module to use.')
    parser.add_argument('--debug', action='store_true', help='debug verbosity')
    args = utils.parse_args(parser)

    define_binning( args )
