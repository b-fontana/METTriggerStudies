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

def key_exists(d, k1, k2):
    return k1 in d and k2 in d[k1]

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
            
def set_quantiles(acc, quant, ntotal, nbatch):
    """
    Accumulates the weighted average of quantiles.
    - 'acc': accumulator
    - 'quant': new quantile being averaged
    """
    for k in acc:
        if acc[k] is None: # first iteration
            acc[k] = quant[k]
        else:
            frac_old = (ntotal-nbatch) / ntotal
            frac_new = nbatch / ntotal
            acc[k] = acc[k]*frac_old + quant[k]*frac_new
            
@utils.set_pure_input_namespace
def define_binning(args):
    """
    Determine histogram quantiles
    """
    quant_down, quant_up = 0., 1.
    cfg = importlib.import_module(args.configuration)

    channel_to_dataset = {'etau'   : 'EGamma',
                          'mutau'  : 'SingleMuon',
                          'tautau' : 'Tau',
                          'mumu'   : 'SingleMuon',
                          'ee'     : 'EGamma',}
    
    # Initialization
    # nTotEntries = 0
    # _maxedge, _minedge = ({} for _ in range(2))
    # for var in args.variables:
    #     _maxedge[var], _minedge[var] = ({} for _ in range(2))
    #     for chn in args.channels:
    #         _maxedge[var][chn], _minedge[var][chn] = ({} for _ in range(2))

    min_max, quants = {}, {}
    for k in args.channels:
        min_max[k], quants[k] = {}, {}
        for v in args.variables:
            min_max[k][v] = [1e10, -1e10]
            quants[k][v] = None
        
    ###############################################
    ############## Data Loop: Start ###############
    ###############################################
    if not skip_data_loop(args, cfg):

        for chn in args.channels:
            sample = channel_to_dataset[chn]
            
            #### Parse input list
            filelist, _ = utils.get_root_inputs(sample, args.indir, include_tree=True)

            treesize = 0
            percentage = 0.1

            quantiles = np.linspace(0., 1., num=args.nbins+1)
            quantiles[-1] = 0.99
                
            branches = args.variables + ('pairType',)
            nfiles = len(filelist)
            for ib,batch in enumerate(up.iterate(files=filelist, expressions=branches,
                                                 step_size='500 MB', library='ak')):
                treesize += len(batch)

                # assume X% of the files are enough to estimate quantiles
                if ib+1 > int(percentage*nfiles) and nfiles > 30:
                    break

                print('{}/{} {} files (only {} will be processed)\r'.format(ib+1, nfiles, sample,
                                                                            int(percentage*nfiles)),
                      end='' if ib+1!=nfiles else '\n', flush=True)
                
                if chn == 'all':
                    sel_chn = batch.pairType < main.sel[chn]['pairType'][1]
                else:
                    sel_chn = batch.pairType == main.sel[chn]['pairType'][1]
                batch_sel = batch[sel_chn]
                if key_exists(cfg.binedges, v, chn) and cfg.binedges[v][chn][0] == "quantiles":
                    for v in args.variables:
                        batch_sel[v] = batch_sel[v][batch_sel[v] > cfg.binedges[v][chn][1] &
                                                    batch_sel[v] < cfg.binedges[v][chn][2]]
                
                if len(batch_sel) != 0:
                    # trim outliers
                    quant1 = {v:np.quantile(batch_sel[v], q=np.array(quantiles)) for v in args.variables}
                    set_quantiles(quants[chn], quant1, ntotal=treesize, nbatch=len(batch))
                    quant2 = {v:np.quantile(batch_sel[v], q=np.array([0., 0.99])) for v in args.variables}
                    set_min_max(min_max[chn], quant2)
                else:
                    mes = ("Channel {} is not present in batch {} of sample {} in folders {}."
                           .format(chn, ib, sample, args.indir))
                    raise RuntimeError(mes)

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
                            elif key_exists(cfg.binedges, v, chn) and cfg.binedges[v][chn][0] == "quantiles":
                                dset = vgroup.create_dataset(chn, dtype=float, shape=(args.nbins+1,))
                                dset[:] = quants[chn][v]
                            else:
                                nb_tmp = len(cfg.binedges[v][chn])-1
                                dset = vgroup.create_dataset(chn, dtype=float, shape=(nb_tmp+1,))
                                dset[:] = cfg.binedges[v][chn]
                        else: # if the user did not specify the binning, equally spaced bins are used
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
