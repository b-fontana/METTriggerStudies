import os
import re
import argparse
import ctypes
import numpy as np
from copy import copy

import ROOT
ROOT.gROOT.SetBatch(True)
from ROOT import (
    TCanvas,
    TEfficiency,
    TFile,
    TGraphAsymmErrors,
    TH1D,
    TH2D,
    TLatex,
    TLegend,
    TLine,
    TPad,
    )
from array import array

import sys
sys.path.append( os.path.join(os.environ['CMSSW_BASE'], 'src', 'METTriggerStudies'))
from utils.utils import (
    create_single_dir,
    get_display_variable_name,
    get_histo_names,
    get_key_list,
    get_root_object,
    load_binning,
    print_configuration,
    redraw_border,
    rewrite_cut_string,
    uniformize_bin_width,
    )
from luigi_conf import (
    _extensions,
    _placeholder_cuts,
    )

def convertGraphToHist(graph):
    """
    Converts 1D TGraphAsymmErrors to TH1D.
    TODO: support 2D graphs [3D histograms hopefully not...]
    The convertion is required due to ROOT's bad support of 2D objects:
    - TGraph2DAsymmErrors must have the same number of points in X and Y!
    - TGraph3DAsymmErrors does not exist
    - Histograms are the most flexible, although do not support asym errors.
    """
    n = graph.GetN()
    name = graph.GetName()

    # define x (bin centers)
    x = array('d')
    x.append((3.0 * graph.GetX()[0] - graph.GetX()[1]) / 2.0)
    for i in range(1,n):
        x.append((graph.GetX()[i-1] + graph.GetX()[i]) / 2.0)
    x.append((3.0 * graph.GetX()[n-1] - graph.GetX()[n-2]) / 2.0)
    histo = ROOT.TH1D(name, name, n, x)

    # fill y values (counts)
    for i in range(n):
        histo.SetBinContent(i+1, graph.GetY()[i])
        # losing up/down asymmetric info...
        toterr = graph.GetErrorYhigh(i)+graph.GetErrorYlow(i)
        histo.SetBinError(i+1, toterr)

    return histo

def aggregateEfficienciesAndScaleFactors(indir, outdir, channel, subtag, prefix, variables, debug):
    extension = 'root'
    _outname = os.path.join(outdir, prefix + channel + '.' + extension)
    fout = ROOT.TFile.Open(_outname , 'RECREATE')
    
    walk_path = os.path.join(indir, channel)
    var_re = '(?:' + '|'.join(variables) + ')' # ?: match but do not capture

    regex_str = ( prefix + '.*' + channel + '.*(' + var_re +
                 ')_(.+)' + args.subtag + '_CUTS_(.+' + ')\.' + extension)
    regex = re.compile( regex_str )
    
    for root, _, files in os.walk( walk_path ):
        for afile in files:
            if afile.endswith('.' + extension):
                froot = ROOT.TFile.Open( os.path.join(root,afile), 'READ' )
                keyList = ROOT.TIter(froot.GetListOfKeys())
                fout.cd()
                for key in keyList:
                    cl = ROOT.gROOT.GetClass(key.GetClassName())
                    if not cl.InheritsFrom('TGraph'):
                        continue
                    h = key.ReadObj()
                    h = convertGraphToHist(h)
                    
                    var, trigger, cut = regex.match(afile).groups()
                    new_name = h.GetName() + '_VAR_' + var + '_TRG_' + trigger + '_CUT_' + cut

                    h.Write(new_name)
            
parser = argparse.ArgumentParser(description='Draw trigger scale factors')

parser.add_argument('--indir', help='Inputs directory', required=True)
parser.add_argument('--outdir', help='Output directory', required=True, )
parser.add_argument('--channel', dest='channel', required=True, type=str,
                    help='Select the channels over which the workflow will be run.' )
parser.add_argument('--file_prefix', help='ROOT file prefix', required=True)
parser.add_argument('--variables', dest='variables', required=True, type=str, nargs='+',
                    help='Select the variables over which the workflow will be run.' )
parser.add_argument('--subtag', dest='subtag', required=True, help='subtag')
parser.add_argument('--debug', action='store_true', help='debug verbosity')
args = parser.parse_args()
print_configuration(args)

aggregateEfficienciesAndScaleFactors(args.indir, args.outdir, args.channel, args.subtag,
                                     args.file_prefix, args.variables, args.debug)
