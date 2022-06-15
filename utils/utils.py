import os
import glob
import inspect
import re
import operator
import argparse
import functools
import itertools as it
import numpy as np
import h5py
from types import SimpleNamespace

import ROOT
from ROOT import (
    TGraphAsymmErrors,
    TH2D,
    TLine,
)

from luigi_conf import (
    _placeholder_cuts as pholdc,
    _sel,
    _triggers_map,
    _triggers_custom,
)
  
def add_slash(s):
    """Adds single slash to path if absent"""
    s = s if s[-1] == '/' else s + '/'
    return s

def add_vnames(*vnames):
    return '_VERSUS_'.join(vnames)

def at_least_two(x1, x2, x3):
    """Checks if at least two out of the three boleans are True."""
    return x1 if (x2 or x3) else (x2 and x3)
    
def build_prog_path(base, script_name):
    script = os.path.join(base, 'scripts')
    script = os.path.join(script, script_name)
    return 'python3 {}'.format(script)

def check_bit(number, bitpos):
    bitdigit = 1
    res = bool(number&(bitdigit<<bitpos))
    return res
            
def create_single_dir(p):
    """Creates a directory if it does not exist"""
    try:
        if not os.path.exists(p): 
            os.makedirs(p)
    except PermissionError:
        m = ( "You tried to create folder {}. Make sure you have the rights!"
              "Are you sure the path is correct?".format(p) )
        print(m)
        raise
    except FileExistsError:
        pass
    
def create_single_file(f):
    """Creates a dummy file if it does not exist"""
    try:
        os.remove(f)
    except OSError as e:
        if e.errno != errno.ENOENT: # errno.ENOENT = no such file or directory
            raise
    with open(f, 'x') as newf:
        newf.write('[utils] Dummy text.')

def debug(message, flag=True):
    decorator = ' ============ '
    if flag:
        print( decorator + message + decorator, flush=True )

class dot_dict(dict):
    """dot.notation access to dictionary attributes"""
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

def find_bin(edges, value, var):
    """
    Find the bin id corresponding to one value, given the bin edges.
    Bin 0 stands for underflow.
    """
    binid = np.digitize(value, edges)
    if binid == len(edges):
        binid -= 1 # include overflow

    # check for out-of-bounds
    if binid==0:
        print(binid, value, var)
        print(edges, len(edges))
        print( '[{}] WARNING: This should never happen in production code. '
               'Are you sure the binning is well defined?'.format(os.path.basename(__file__)) )
        binid = 1
    elif binid>len(edges):
        raise ValueError('[{}] WARNING: This cannot happen in production code.')

    return binid

def generate_trigger_combinations(channel, trigs):
    """
    Set all possible trigger intersection combinations, per channel.
    Does not consider intersections between incompatible triggers and channels
    (for instance, `etau` channel with `IsoMu24` trigger).
    Each intersection is sorted alphabetically
    (useful for matching with the KLUB framework).
    """
    exclusive = { 'etau'   : ('Ele32', 'Ele35', 'EleIsoTauCustom'),
                  'mutau'  : ('IsoMu24', 'IsoMu27', 'IsoMuIsoTauCustom'),
                  'tautau' : ('IsoDoubleTauCustom',) }

    # look only at combinations where the channel is imcompatible with the trigger
    pruntrigs = [ exclusive[x] for x in exclusive if x != channel ]
    pruntrigs = set(it.chain(*pruntrigs))
    pruntrigs = set(trigs) - pruntrigs

    length1 = list(it.chain.from_iterable(it.combinations(sorted(pruntrigs), 1)))
    for elem in length1:
        if elem not in _triggers_map.keys():
            mess = '[utils.generate_trigger_combinations] '
            mess += 'Trigger {} is not supported'.format(elem)
            raise ValueError(mess)

    complete = list( it.chain.from_iterable(it.combinations(sorted(pruntrigs), x)
                                            for x in range(1,len(pruntrigs)+1)) )
    return complete

def get_display_variable_name(channel, var):
    if channel == 'mutau':
        var_custom = var.replace('dau1', 'mu').replace('dau2', 'tau')
    elif channel == 'etau':
        var_custom = var.replace('dau1', 'electron').replace('dau2', 'electron')
    elif channel == 'tautau':
        var_custom = var.replace('dau1', 'tau').replace('dau2', 'tau')
    elif channel == 'mumu':
        var_custom = var.replace('dau1', 'mu').replace('dau2', 'mu')
    else:
        var_custom = var
    return var_custom

def get_key_list(afile, inherits=['TH1']):
    tmp = []
    keylist = ROOT.TIter(afile.GetListOfKeys())
    for key in keylist:
        cl = ROOT.gROOT.GetClass(key.GetClassName())

        inherited = functools.reduce( lambda x, y: x or y,
                                      [ cl.InheritsFrom(x) for x in inherits ] )
        if not inherited:
            continue
    
        h = key.ReadObj()
        tmp.append( h.GetName() )
    return tmp

def get_hnames(opt):
    if opt == 'Ref1D':
        return lambda a,b : 'Ref1D_{}_{}'.format(a,b)
    elif opt == 'Trig1D':
        return lambda a,b,c : 'Trig1D_{}_{}_{}{}'.format(a,b,c,pholdc)
    elif opt == 'Ref2D':
        return lambda a,b : 'Ref2D_{}_{}'.format(a,b)
    elif opt == 'Trig2D':
        return lambda a,b,c : 'Trig2D_{}_{}_{}{}'.format(a,b,c,pholdc)
    elif opt == 'Canvas2D':
        return lambda a,b,c : 'Canvas2D_{}_{}_{}{}'.format(a,b,c,pholdc)
    elif opt == 'Closure':
        return lambda a,b,c,d : 'Closure{}_{}_{}_{}'.format(a,b,c,d)
    else:
        import inspect
        currentFunction = inspect.getframeinfo(frame).function
        raise ValueError('[{}] option not supported.'.format(currentFunction))

def get_obj_max_min(graph, is_histo):
    vmax, vmin = 0, 1e10
    npoints = graph.GetN()
    for point in range(npoints):
        if is_histo:
            val = graph.GetBinContent(point+1)
        else:
            val = graph.GetPointY(point)
        if val > vmax:
            vmax = val
        if val < vmin:
            vmin = val
    return vmax, vmin

def get_root_input_files(proc, indir):
    #### Check input folder
    if not isinstance(indir, (tuple,list)):
        indir = [indir]
    fexists = []
    for idir in indir:
        g = glob.glob(os.path.join(idir, proc + '*/goodfiles.txt'))
        if len(g)==0:
            fexists.append(False)
        else:
            fexists.append(True)

    if sum(fexists) != 1: #check one and only one is True
        m = '[utils.py] WARNING: More than one file exists for the {} sample.'.format(proc)
        m += ' Selecting directory {} '.format(indir[fexists.index(True)])
        m += 'from the following list: {}.'.format(indir)

    #this is the only correct input directory
    inputdir = indir[ fexists.index(True) ]

    inputfiles = glob.glob( os.path.join(indir[fexists.index(True)],
                                         proc + '*/goodfiles.txt') )

    #### Parse input list
    filelist=[]
    for inp in inputfiles:
        with open(inp) as fIn:
            for line in fIn:
                if '.root' in line:
                    filelist.append(line)

    return filelist, inputdir

def get_root_object(name, afile):
    try:
        _keys = afile.GetListOfKeys()
    except ReferenceError:
        print('File {} does not exist!'.format(afile))
        raise
    if name not in _keys:
        msg =  'Wrong ROOT object name!\n'
        msg += 'File name: {}\n'.format(afile.GetName())
        msg += 'Object name: {}\n'.format(name)
        msg += '{} keys:\n'.format(len(_keys))
        for k in _keys[:5]:
            msg += ' - {}\n'.format(k.GetName())
        msg += ' - ...'
        raise ValueError(msg)
    return afile.Get(name)

def get_trigger_bit(trigger_name, isdata):
    """
    Returns the trigger bit corresponding to '_triggers_map'
    """
    s = 'data' if isdata else 'mc'
    res = _triggers_map[trigger_name]
    try:
        res = res[s]
    except KeyError:
        print('You likely forgot to add your custom trigger to _triggers_custom.')
        raise
    return res

def hadd_subpaths(args, channel=''):
    channel_str = '' if channel=='' else channel+'_'
    _tbase1 = args.tprefix + channel_str + args.dataset_name
    _tbase2 = '_Sum' + args.subtag
    return _tbase1, _tbase2

def is_channel_consistent(chn, pairtype):
    opdict = { '<':  operator.lt,
               '>':  operator.gt,
               '==': operator.eq }

    op, val = _sel[chn]['pairType']
    return opdict[op](pairtype, val)
  
def join_name_trigger_intersection(tuple_element):
    inters = '_PLUS_'
    return inters.join(tuple_element)

class LeafManager():
    """
    Class to manage TTree branch leafs, making sure they exist.
    """
    def __init__(self, fname, t_in):
        self.fname = fname
        self.tree = t_in
        self.absent_leaves = set()
        self.error_prefix = '[LeafManager]: '
        
    def get_leaf(self, leaf):
        if not isinstance(leaf, str):
            m = 'The leaf must be a string.'
            raise TypeError(self.error_prefix + m)
        try:
            obj = self.tree.GetListOfBranches().FindObject(leaf)
            name = obj.GetName()
            getAttr = lambda x : getattr(self.tree, x)
            return getAttr(leaf)
        except ReferenceError:
            if leaf not in self.absent_leaves:
                m = 'WARNING: leaf ' + leaf + ' does not exist in file ' + self.fname
                print(self.error_prefix + m)
                self.absent_leaves.add(leaf)
            return 0.

def load_binning(afile, key, variables, channels):
    """
    Load the Binning stored in a previous task.
    """
    binedges, nbins = ({} for _ in range(2))
    with h5py.File(afile, 'r') as f:
        try:
            group = f[key]
        except KeyError:
            print('{} does not have key {}.'.format(afile, key))
            print('Available keys: {}'.format(f.keys()))
            raise
          
        for var in variables:
            subgroup = group[var]
            binedges[var], nbins[var] = ({} for _ in range(2))
            for chn in channels:
                binedges[var][chn] = np.array(subgroup[chn][:])
                nbins[var][chn] = len(binedges[var][chn]) - 1

    return binedges, nbins

def pass_selection_cuts(leaf_manager, invert_mass_cut=True):
    """
    Applies selection cut per TTree entry.
    Returns `True` only if all selection is passed.
    """
    mhh = leaf_manager.get_leaf( 'HHKin_mass' )
    if mhh<1:
        return False

    pairtype = leaf_manager.get_leaf( 'pairType' )
    dau1_eleiso = leaf_manager.get_leaf( 'dau1_eleMVAiso'    )
    dau1_muiso  = leaf_manager.get_leaf( 'dau1_iso'          )
    dau1_tauiso = leaf_manager.get_leaf( 'dau1_deepTauVsJet' )
    dau2_tauiso = leaf_manager.get_leaf( 'dau2_deepTauVsJet' )

    # Loose / Medium / Tight
    bool0 = pairtype==0 and (dau1_muiso>=0.15 or dau2_tauiso<5)
    bool1 = pairtype==1 and (dau1_eleiso!=1 or dau2_tauiso<5)
    bool2 = pairtype==2 and (dau1_tauiso<5 or dau2_tauiso<5)
    if bool0 or bool1 or bool2:
        return False

    #((tauH_SVFIT_mass-116.)*(tauH_SVFIT_mass-116.))/(35.*35.) + ((bH_mass_raw-111.)*(bH_mass_raw-111.))/(45.*45.) <  1.0
    svfit_mass = leaf_manager.get_leaf('tauH_SVFIT_mass')
    bh_mass    = leaf_manager.get_leaf('bH_mass_raw')

    mcut = ( (svfit_mass-129.)*(svfit_mass-129.) / (53.*53.) +
             (bh_mass-169.)*(bh_mass-169.) / (145.*145.) ) <  1.0
    if mcut and invert_mass_cut: # inverted elliptical mass cut (-> ttCR)
        return False

    #pass_met = leaf_manager.get_leaf('isMETtrigger')
    #pass_tau = leaf_manager.get_leaf('isSingleTautrigger')
    #pass_taumet = leaf_manager.get_leaf('isTauMETtrigger')
    pass_lep = leaf_manager.get_leaf('isLeptrigger')
    if not pass_lep:
        return False

    return True

def redraw_border():
    """
    this little macro redraws the axis tick marks and the pad border lines.
    """
    ROOT.gPad.Update()
    ROOT.gPad.RedrawAxis()
    l = TLine()
    l.SetLineWidth(2)

    l.DrawLine(ROOT.gPad.GetUxmin(), ROOT.gPad.GetUymax(), ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymax()) #top border
    l.DrawLine(ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymin(), ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymax()) #right border
    l.DrawLine(ROOT.gPad.GetUxmin(), ROOT.gPad.GetUymin(), ROOT.gPad.GetUxmin(), ROOT.gPad.GetUymax()) #left border
    l.DrawLine(ROOT.gPad.GetUxmin(), ROOT.gPad.GetUymin(), ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymin()) #bottom border

def remove(f):
    if os.path.exists( f ):
        os.remove( f )

def rewrite_cut_string(oldstr, newstr, regex=False):
    if regex:
        _regex = re.findall(r'^.*CUTS_(.+)$', newstr)
        assert(len(_regex)==1)
        _regex = _regex[0]
        newstr = _regex.replace('>', 'L').replace('<', 'S').replace('.', 'p')
    
    res = oldstr.replace(pholdc, '_CUTS_'+newstr)
    return res

def set_pure_input_namespace(func):
    """
    Decorator which forces the input namespace to be a "bare" one.
    Used when luigi calls a function with a luigi.DictParameter().
    'param' can be used to pass additional arguments to the function being decorated.
    """
    def wrapper(args, param=None):
        if not isinstance(args, (argparse.Namespace, SimpleNamespace)):
            args = SimpleNamespace(**args)
        if param:
            return func(args, param)
        else:
            return func(args)

    return wrapper

def set_custom_trigger_bit(trigger, trigBit, run, isData):
    """
    The VBF trigger was updated during data taking, adding HPS
    https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauTrigger
    """
    if trigger not in _triggers_custom:
        import inspect
        currentFunction = inspect.getframeinfo(frame).function
        raise ValueError('[{}] option not supported.'.format(currentFunction))

    if run < 317509 and isData:
        if trigger == 'VBFTauCustom':
            bits = check_bit(trigBit, _triggers_map[trigger]['VBFTau']['data'])
        elif trigger == 'IsoDoubleTauCustom':
            bits = ( check_bit(trigBit, _triggers_map[trigger]['IsoDoubleTau']['data'][0]) or
                     check_bit(trigBit, _triggers_map[trigger]['IsoDoubleTau']['data'][1]) or
                     check_bit(trigBit, _triggers_map[trigger]['IsoDoubleTau']['data'][2]) )
        elif trigger == 'IsoMuIsoTauCustom':
            bits = check_bit(trigBit, _triggers_map[trigger]['IsoMuIsoTau']['data'])
        elif trigger == 'EleIsoTauCustom':
            bits = check_bit(trigBit, _triggers_map[trigger]['EleIsoTau']['data'])

    else:
        s = 'data' if isData else 'mc'
        if trigger == 'VBFTauCustom':
            bits = check_bit(trigBit, _triggers_map[trigger]['VBFTauHPS'][s])
        elif trigger == 'IsoDoubleTauCustom':
            bits = check_bit(trigBit, _triggers_map[trigger]['IsoDoubleTauHPS'][s])
        elif trigger == 'IsoMuIsoTauCustom':
            bits = check_bit(trigBit, _triggers_map[trigger]['IsoMuIsoTauHPS'][s])
        elif trigger == 'EleIsoTauCustom':
            bits = check_bit(trigBit, _triggers_map[trigger]['EleIsoTauHPS'][s])

    return bits

def split_vnames(joinvars):
    return joinvars.split('_VERSUS_')

def pass_any_trigger(trigs, bit, run, isdata):
    # checks that at least one trigger was fired
    flag = False
    for trig in trigs:
        if trig in _triggers_custom:
            flag = set_custom_trigger_bit(trig, bit, run, isdata)
        else:
            flag = check_bit(bit, get_trigger_bit(trig, isdata))
        if flag:
            return True
    return False

def pass_trigger_bits(trig, trig_bit, run, isdata):
    if trig in _triggers_custom:
        return set_custom_trigger_bit(trig, trig_bit, run, isdata)
    else:
        return check_bit(trig_bit, get_trigger_bit(trig, isdata))

def print_configuration(parse_args):
    filename = inspect.stack()[1].filename 
    
    print('----------------------------------------', flush=True)
    print('Script Name: {}'.format(filename))
    print('Script Configuration:', flush=True)
    options = vars(parse_args)
    maxlkey = max(len(x) for x in options.keys())
    for k,v in options.items():
        k = '--' + k
        if isinstance(v, (tuple,list)):
            v = ' '.join(v)
        print('{0:>{d1}}   {1}'.format(k, v, d1=maxlkey+3), flush=True)
    print('----------------------------------------', flush=True)

def parse_args(parser):
    args = parser.parse_args()
    print_configuration(args)
    return args
    
def slash_to_underscore_and_keep(s, n=4):
    """Replaces slashes by underscores, keeping only the last 'n' slash-separated strings"""
    return '_'.join( s.split('/')[-n:] )

def apply_equal_bin_width(old, roundx=2, roundy=2):
    """
    Replaces the object by changing the labels and
    ensuring all bins have the same visual width.
    Internally the bins are numbered from 0 to nBins.
    """
    darr = lambda x : np.array(x).astype(dtype=np.double)
    def xfunc(l, r):
        lstr = str(round(l, roundx)) if roundx!=0 else str(int(l))
        rstr = str(round(r, roundx)) if roundx!=0 else str(int(r))
        return '[' + lstr + ';' + rstr + '['
    def yfunc(l, r):
        lstr = str(round(l, roundy)) if roundy!=0 else str(int(l))
        rstr = str(round(r, roundy)) if roundy!=0 else str(int(r))
        return '#splitline{[' + lstr + ';}{' + rstr + '[}'
    
    if old.InheritsFrom(TGraphAsymmErrors.Class()):
        h = TGraphAsymmErrors( old.GetN() )
        for ip in range(old.GetN()):
            h.SetPoint(ip, ip, old.GetPointY(ip) )
            h.SetPointError(ip, .5, .5,
                            old.GetErrorYlow(ip), old.GetErrorYhigh(ip) )

    elif old.InheritsFrom(TH2D.Class()): 
        name = old.GetName() + '_equal_width'
        nx, ny = old.GetNbinsX(), old.GetNbinsY()
        h = TH2D(name, name, nx, 0, nx, ny, 0, ny)
        for by in range(1, old.GetNbinsY()+1):
            for bx in range(1, old.GetNbinsX()+1):
                h.SetBinContent(bx, by, old.GetBinContent(bx, by))
                h.SetBinError(bx, by, old.GetBinError(bx, by))

        # Change bin labels
        for bx in range(1, old.GetNbinsX()+1):
            ledge = old.GetXaxis().GetBinLowEdge(bx)
            redge = old.GetXaxis().GetBinUpEdge(bx)
            h.GetXaxis().SetBinLabel(bx, xfunc(ledge,redge))
        for by in range(1, old.GetNbinsY()+1):
            dedge = old.GetYaxis().GetBinLowEdge(by)
            uedge = old.GetYaxis().GetBinUpEdge(by)
            h.GetYaxis().SetBinLabel(by, yfunc(dedge,uedge))
    else:
        mess = '[apply_equal_bin_width] '
        mess += 'The object should either be a TGraphasymmErrors or a TH2D'
        raise ValueError(mess)
    return h

def upify(s):
    """capitalizes the first letter of the passed string"""
    return s[0].upper() + s[1:]

def write_trigger_string(trig_comb, inters_str, items_per_line=1, join='+'):
    c1, c2 = 0, 0
    trig_str = trig_comb.split(inters_str)

    if len(trig_str)==1:
        return 'Trigger: ' + trig_str[0]
    
    loopstr = 'Triggers: '
    it = [iter(trig_str)]*items_per_line
    raw = list(zip(*it))
    raw = [' + '.join(x) for x in raw]
    rest = len(trig_str) % items_per_line
    end = trig_str[-rest:]
    
    if end != trig_str:
        raw += end

    for elem in raw:
        if elem == raw[-1]:
            loopstr += elem + '}'*(len(raw)-1)
        else:
            loopstr += '#splitline{'
            loopstr += elem + ' ' + join
            loopstr += '}{'

    return loopstr
