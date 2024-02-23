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

import inclusion
from inclusion.config import main
from inclusion.utils import utils

import ROOT

def add_slash(s):
    """Adds single slash to path if absent"""
    s = s if s[-1] == '/' else s + '/'
    return s

def add_vnames(*vnames):
    return '_VERSUS_'.join(vnames)

def at_least_two(x1, x2, x3):
    """Checks if at least two out of the three boleans are True."""
    return x1 if (x2 or x3) else (x2 and x3)

def build_script_command(name, sep, **kw):
    if name:
        p = build_script_path(name)
        comm = 'python3 {} '.format(p)
    else:
        comm = ' '
    for k,v in kw.items():
        if len(k)==1:
            comm += '-{} {}'.format(k,v)
        else:
            comm += '--{} {}'.format(k,v)
        comm += sep
    return comm
    
def build_script_path(name):
    path = os.path.join(main.local_folder,
                        main.folders['scripts'],
                        name)
    return path

def check_discr_vars_correctness(triggers, dset, channel, exclusive):
    """
    Check if the intersections exist, are complete, and do not have duplicates.
    Input dictionary should have the following signature:
    """
    chn_inters = generate_trigger_combinations(channel, triggers, exclusive)
    flatten = list(dset.keys())

    # type check
    for x in flatten:
        if not isinstance(x, tuple):
            mes = 'Intersection {} must be defined as a tuple'.format(x)
            raise TypeError(mes)
    
    # existence check
    for x in flatten:
        if x not in chn_inters and len(x)!=0:
            mes = 'Intersection {} is not required by channel {}.'
            mes = mes.format(f, channel)
            raise ValueError(mes)

    # completness check
    diff = set(chn_inters)-set(flatten)
    if diff:
        mes = 'Some intersections were not included in the configuration:\n'
        for elem in diff:
            mes += ' - ' + utils.join_name_trigger_intersection(elem) + '\n'
        raise ValueError(mes)

    # duplicate check
    dup = set([x for x in flatten if flatten.count(x) > 1])
    if dup:
        mes = 'Some intersections in the configuration are duplicated:\n'
        for elem in dup:
            mes += ' - ' + utils.join_name_trigger_intersection(elem) + '\n'
        raise ValueError(mes)
    
    return

def check_inters_correctness(triggers, dchn, dgen, channel, exclusive):
    """
    Check if the intersections exist, are complete, and do not have duplicates.
    Input dictionary should have the following signature:
    d = {'dataset1': (tuple1, tuple2,), 'dataset2': (tuple3, tuple4,), ...}
    """
    chn_inters = generate_trigger_combinations(channel, triggers, exclusive)
    flatten =  [ tuple(sorted(w)) for x in dchn for w in dchn[x] ]
    flatten += [ tuple(sorted(w)) for x in dgen for w in dgen[x] ]
    
    # type check
    for x in flatten:
        if not isinstance(x, tuple):
            mes = 'Intersection {} must be defined as a tuple'.format(x)
            raise TypeError(mes)

    # existence check
    for x in flatten:
        if x not in chn_inters and len(x)!=0:
            mes = 'Intersection {} is not required by channel {}.'
            mes = mes.format(x, channel)
            raise ValueError(mes)

    # completness check
    diff = set(chn_inters)-set(flatten)
    if diff:
        mes = 'Some intersections were not included in the configuration:\n'
        for elem in diff:
            mes += ' - ' + utils.join_name_trigger_intersection(elem) + '\n'
        raise ValueError(mes)

    # duplicate check
    dup = set([x for x in flatten if flatten.count(x) > 1])
    if dup:
        mes = 'Some intersections in the configuration are duplicated:\n'
        for elem in dup:
            mes += ' - ' + utils.join_name_trigger_intersection(elem) + '\n'
        raise ValueError(mes)
    
    return

def check_triggers_exist(l, triggers):
    """Check individual triggers match the ones in the configuration file."""
    for elem in l:
        if elem not in triggers:
            mess = '[utils.check_triggers_exist] '
            mess += 'Trigger {} is not supported'.format(elem)
            raise ValueError(mess)
    return
    
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

def debug(message, flag=True, fname=None):
    decorator = ' ============ '
    if flag:
        mes = '[' + os.path.basename(fname) + '] ' if fname is not None else ''
        mes += decorator + message + decorator
        print(mes, flush=True)

def define_used_tree_variables(cut):
    """
    Return the list of all required TTree variables.
    This is the sum of default ones plus the ones contained
    in the user-provided custom cut.
    Repeated variables are deleted.
    """
    _entries = ('triggerbit', 'RunNumber', 'HHKin_mass', 'isLeptrigger', 'pairType', 'isOS',
                'MC_weight', 'IdSF_deep_2d', 'PUReweight', 'L1pref_weight', 'trigSF', 'PUjetID_SF', 'bTagweightReshape',
                'dau1_eleMVAiso', 'dau1_iso', 'dau2_iso', 'dau1_deepTauVsJet', 'dau2_deepTauVsJet',
                'nleps', 'nbjetscand', 'tauH_mass', 'bH_mass_raw',
                'bjet1_bID_deepFlavor', 'bjet2_bID_deepFlavor',
                'isVBF', 'VBFjj_mass', 'VBFjj_deltaEta')
    if cut is not None:
        _regex = tuple(set(re.findall(r'self\.entries\.(.+?)\s', cut)))
    else:
        _regex = ()
    return tuple(set(_entries + _regex))
    
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

def flatten_nested_dict(d):
    """
    Splits keys and values.
    Dictionaries are not straightforward to pass as arguments.
    """
    keys, vals = ([] for _ in range(2))
    for k,v in d.items():
        for x in v:
            keys.append(k)
        vals.extend(v)
    return keys, vals

def generate_trigger_combinations(channel, trigs, exclusive):
    """
    Set all possible trigger intersection combinations, per channel.
    Does not consider intersections between incompatible triggers and channels
    (for instance, `etau` channel with `IsoMu24` trigger).
    Each intersection is sorted alphabetically
    (useful for matching with the KLUB framework).
    """
    # look only at combinations where the channel is imcompatible with the trigger
    pruntrigs = [ exclusive[x] for x in exclusive if x not in (channel, 'general') ]
    pruntrigs = set(it.chain(*pruntrigs))
    pruntrigs = set(trigs) - pruntrigs

    length1 = list(it.chain.from_iterable(it.combinations(sorted(pruntrigs), 1)))
    check_triggers_exist(length1, trigs)
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

    if var == "metnomu_et":
        var_custom = r"MET (no-\mu)"
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
    ph = main.placeholder_cuts
    if opt == 'Ref1D':
        return lambda a,b,c : 'Ref1D_{}_{}_{}'.format(a,b,c)
    elif opt == 'Trig1D':
        return lambda a,b,c : 'Trig1D_{}_{}_{}{}'.format(a,b,c,ph)
    elif opt == 'Ref2D':
        return lambda a,b,c : 'Ref2D_{}_{}_{}'.format(a,b,c)
    elif opt == 'Trig2D':
        return lambda a,b,c : 'Trig2D_{}_{}_{}{}'.format(a,b,c,ph)
    elif opt == 'Canvas2D':
        return lambda a,b,c : 'Canvas2D_{}_{}_{}{}'.format(a,b,c,ph)
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

def get_root_inputs(proc, indir, include_tree=False):
    #### Check input folder
    if not isinstance(indir, (tuple,list)):
        indir = [indir]
    fexists = []
    for idir in indir:
        g = glob.glob(os.path.join(idir, proc + '*/goodfiles.txt'))
        fexists.append(False if len(g)==0 else True)

    # check one and only one is True
    if sum(fexists) == 0:
        m = '[get_root_input_files] No file was found for the {} sample. '.format(proc)
        m += 'The inspect path was {}'.format('INDIR/' + proc + '*/goodfiles.txt')
        m += ' where INDIR stands for the following: {}.'.format(indir)
        raise ValueError(m)
    elif sum(fexists) > 1:
        m = '[get_root_input_files] More than one file exists for the {} sample.'.format(proc)
        m += ' Selecting directory {} '.format(indir[fexists.index(True)])
        m += 'from the following list: {}.'.format(indir)
        raise ValueError(m)

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
                    if include_tree:
                        line = line[:-1] # remove new line
                        if not os.path.exists(line): 
                            mes = '[' + os.path.basename(__file__) + '] '
                            mes += ' The input file does not exist: {}'.format(line)
                            raise ValueError(mes)
                        filelist.append(line + ':HTauTauTree')
                    else:
                        if line[:-1] not in main.corrupted_files:
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

def hadd_subpaths(args, channel=''):
    channel_str = '' if channel=='' else channel+'_'
    _tbase1 = args.tprefix + channel_str + args.dataset_name
    _tbase2 = '_Sum' + args.subtag
    return _tbase1, _tbase2

def is_channel_consistent(chn, pairtype):
    opdict = { '<':  operator.lt,
               '>':  operator.gt,
               '==': operator.eq }

    op, val = main.sel[chn]['pairType']
    return opdict[op](pairtype, val)

def is_trigger_comb_in_channel(chn, tcomb, triggers, exclusive):
    possible_trigs = generate_trigger_combinations(chn, triggers, exclusive)
    split = tuple(tcomb.split(main.inters_str))
    return split in possible_trigs
                               
def is_nan(num):
    return num!= num

def join_name_trigger_intersection(tuple_element):
    inters = main.inters_str
    return inters.join(tuple_element)

def join_strings(*args, sep=''):
    return sep.join(args)

def key_exists(d, k1, k2):
    return k1 in d and k2 in d[k1]

def load_binning(afile, key, variables, channels):
    """
    Load the Binning stored in a previous task.
    """
    binedges, nbins = ({} for _ in range(2))
    with h5py.File(afile, 'r') as f:
        try:
            group = f[key]
        except KeyError:
          missing_key_print(afile, key)
          
        for var in variables:
            try:
                subgroup = group[var]
            except KeyError:
                missing_key_print(afile, key)
                
            binedges[var], nbins[var] = ({} for _ in range(2))
            for chn in channels:
                binedges[var][chn] = np.array(subgroup[chn][:])
                nbins[var][chn] = len(binedges[var][chn]) - 1

    return binedges, nbins
    
def missing_key_print(afile, key):
    print('{} does not have key {}.'.format(afile, key))
    print('Available keys: {}'.format(f.keys()))
    raise

def redraw_border():
    """
    this little macro redraws the axis tick marks and the pad border lines.
    """
    ROOT.gPad.Update()
    ROOT.gPad.RedrawAxis()
    l = ROOT.TLine()
    l.SetLineWidth(2)

    l.DrawLine(ROOT.gPad.GetUxmin(), ROOT.gPad.GetUymax(), ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymax()) #top border
    l.DrawLine(ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymin(), ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymax()) #right border
    l.DrawLine(ROOT.gPad.GetUxmin(), ROOT.gPad.GetUymin(), ROOT.gPad.GetUxmin(), ROOT.gPad.GetUymax()) #left border
    l.DrawLine(ROOT.gPad.GetUxmin(), ROOT.gPad.GetUymin(), ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymin()) #bottom border

def remove(f):
    if os.path.exists(f):
        os.remove(f)

def rewrite_cut_string(oldstr, newstr, regex=False):
    if regex:
        _regex = re.findall(r'^.*CUTS_(.+)$', newstr)
        assert(len(_regex)==1)
        _regex = _regex[0]
        newstr = _regex.replace('>', 'L').replace('<', 'S').replace('.', 'p')
    
    res = oldstr.replace(main.placeholder_cuts, '_CUTS_'+newstr)
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

def split_vnames(joinvars):
    return joinvars.split('_VERSUS_')

def stoi(s):
    return int(float(s))

def parse_args(parser):
    args = parser.parse_args()
    print_configuration(args)
    return args

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
            v = ' '.join((str(x) for x in v))
        print('{0:>{d1}}   {1}'.format(k, v, d1=maxlkey+3), flush=True)
    print('----------------------------------------', flush=True)

def get_ptcuts(channel, year):
    """
    Trigger pT cuts.
    Order: single lepton, cross lepton leg, cross tau leg.
    """
    if channel == "etau":
        if year == "2016":
            ptcuts = (26,)
        elif year == "2017" or year == "2018":
            ptcuts = (33, 25, 35)

    elif channel == "mutau":
        if year == "2016":
            ptcuts = (25, 20, 25)
        elif year == "2017":
            ptcuts = (28, 21, 32)
        elif year == "2018":
            ptcuts = (25, 21, 32)

    elif channel == "tautau":
        ptcuts = (40, 40)

    return ptcuts
    
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
    
    if old.InheritsFrom(ROOT.TGraphAsymmErrors.Class()):
        h = ROOT.TGraphAsymmErrors( old.GetN() )
        for ip in range(old.GetN()):
            h.SetPoint(ip, ip, old.GetPointY(ip) )
            h.SetPointError(ip, .5, .5,
                            old.GetErrorYlow(ip), old.GetErrorYhigh(ip) )

    elif old.InheritsFrom(ROOT.TH2D.Class()): 
        name = old.GetName() + '_equal_width'
        nx, ny = old.GetNbinsX(), old.GetNbinsY()
        h = ROOT.TH2D(name, name, nx, 0, nx, ny, 0, ny)
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

def total_sum_weights(f, isdata):
    if isdata:
        return 1.
    else:
        search_str = os.path.join(os.path.dirname(f), 'goodfiles.txt')
        xsec_norm = 0.
        with open(search_str, 'r') as afile:
            for elem in afile:
                ftmp = ROOT.TFile(elem.replace('\n', ''), "READ")
                xsec_norm += ftmp.Get('h_eff').GetBinContent(1)     
        return xsec_norm
    return None

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
