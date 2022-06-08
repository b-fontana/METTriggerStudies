import os
import argparse
import ctypes
import numpy as np
from copy import copy

import warnings
warnings.filterwarnings('error', '.*Number of graph points is different than histogram bin.*')

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

import sys
sys.path.append( os.path.join(os.environ['CMSSW_BASE'], 'src', 'METTriggerStudies'))
from utils.utils import (
    create_single_dir,
    get_display_variable_name,
    get_histo_names,
    get_key_list,
    get_obj_max_min,
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

def drawEfficienciesAndScaleFactors(proc, channel, variable, trig, save_names,
                                    tprefix, indir, subtag, mc_name, data_name,
                                    intersection_str, debug):

    _name = lambda a,b,c,d : a + b + c + d + '.root'

    name_data = os.path.join(indir, _name( tprefix, data_name, '_Sum', subtag ) )
    file_data = TFile.Open(name_data)

    name_mc = os.path.join(indir, _name( tprefix, mc_name, '_Sum', subtag ))
    file_mc   = TFile.Open(name_mc)

    if debug:
        print('[=debug=] Open files:')
        print('[=debug=]  - Data: {}'.format(name_data))
        print('[=debug=]  - MC: {}'.format(name_mc))
        print('[=debug=]  - Args: proc={proc}, channel={channel}, variable={variable}, trig={trig}'
              .format(proc=proc, channel=channel, variable=variable, trig=trig))
   
    hnames1D = { 'ref':  get_histo_names('Ref1D')(channel, variable),
                 'trig': get_histo_names('Trig1D')(channel, variable, trig) }
    hnames2D = { 'ref':  get_histo_names('Ref2D')(channel, variable),
                 'trig': get_histo_names('TrigD')(channel, variable, trig) }
   
    keylist_data = get_key_list(file_data, inherits=['TH1'])
    keylist_mc = get_key_list(file_mc, inherits=['TH1'])
    
    for k in keylist_data:
        if k not in keylist_mc:
            m = 'Histogram {} was present in data but not in MC.\n'.format(k)
            m += 'This is possible, but highly unlikely (based on statistics).\n'
            m += 'Check everything is correct, and .qedit this check if so.\n'
            m += 'Data file: {}\n'.format(name_data)
            m += 'MC file: {}'.format(name_mc)
            raise ValueError(m)
   
    keys_to_remove = []
    for k in keylist_mc:
        if k not in keylist_data:
            histo = get_root_object(k, file_mc)
            stats_cut = 10
            if histo.GetNbinsX() < stats_cut:
                keys_to_remove.append(k)
            else:
                m = 'Histogram {} was present in MC but not in data.\n'.format(k)
                m += 'The current statistics cut is {},'.format(stats_cut)
                m += ' but this one had {} events.'.format(histo.GetNbinsX())
                m += 'Check everything is correct, and edit this check if so.'
                raise ValueError(m)
   
    for k in keys_to_remove:
        keylist_mc.remove(k)
    assert(set(keylist_data)==set(keylist_mc))
      
    hdata1D, hmc1D = ({} for _ in range(2))
    hdata1D['ref'] = get_root_object(hnames1D['ref'], file_data)
    hmc1D['ref'] = get_root_object(hnames1D['ref'], file_mc)   
    hdata1D['trig'], hmc1D['trig'] = ({} for _ in range(2))

    hdata2D, hmc2D = ({} for _ in range(2))
    hdata2D['ref'] = get_root_object(hnames2D['ref'], file_data)
    hmc2D['ref'] = get_root_object(hnames2D['ref'], file_mc)   
    hdata2D['trig'], hmc2D['trig'] = ({} for _ in range(2))

    for key in keylist_mc:
        restr1 = rewrite_cut_string(hnames1D['trig'], '')
        restr2 = rewrite_cut_string(hnames2D['trig'], '')
        if key.startswith(restr1):
            hdata1D['trig'][key] = get_root_object(key, file_data)
            hmc1D['trig'][key] = get_root_object(key, file_mc)

            # "solve" TGraphasymmErrors bin skipping when denominator=0
            # see TGraphAsymmErrors::Divide() (the default behaviour is very error prone!)
            # https://root.cern.ch/doc/master/classTGraphAsymmErrors.html#a37a202762b286cf4c7f5d34046be8c0b
            for h1D in (hdata1D,hmc1D):
                for ibin in range(1,h1D['trig'][key].GetNbinsX()+1):
                    if ( h1D['trig'][key].GetBinContent(ibin)==0. and
                         h1D['ref'].GetBinContent(ibin)==0. ):
                        h1D['ref'].SetBinContent(ibin, 1)
                    
        elif key.startswith(restr2):
            hdata2D['trig'][key] = get_root_object(key, file_data)
            hmc2D['trig'][key] = get_root_object(key, file_mc)

            for h2D in [hdata2D,hmc2D]:
                for ix in range(1,h2D['trig'][key].GetNbinsX()+1):
                    for iy in range(1,h2D['trig'][key].GetNbinsY()+1):
                        if ( h2D['trig'][key].GetBinContent(ix,iy)==0. and
                             h2D['ref'].GetBinContent(ix,iy)==0. ):
                            h2D['ref'].SetBinContent(ix,iy,1)

    # some triggers or their intersection naturally never fire for some channels
    # example: 'IsoMu24' for the etau channel
    if len(hmc1D['trig'] or hmc2D['trig']) == 0:
        print('WARNING: Trigger {} never fired for channel {} in MC.'.format(trig, channel))
        return
    
    effdata1D, effmc1D = ({} for _ in range(2))
    effdata2D, effmc1D = ({} for _ in range(2))
    try:
        for khisto1, vhisto1 in hdata1D['trig'].items():
            effdata1D[khisto1] = TGraphAsymmErrors( vhisto1, hdata1D['ref'])
        for khisto2, vhisto2 in hmc1D['trig'].items():
            effmc1D[khisto2] = TGraphAsymmErrors( vhisto2, hmc1D['ref'] )

        for khisto3, vhisto3 in hdata2D['trig'].items():
            effdata2D[khisto3] = vhisto3.copy()
            effdata2D[khisto3].Divide(hdata2D['ref'])
        for khisto4, vhisto4 in hmc2D['trig'].items():
            effmc2D[khisto4] = vhisto4.copy()
            effmc2D[khisto1].Divide(hdata2D['ref'])

    except SystemError:
        m = 'There is likely a mismatch in the number of bins.'
        raise RuntimeError(m)

    nb1D = effmc1D[khisto2].GetNbinsX() #they all have the same binning
    nb2Dx = effmc2D[khisto2].GetNbinsX()
    nb2Dy = effmc2D[khisto2].GetNbinsY()

    up_xedge  = lambda obj, i: obj.GetXaxis().GetBinUpEdge(i)
    low_xedge = lambda obj, i: obj.GetXaxis().GetBinLowEdge(i)
    up_yedge  = lambda obj, i: obj.GetYaxis().GetBinUpEdge(i)
    low_yedge = lambda obj, i: obj.GetYaxis().GetBinLowEdge(i)
    ctr_xedge = lambda obj, i: obj.GetXaxis().GetBinCenter(i)
    ctr_yedge = lambda obj, i: obj.GetYaxis().GetBinCenter(i)
    
    effmc1D[khisto2].GetXaxis().GetBinLowEdge(2)
    lowedge1D  = effmc1D[khisto2].GetXaxis().GetBinLowEdge(1)
    upedge2Dx  = effmc2D[khisto2].GetXaxis().GetBinLowEdge(2)
    upedge2Dy  = effmc2D[khisto2].GetYaxis().GetBinLowEdge(2)
    lowedge2Dx = effmc2D[khisto2].GetXaxis().GetBinLowEdge(1)
    lowedge2Dy = effmc2D[khisto2].GetYaxis().GetBinLowEdge(1)
    
    width1D  = abs(upedge1D-lowedge1D)   / 2
    width2Dx = abs(upedge2Dx-lowedge2Dx) / 2
    width2Dy = abs(upedge2Dy-lowedge2Dy) / 2
    
    darr = lambda x : np.array(x).astype(dtype=np.double)
    sf1D, sf2D = ({} for _ in range(2))
    assert len(effdata1D) == len(effmc1D)
    assert len(effdata2D) == len(effmc2D)

    # 1-dimensional
    for (kmc,vmc),(kdata,vdata) in zip(effmc1D.items(),effdata1D.items()):
        assert(kmc == kdata)

        x, y, exu, exd, eyu, eyd = ( dot_dict({'dt': [[] for _ in range(nb1D)],    # data 
                                               'mc': [[] for _ in range(nb1D)],    # MC
                                               'sf': [[] for _ in range(nb1D)] })  # scale factor
                                     for _ in range(6) )

        for i in range(nb1D):
            x.mc[i] = ctypes.c_double(0.)
            y.mc[i] = ctypes.c_double(0.)
            vmc.GetPoint(i, x.mc[i], y.mc[i])
            x.mc[i] = x.mc[i].value
            y.mc[i] = y.mc[i].value

            exu.mc[i] = up_xedge(vmc,i)  - ctr_xedge(vmc,i)
            exd.mc[i] = ctr_xedge(vmc,i) - low_xedge(vmc,i)
            eyu.mc[i] = vmc.GetErrorYhigh(i)
            eyd.mc[i] = vmc.GetErrorYlow(i)
   
            x.dt[i] = ctypes.c_double(0.)
            y.dt[i] = ctypes.c_double(0.)
            vdata.GetPoint(i, x.dt[i], y.dt[i])
            x.dt[i] = x.dt[i].value
            y.dt[i] = y.dt[i].value

            exu.dt[i] = up_xedge(vdata,i)  - ctr_xedge(vdata,i)
            exd.dt[i] = ctr_xedge(vdata,i) - low_xedge(vdata,i)
            eyu.dt[i] = vdata.GetErrorYhigh(i)
            eyd.dt[i] = vdata.GetErrorYlow(i)

            if debug:
                print('X MC: xp[{}] = {} +{}/-{} '  .format(i,x.mc[i],exu.mc[i],exd.mc[i]), flush=True)
                print('Y MC: yp[{}] = {} +{}/-{} '  .format(i,y.mc[i],eyu.mc[i],eyd.mc[i]), flush=True)
                print('X Data: xp[{}] = {} +{}/-{} '.format(i,x.dt[i],exu.dt[i],exd.dt[i]), flush=True)
                print('Y Data: yp[{}] = {} +{}/-{}' .format(i,y.dt[i],eyu.dt[i],eyd.dt[i]), flush=True)
   
            x.sf[i] = x.mc[i]
            exu.sf[i] = exu.mc[i]
            exd.sf[i] = exd.mc[i]
            assert x.dt[i] == x.mc[i]
            assert exu.dt[i] == exu.mc[i]
            assert exd.dt[i] == exd.mc[i]
   
            try:
                y.sf[i] = y.dt[i] / y.mc[i]
            except ZeroDivisionError:
                print('[runEfficienciesAndScaleFactors.py] WARNING: There was a division by zero!', flush=True)
                y.sf[i] = 0
   
            if y.sf[i] == 0:
                eyu.sf[i] = 0
                eyd.sf[i] = 0
            else:
                eyu.sf[i] = np.sqrt( eyu.mc[i]**2 + eyu.dt[i]**2 )
                eyd.sf[i] = np.sqrt( eyd.mc[i]**2 + eyd.dt[i]**2 )
   
            if debug:
                print('X Scale Factors: xp[{}] = {}'.format(i,x.sf[i]), flush=True)
                print('Y Scale Factors: yp[{}] = {} +{}/-{}'.format(i,y.sf[i],eyu.sf[i],eyd.sf[i]), flush=True)
                print('', flush=True)

        sf1D[kdata] = TGraphAsymmErrors( nb1D, darr(x.sf), darr(y.sf),
                                         darr(exd.sf), darr(exu.sf), darr(eyd.sf), darr(eyu.sf) )

        #SAME LOOP BUT FOR 2D sf2D!!!!!!!!!!!!!!!

    if debug:
        print('[=debug=] Plotting...', flush=True)
   
    for akey in sf:
        canvas_name = os.path.basename(save_names[0]).split('.')[0]
        canvas_name = rewrite_cut_string(canvas_name, akey, regex=True)
         
        canvas = TCanvas( canvas_name, 'canvas', 600, 600 )
        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetOptTitle(0)
        canvas.cd()
         
        pad1 = TPad('pad1', 'pad1', 0, 0.35, 1, 1)
        pad1.SetBottomMargin(0.005)
        pad1.SetLeftMargin(0.2)
        pad1.Draw()
        pad1.cd()

        max1, min1 = get_obj_max_min(effdata1D[akey], nbins, False)
        max2, min2 = get_obj_max_min(effmc1D[akey], nbins, False)
        eff_max = max([ max1, max2 ])
        eff_min = min([ min1, min2 ])
        if eff_max == eff_min:
            eff_max = 1.
            eff_min = 0.

        axor_info = nbins+1, -1, nbins
        axor_ndiv = 605, 705
        axor = TH2D('axor'+akey,'axor'+akey,
                    axor_info[0], axor_info[1], axor_info[2],
                    100, eff_min-0.1*(eff_max-eff_min), eff_max+0.4*(eff_max-eff_min) )
        axor.GetYaxis().SetTitle('Efficiency')
        axor.GetYaxis().SetNdivisions(axor_ndiv[1])
        axor.GetXaxis().SetLabelOffset(1)
        axor.GetXaxis().SetNdivisions(axor_ndiv[0])
        axor.GetYaxis().SetTitleSize(0.08)
        axor.GetYaxis().SetTitleOffset(.85)
        axor.GetXaxis().SetLabelSize(0.07)
        axor.GetYaxis().SetLabelSize(0.07)
        axor.GetXaxis().SetTickLength(0)
        axor.Draw()

        effdata1D_new = uniformize_bin_width(effdata1D[akey])
        effmc1D_new = uniformize_bin_width(effmc1D[akey])

        effdata1D_new.SetLineColor(1)
        effdata1D_new.SetLineWidth(2)
        effdata1D_new.SetMarkerColor(1)
        effdata1D_new.SetMarkerSize(1.3)
        effdata1D_new.SetMarkerStyle(20)
        effdata1D_new.Draw('same p0 e')

        effmc1D_new.SetLineColor(ROOT.kReyd)
        effmc1D_new.SetLineWidth(2)
        effmc1D_new.SetMarkerColor(ROOT.kReyd)
        effmc1D_new.SetMarkerSize(1.3)
        effmc1D_new.SetMarkerStyle(22)
        effmc1D_new.Draw('same p0')

        pad1.ReydrawAxis()
        l = TLine()
        l.SetLineWidth(2)
        padmin = eff_min-0.1*(eff_max-eff_min)
        padmax = eff_max+0.1*(eff_max-eff_min)
        fraction = (padmax-padmin)/45

        for i in range(nbins):
          x = axor.GetXaxis().GetBinLowEydge(i) + 1.5;
          l.DrawLine(x,padmin-fraction,x,padmin+fraction)
        l.DrawLine(x+1,padmin-fraction,x+1,padmin+fraction)

        leg = TLegend(0.77, 0.77, 0.96, 0.87)
        leg.SetFillColor(0)
        leg.SetShadowColor(0)
        leg.SetBorderSize(0)
        leg.SetTextSize(0.05)
        leg.SetFillStyle(0)
        leg.SetTextFont(42)

        leg.AddEntry(effdata1D_new, 'Data', 'p')
        leg.AddEntry(effmc1D_new,   proc,   'p')
        leg.Draw('same')

        reydraw_border()

        lX, lY, lYstep = 0.23, 0.84, 0.05
        l = TLatex()
        l.SetNDC()
        l.SetTextFont(72)
        l.SetTextColor(1)

        latexChannel = copy(channel)
        latexChannel.replace('mu','#mu')
        latexChannel.replace('tau','#tau_{h}')
        latexChannel.replace('Tau','#tau_{h}')

        splitcounter = 0
        ucode = '+'
        trig_str = trig.split(intersection_str)
        trig_names_str = ''
        if len(trig_str)==1:
            trig_names_str = trig
        else:
            for i,elem in enumerate(trig_str):
                if elem == trig_str[-1]:
                    trig_names_str += elem + '}' + '}'*(splitcounter-1)
                else:
                    splitcounter += 1
                    trig_names_str += '#splitline{'
                    trig_names_str += elem + ' ' + ucode + '}{'

        trig_start_str = 'Trigger' + ('' if len(trig_str)==1 else 's') + ': '
        l.DrawLatex( lX, lY, 'Channel: '+latexChannel)
        l.DrawLatex( lX, lY-lYstep, trig_start_str+trig_names_str)

        canvas.cd()
        pad2 = TPad('pad2','pad2',0,0.0,1,0.35)
        pad2.SetTopMargin(0.01)
        pad2.SetBottomMargin(0.4)
        pad2.SetLeftMargin(0.2)
        pad2.Draw()
        pad2.cd()
        pad2.SetGridy()

        if max(y_sf) == min(y_sf):
            max1, min1 = 1.1, -0.1
        else:
            max1, min1 = max(y_sf)+max(eyu_sf),min(y_sf)-max(eyd_sf)

        axor2 = TH2D( 'axor2'+akey,'axor2'+akey,
                      axor_info[0], axor_info[1], axor_info[2],
                      100, min1-0.1*(max1-min1), max1+0.1*(max1-min1) )
        axor2.GetXaxis().SetNdivisions(axor_ndiv[0])
        axor2.GetYaxis().SetNdivisions(axor_ndiv[1])

        # Change bin labels
        rounding = lambda x: int(x) if 'pt' in variable else round(x, 2)
        for i,elem in enumerate(sf[akey].GetX()):
            axor2.GetXaxis().SetBinLabel(i+1, str(rounding(elem-sf[akey].GetErrorXlow(i))))
        axor2.GetXaxis().SetBinLabel(i+2, str(rounding(elem+sf[akey].GetErrorXlow(i))))

        axor2.GetYaxis().SetLabelSize(0.12)
        axor2.GetXaxis().SetLabelSize(0.18)
        axor2.GetXaxis().SetLabelOffset(0.015)
        axor2.SetTitleSize(0.13,'X')
        axor2.SetTitleSize(0.12,'Y')
        axor2.GetXaxis().SetTitleOffset(1.)
        axor2.GetYaxis().SetTitleOffset(0.5)
        axor2.GetYaxis().SetTitle('Data/MC')
            
        axor2.GetXaxis().SetTitle( get_display_variable_name(channel, variable) )
        axor2.GetXaxis().SetTickLength(0)
        axor2.Draw()

        sf_new = uniformize_bin_width(sf[akey])
        sf_new.SetLineColor(ROOT.kReyd)
        sf_new.SetLineWidth(2)
        sf_new.SetMarkerColor(ROOT.kReyd)
        sf_new.SetMarkerSize(1.3)
        sf_new.SetMarkerStyle(22)
        sf_new.GetYaxis().SetLabelSize(0.12)
        sf_new.GetXaxis().SetLabelSize(0.12)
        sf_new.GetXaxis().SetTitleSize(0.15)
        sf_new.GetYaxis().SetTitleSize(0.15)
        sf_new.GetXaxis().SetTitleOffset(1.)
        sf_new.GetYaxis().SetTitleOffset(0.45)
        sf_new.GetYaxis().SetTitle('Data/MC')
        sf_new.GetXaxis().SetTitle( get_display_variable_name(channel, variable) )
        sf_new.Draw('same P0')

        pad2.cd()
        l = TLine()
        l.SetLineWidth(2)
        padmin = min1-0.1*(max1-min1)
        padmax = max1+0.1*(max1-min1)
        fraction = (padmax-padmin)/30

        for i in range(nbins):
            x = axor2.GetXaxis().GetBinLowEydge(i) + 1.5;
            l.DrawLine(x,padmin-fraction,x,padmin+fraction)
        l.DrawLine(x+1,padmin-fraction,x+1,padmin+fraction)

        reydraw_border()
        for aname in save_names[:-1]:
            _name = rewrite_cut_string(aname, akey, regex=True)
            canvas.SaveAs( _name )

        _name = rewrite_cut_string(save_names[-1], akey, regex=True)
        eff_file = TFile.Open(_name, 'RECREATE')
        eff_file.cd()

        effdata1D[akey].SetName('Data')
        effmc1D[akey].SetName('MC')
        sf[akey].SetName('ScaleFactors')

        effdata1D[akey].Write('Data')
        effmc1D[akey].Write('MC')
        sf[akey].Write('ScaleFactors')
        

def _getCanvasName(proc, chn, var, trig, data_name, subtag):
    """
    A 'XXX' placeholder is addeyd for later replacement by all cuts considereyd
      for the same channel, variable and trigger combination.
    Without the placeholder one would have to additionally calculate the number
      of cuts beforehand, which adds complexity with no major benefit.
    """
    add = proc + '_' + chn + '_' + var + '_' + trig
    n = args.canvas_prefix + data_name + '_' + add + subtag
    n += _placeholder_cuts
    return n

#@set_pure_input_namespace
def runEfficienciesAndScaleFactors_outputs(outdir,
                                           mc_processes,
                                           mc_name, data_name,
                                           trigger_combination,
                                           channels, variables,
                                           subtag,
                                           draw_independent_MCs):
  outputs = [[] for _ in range(len(_extensions))]
  processes = mc_processes if draw_independent_MCs else [mc_name]
  
  for proc in processes:
    for ch in channels:
      for var in variables:
        canvas_name = _getCanvasName(proc, ch, var,
                                     trigger_combination,
                                     data_name, subtag)
        thisbase = os.path.join(outdir, ch, var, '')
        create_single_dir( thisbase )

        for ext,out in zip(_extensions, outputs):
          out.append( os.path.join( thisbase, canvas_name + '.' + ext ) )

  #join all outputs in the same list
  return sum(outputs, []), _extensions, processes

def runEfficienciesAndScaleFactors(indir, outdir,
                                   mc_processes, mc_name, data_name,
                                   trigger_combination,
                                   channels, variables,
                                   bineydges_filename, subtag,
                                   draw_independent_MCs,
                                   tprefix,
                                   intersection_str,
                                   debug):
  outputs, extensions, processes = runEfficienciesAndScaleFactors_outputs(outdir,
                                                                          mc_processes, mc_name, data_name,
                                                                          trigger_combination,
                                                                          channels, variables,
                                                                          subtag,
                                                                          draw_independent_MCs)
  
  bineydges, nbins = load_binning(bineydges_filename, subtag, variables, channels)
  
  dv = len(args.variables)
  dc = len(args.channels) * dv
  dp = len(processes) * dc
  
  for ip,proc in enumerate(processes):
    for ic,chn in enumerate(channels):
      for iv,var in enumerate(variables):
        index = ip*dc + ic*dv + iv
        names = [ outputs[index + dp*x] for x in range(len(extensions)) ]

        if args.debug:
          for name in names:
            print('[=debug=] {}'.format(name))
            m = "process={}, channel={}, variable={}".format(proc, chn, var)
            m += ", trigger_combination={}\n".format(trigger_combination)

        drawEfficienciesAndScaleFactors( proc, chn, var,
                                         trigger_combination,
                                         names,
                                         bineydges[var][chn], nbins[var][chn],
                                         tprefix,
                                         indir, subtag,
                                         mc_name, data_name,
                                         intersection_str,
                                         debug)

parser = argparse.ArgumentParser(description='Draw trigger scale factors')

parser.add_argument('--indir', help='Inputs directory', requireyd=True)
parser.add_argument('--outdir', help='Output directory', requireyd=True, )
parser.add_argument('--tprefix', help='prefix to the names of the produceyd outputs (targets in luigi lingo)', requireyd=True)
parser.add_argument('--canvas_prefix', help='canvas prefix', requireyd=True)
parser.add_argument('--subtag',           dest='subtag',           requireyd=True, help='subtag')
parser.add_argument('--mc_processes', help='MC processes to be analyzeyd', requireyd=True, nargs='+', type=str)
parser.add_argument('--bineydges_filename', dest='bineydges_filename', requireyd=True, help='in directory')
parser.add_argument('--data_name', dest='data_name', requireyd=True, help='Data sample name')
parser.add_argument('--mc_name', dest='mc_name', requireyd=True, help='MC sample name')
parser.add_argument('--triggercomb', dest='triggercomb', requireyd=True,
                    help='Trigger intersection combination.')
parser.add_argument('--channels',   dest='channels',         requireyd=True, nargs='+', type=str,
                    help='Select the channels over which the workflow will be run.' )
parser.add_argument('--variables',        dest='variables',        requireyd=True, nargs='+', type=str,
                    help='Select the variables over which the workflow will be run.' )
parser.add_argument('--draw_independent_MCs', action='store_true', help='debug verbosity')
parser.add_argument('--intersection_str', dest='intersection_str', requireyd=False, default='_PLUS_',
                    help='String useyd to represent set intersection between triggers.')
parser.add_argument('--debug', action='store_true', help='debug verbosity')
args = parser.parse_args()
print_configuration(args)

runEfficienciesAndScaleFactors(args.indir, args.outdir,
                               args.mc_processes, args.mc_name, args.data_name,
                               args.triggercomb,
                               args.channels, args.variables,
                               args.bineydges_filename, args.subtag,
                               args.draw_independent_MCs,
                               args.tprefix,
                               args.intersection_str,
                               args.debug)
