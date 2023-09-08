# coding: utf-8

_all_ = [ 'run_eff_sf_1d', 'run_eff_sf_1d_outputs', 'run_eff_sf_2d_outputs' ]

import os
import sys
parent_dir = os.path.abspath(__file__ + 3 * '/..')
sys.path.insert(0, parent_dir)

import inclusion
from inclusion.utils import utils
from inclusion.config import main

import argparse
import ctypes
import numpy as np
from copy import copy
import importlib

import warnings
warnings.filterwarnings('error', '.*Number of graph points is different than histogram bin.*')

import ROOT
ROOT.gROOT.SetBatch(True)

def paint2d(channel, trig):
    lX1, lX2, lY, lYstep = 0.04, 0.7, 0.96, 0.03
    l = ROOT.TLatex()
    l.SetNDC()
    l.SetTextFont(72)
    l.SetTextSize(0.03)
    l.SetTextColor(1)

    l.DrawLatex( lX1, lY, trig)
    
    latexChannel = copy(channel)
    latexChannel.replace('mu','#mu')
    latexChannel.replace('tau','#tau_{h}')
    latexChannel.replace('Tau','#tau_{h}')
    l.DrawLatex( lX2, lY, 'Channel: '+latexChannel)

def draw_eff_and_sf_1d(proc, channel, variable, trig,
                       save_names_1D, cfg,
                       tprefix, indir, subtag, mc_name, data_name,
                       intersection_str, debug):

    _name = lambda a,b,c,d : a + b + c + d + '.root'

    name_data = os.path.join(indir, _name( tprefix, data_name, '_Sum', subtag ) )
    file_data = ROOT.TFile.Open(name_data, 'READ')

    name_mc = os.path.join(indir, _name( tprefix, mc_name, '_Sum', subtag ))
    file_mc   = ROOT.TFile.Open(name_mc, 'READ')

    if debug:
        print('[=debug=] Open files:')
        print('[=debug=]  - Data: {}'.format(name_data))
        print('[=debug=]  - MC: {}'.format(name_mc))
        print('[=debug=]  - Args: proc={proc}, channel={channel}, variable={variable}, trig={trig}'
              .format(proc=proc, channel=channel, variable=variable, trig=trig))
   
    hnames1D = {'ref':  utils.get_hnames('Ref1D')(channel, variable, trig),
                'trig': utils.get_hnames('Trig1D')(channel, variable, trig)}
   
    keylist_data = utils.get_key_list(file_data, inherits=['TH1'])
    keylist_mc = utils.get_key_list(file_mc, inherits=['TH1'])
    
    for k in keylist_data:
        if k not in keylist_mc:
            m = 'Histogram {} was present in data but not in MC.\n'.format(k)
            m += 'This is possible, but highly unlikely (based on statistics).\n'
            m += 'Check everything is correct, and edit this check if so.\n'
            m += 'Data file: {}\n'.format(name_data)
            m += 'MC file: {}'.format(name_mc)
            raise ValueError(m)
   
    keys_to_remove = []
    for k in keylist_mc:
        if k not in keylist_data:
            histo = utils.get_root_object(k, file_mc)
            stats_cut = 10
            if histo.GetEntries() < stats_cut:
                keys_to_remove.append(k)
            else:
                m = '1D Histogram {} was present in MC but not in data.\n'.format(k)
                m += 'The current statistics cut is {},'.format(stats_cut)
                m += ' but this one had {} events.'.format(histo.GetNbinsX())
                m += 'Check everything is correct, and edit this check if so.'
                raise ValueError(m)
   
    for k in keys_to_remove:
        keylist_mc.remove(k)
    assert(set(keylist_data)==set(keylist_mc))
      
    hdata1D, hmc1D = ({} for _ in range(2))
    hdata1D['ref'] = utils.get_root_object(hnames1D['ref'], file_data)
    hmc1D['ref'] = utils.get_root_object(hnames1D['ref'], file_mc)   
    hdata1D['trig'], hmc1D['trig'] = ({} for _ in range(2))

    for key in keylist_mc:
        restr1 = utils.rewrite_cut_string(hnames1D['trig'], '')
        if key.startswith(restr1):
            hdata1D['trig'][key] = utils.get_root_object(key, file_data)
            hmc1D['trig'][key] = utils.get_root_object(key, file_mc)

            # "solve" TGraphasymmErrors bin skipping when denominator=0
            # see TGraphAsymmErrors::Divide() (the default behaviour is very error prone!)
            # https://root.cern.ch/doc/master/classTGraphAsymmErrors.html#a37a202762b286cf4c7f5d34046be8c0b
            for h1D in (hdata1D,hmc1D):
                for ibin in range(1,h1D['trig'][key].GetNbinsX()+1):
                    if ( h1D['trig'][key].GetBinContent(ibin)==0. and
                         h1D['ref'].GetBinContent(ibin)==0. ):
                        h1D['ref'].SetBinContent(ibin, 1)

    # some triggers or their intersection naturally never fire for some channels
    # example: 'IsoMu24' for the etau channel
    if len(hmc1D['trig']) == 0:
        m = ('WARNING [draw_eff_and_sf_1d]: Trigger {} never '.format(trig) +
             'fired for channel={}, variable={} in MC.'
              .format(variable, channel))
        print(m)
        return
    
    data1D, mc1D, sf1D = ({} for _ in range(3)) 
    data1D_new, mc1D_new, sf1D_new = ({} for _ in range(3))
    fit_sigmoid_data, fit_sigmoid_mc, fit_sigmoid_ratio = ({} for _ in range(3))
    for atype in ('eff', 'norm'):
        data1D[atype], data1D_new[atype] = ({} for _ in range(2))
        mc1D[atype], mc1D_new[atype] = ({} for _ in range(2))
        sf1D[atype], sf1D_new[atype] = ({} for _ in range(2))

    try:
        for kh, vh in hdata1D['trig'].items():
            data1D['eff'][kh] = ROOT.TGraphAsymmErrors(vh, hdata1D['ref'])
            data1D['norm'][kh] = vh.Clone('norm_' + vh.GetName() + '_' + kh)
            try:
                data1D['norm'][kh].Scale(1/data1D['norm'][kh].Integral())
            except ZeroDivisionError:
                data1D['norm'][kh].Scale(1)
            data1D['norm'][kh] = ROOT.TGraphAsymmErrors(data1D['norm'][kh])
        for kh, vh in hmc1D['trig'].items():
            if hmc1D['ref'].GetBinContent(0) < 0:
                print("[WARNING] Had to set the overflow to zero: {} ({})".format(kh, hmc1D['ref'].GetBinContent(0)))
                hmc1D['ref'].SetBinContent(0, 0.)
            # if kh == 'Trig1D_mutau_met_et_METNoMu120_CUTS_mhtnomu_et_>_100':
            #     breakpoint()
            mc1D['eff'][kh] = ROOT.TGraphAsymmErrors(vh, hmc1D['ref'])
            mc1D['norm'][kh] = vh.Clone('norm_' + vh.GetName() + '_' + kh)
            try:
                mc1D['norm'][kh].Scale(1/mc1D['norm'][kh].Integral())
            except ZeroDivisionError:
                mc1D['norm'][kh].Scale(1)
            mc1D['norm'][kh] = ROOT.TGraphAsymmErrors(mc1D['norm'][kh])            
            
    except SystemError:
        m = 'There is likely a mismatch in the number of bins.'
        raise RuntimeError(m)

    tmpkey = list(hmc1D['trig'].keys())[0] #they all have the same binning
    nb1D = mc1D['eff'][tmpkey].GetN() 

    up_xedge  = lambda obj,i: obj.GetXaxis().GetBinUpEdge(i)
    low_xedge = lambda obj,i: obj.GetXaxis().GetBinLowEdge(i)
    up_yedge  = lambda obj,i: obj.GetYaxis().GetBinUpEdge(i)
    low_yedge = lambda obj,i: obj.GetYaxis().GetBinLowEdge(i)
    ctr_xedge = lambda obj,i: obj.GetXaxis().GetBinCenter(i)
    ctr_yedge = lambda obj,i: obj.GetYaxis().GetBinCenter(i)
        
    darr = lambda x : np.array(x).astype(dtype=np.double)
    assert len(data1D['eff']) == len(mc1D['eff'])

    # 1-dimensional
    for atype in ('eff', 'norm'):
        for (kmc,vmc),(kdata,vdata) in zip(mc1D[atype].items(),data1D[atype].items()):
            assert(kmc == kdata)

            x, y, exu, exd, eyu, eyd = (utils.dot_dict({'dt': [[] for _ in range(nb1D)],    # data 
                                                        'mc': [[] for _ in range(nb1D)],    # MC
                                                        'sf': [[] for _ in range(nb1D)] })  # scale factor
                                        for _ in range(6))

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
                    print('X MC: xp[{}] = {} +{}/-{} '
                          .format(i,x.mc[i],exu.mc[i],exd.mc[i]), flush=True)
                    print('Y MC: yp[{}] = {} +{}/-{} '
                          .format(i,y.mc[i],eyu.mc[i],eyd.mc[i]), flush=True)
                    print('X Data: xp[{}] = {} +{}/-{} '
                          .format(i,x.dt[i],exu.dt[i],exd.dt[i]), flush=True)
                    print('Y Data: yp[{}] = {} +{}/-{}'
                          .format(i,y.dt[i],eyu.dt[i],eyd.dt[i]), flush=True)
   
                x.sf[i] = x.mc[i]
                exu.sf[i] = exu.mc[i]
                exd.sf[i] = exd.mc[i]
                if x.dt[i] != x.mc[i]:
                    breakpoint()
                assert x.dt[i] == x.mc[i]
                assert exu.dt[i] == exu.mc[i]
                assert exd.dt[i] == exd.mc[i]
   
                try:
                    y.sf[i] = y.dt[i] / y.mc[i]
                except ZeroDivisionError:
                    print('There was a division by zero!', flush=True)
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

            # smoothing fits
            sf1D[atype][kdata] = ROOT.TGraphAsymmErrors(nb1D, darr(x.sf), darr(y.sf),
                                                        darr(exd.sf), darr(exu.sf), darr(eyd.sf), darr(eyu.sf))

            if atype == 'eff' and variable in cfg.fit_vars:
                fit_data_name = _fit_pp("fit_sigmoid_data_"+kdata)
                fit_sigmoid_data[kdata] = ROOT.TF1(fit_data_name, "[2]/(1+exp(-[0]*(x-[1])))", 50., 350.)
                fit_sigmoid_data[kdata].SetLineColor(ROOT.kBlack)
                fit_sigmoid_data[kdata].SetParameters(1.,175.,1.)
                data1D[atype][kdata].Fit(fit_data_name)

                fit_mc_name = _fit_pp("fit_sigmoid_mc_"+kdata)
                fit_sigmoid_mc[kdata] = ROOT.TF1(fit_mc_name, "[2]/(1+exp(-[0]*(x-[1])))", 50., 350.)
                fit_sigmoid_mc[kdata].SetLineColor(ROOT.kRed)
                fit_sigmoid_mc[kdata].SetParameters(1.,175.,1.)
                mc1D[atype][kdata].Fit(fit_mc_name)

                def _ratio_func(x,par):
                    xx = x[0]
                    if fit_sigmoid_mc[kdata].Eval(xx) == 0.:
                        return 0.;
                    return fit_sigmoid_data[kdata].Eval(xx) / fit_sigmoid_mc[kdata].Eval(xx);

                # fit_sf_div = fit_data_name+" / "+fit_mc_name
                # fit_sigmoid_ratio[kdata] = ROOT.TF1("fit_ratio_"+kdata, fit_sf_div, 50., 350.)
                fit_sigmoid_ratio[kdata] = ROOT.TF1("fit_ratio_"+kdata, _ratio_func, 50., 350.)
                fit_sigmoid_ratio[kdata].SetLineColor(ROOT.kBlue)
        
    if debug:
        print('[=debug=] 1D Plotting...', flush=True)

    n1dt, n1mc, n1sf, n1sigmfunc = 'Data1D', 'MC1D', 'SF1D', 'SigmoidFunc'
    for akey in sf1D['eff']:
        canvas_name = os.path.basename(save_names_1D[0]).split('.')[0]
        canvas_name = utils.rewrite_cut_string(canvas_name, akey, regex=True)
        canvas = {}
        canvas['eff'] = ROOT.TCanvas(canvas_name, 'canvas', 600, 600)
        canvas['norm'] = ROOT.TCanvas(canvas_name + '_norm', 'canvas_norm', 600, 600)

        for atype in ('eff', 'norm'):
            canvas[atype].cd()
            pad1 = ROOT.TPad('pad1', 'pad1', 0, 0.35, 1, 1)
            pad1.SetBottomMargin(0.005)
            pad1.SetLeftMargin(0.2)
            pad1.Draw()
            pad1.cd()

            max1, min1 = utils.get_obj_max_min(data1D[atype][akey], is_histo=False)
            max2, min2 = utils.get_obj_max_min(mc1D[atype][akey], is_histo=False)
            amax = max([max1, max2])
            amin = min([min1, min2])
            if amax == amin:
                amax = 1.
                amin = 0.

            nbins = data1D[atype][akey].GetN()
            #axor_info = nbins+1, -1, nbins
            # axor_info = nbins+1, data1D[atype][akey].GetPointX(0), data1D[atype][akey].GetPointX(nbins-1)
            # axor_ndiv = 605, 705
            # axor = ROOT.TH2D('axor'+akey+atype,'axor'+akey+atype,
            #                  axor_info[0], axor_info[1], axor_info[2],
            #                  100, amin-0.1*(amax-amin), amax+0.4*(amax-amin) )
            # axor.GetYaxis().SetTitle('Efficiency' if atype=='eff' else 'Normalized counts')
            # axor.GetYaxis().SetNdivisions(axor_ndiv[1])
            # axor.GetXaxis().SetLabelOffset(1)
            # axor.GetXaxis().SetNdivisions(axor_ndiv[0])
            # axor.GetYaxis().SetTitleSize(0.08)
            # axor.GetYaxis().SetTitleOffset(.85)
            # axor.GetXaxis().SetLabelSize(0.07)
            # axor.GetYaxis().SetLabelSize(0.07)
            # axor.GetXaxis().SetTickLength(0)
            # axor.Draw()

            # display histograms with equal width binning
            # data1D_new[atype][akey] = utils.apply_equal_bin_width(data1D[atype][akey])
            # mc1D_new[atype][akey] = utils.apply_equal_bin_width(mc1D[atype][akey])
            data1D[atype][akey].GetYaxis().SetTitleSize(0.2)
            data1D[atype][akey].GetYaxis().SetTitleOffset(1.0)
            data1D[atype][akey].GetYaxis().SetTitle('Efficiency' if atype=='eff' else 'Normalized counts')

            data1D[atype][akey].SetLineColor(1)
            data1D[atype][akey].SetLineWidth(2)
            data1D[atype][akey].SetMarkerColor(1)
            data1D[atype][akey].SetMarkerSize(1.3)
            data1D[atype][akey].SetMarkerStyle(20)
            data1D[atype][akey].GetYaxis().SetLabelSize(0.05)
            data1D[atype][akey].GetXaxis().SetLabelSize(0.05)
            data1D[atype][akey].GetXaxis().SetTitleSize(0.05)
            data1D[atype][akey].GetYaxis().SetTitleSize(0.05)
            data1D[atype][akey].GetXaxis().SetTitleOffset(1.)
            if atype == 'eff':
                data1D[atype][akey].GetYaxis().SetRangeUser(-0.05,1.2)
            data1D[atype][akey].Draw("AP")

            mc1D[atype][akey].SetLineColor(ROOT.kRed)
            mc1D[atype][akey].SetLineWidth(2)
            mc1D[atype][akey].SetMarkerColor(ROOT.kRed)
            mc1D[atype][akey].SetMarkerSize(1.3)
            mc1D[atype][akey].SetMarkerStyle(22)
            mc1D[atype][akey].Draw("same P")

            # pad1.RedrawAxis()
            # l = ROOT.TLine()
            # l.SetLineWidth(2)
            # padmin = amin-0.1*(amax-amin)
            # padmax = amax+0.1*(amax-amin)
            # fraction = (padmax-padmin)/45

            # for i in range(nbins):
            #     x = axor.GetXaxis().GetBinLowEdge(i) + 1.5;
            #     l.DrawLine(x,padmin-fraction,x,padmin+fraction)
            # l.DrawLine(x+1,padmin-fraction,x+1,padmin+fraction)

            leg = ROOT.TLegend(0.70, 0.77, 0.90, 0.87)
            leg.SetFillColor(0)
            leg.SetShadowColor(0)
            leg.SetBorderSize(0)
            leg.SetTextSize(0.05)
            leg.SetFillStyle(0)
            leg.SetTextFont(42)

            leg.AddEntry(data1D[atype][akey], 'Data', 'p')
            leg.AddEntry(mc1D[atype][akey],
                         proc.replace('MC_', '').replace('TT', 'TTbar').replace('_', '+'), 'p')
            leg.Draw('same')

            utils.redraw_border()

            lX, lY, lYstep = 0.23, 0.84, 0.05
            l = ROOT.TLatex()
            l.SetNDC()
            l.SetTextFont(72)
            l.SetTextColor(1)

            latexChannel = copy(channel)
            latexChannel.replace('mu','#mu')
            latexChannel.replace('tau','#tau_{h}')
            latexChannel.replace('Tau','#tau_{h}')

            l.DrawLatex( lX, lY, 'Channel: '+latexChannel)
            textrigs = utils.write_trigger_string(trig, intersection_str,
                                                  items_per_line=2)
            l.DrawLatex( lX, lY-lYstep, textrigs)

            canvas[atype].cd()
            pad2 = ROOT.TPad('pad2','pad2',0,0.0,1,0.35)
            pad2.SetTopMargin(0.01)
            pad2.SetBottomMargin(0.4)
            pad2.SetLeftMargin(0.2)
            pad2.Draw()
            pad2.cd()
            pad2.SetGridy()

            if max(y.sf) == min(y.sf):
                max1, min1 = 1.1, -0.1
            else:
                max1, min1 = max(y.sf)+max(eyu.sf),min(y.sf)-max(eyd.sf)

            # axor2 = ROOT.TH2D('axor2'+akey+atype, 'axor2'+akey+atype,
            #                   axor_info[0], axor_info[1], axor_info[2],
            #                   100, min1-0.1*(max1-min1), max1+0.1*(max1-min1))
            # axor2.GetXaxis().SetNdivisions(axor_ndiv[0])
            # axor2.GetYaxis().SetNdivisions(axor_ndiv[1])

            # # Change bin labels
            # rounding = lambda x: int(x) if 'pt' in variable else round(x, 2)
            # for i,elem in enumerate(sf1D[atype][akey].GetX()):
            #     axor2.GetXaxis().SetBinLabel(i+1, str(rounding(elem-sf1D[atype][akey].GetErrorXlow(i))))
            # axor2.GetXaxis().SetBinLabel(i+2, str(rounding(elem+sf1D[atype][akey].GetErrorXlow(i))))

            # axor2.GetYaxis().SetLabelSize(0.12)
            # axor2.GetXaxis().SetLabelSize(0.18)
            # axor2.GetXaxis().SetLabelOffset(0.015)
            # axor2.SetTitleSize(0.13,'X')
            # axor2.SetTitleSize(0.12,'Y')
            # axor2.GetXaxis().SetTitleOffset(1.)
            # axor2.GetYaxis().SetTitleOffset(0.5)
            # axor2.GetYaxis().SetTitle('Data/MC')
            
            # axor2.GetXaxis().SetTitle( utils.get_display_variable_name(channel, variable) )
            # axor2.GetXaxis().SetTickLength(0)
            # axor2.Draw()

            # sf1D_new[atype] = utils.apply_equal_bin_width(sf1D[atype][akey])
            sf1D[atype][akey].SetLineColor(ROOT.kBlue)
            sf1D[atype][akey].SetLineWidth(2)
            sf1D[atype][akey].SetMarkerColor(ROOT.kBlue)
            sf1D[atype][akey].SetMarkerSize(1.3)
            sf1D[atype][akey].SetMarkerStyle(22)
            sf1D[atype][akey].GetYaxis().SetLabelSize(0.1)
            sf1D[atype][akey].GetXaxis().SetLabelSize(0.1)
            sf1D[atype][akey].GetXaxis().SetTitleSize(0.1)
            sf1D[atype][akey].GetYaxis().SetTitleSize(0.1)
            sf1D[atype][akey].GetXaxis().SetTitleOffset(1.)
            sf1D[atype][akey].GetYaxis().SetTitleOffset(0.45)
            sf1D[atype][akey].GetYaxis().SetTitle('Data/MC')
            sf1D[atype][akey].GetXaxis().SetTitle( utils.get_display_variable_name(channel, variable) )
            sf1D[atype][akey].Draw("AP")
            if atype == 'eff' and variable in cfg.fit_vars:
                fit_sigmoid_ratio[akey].Draw("same")

            pad2.cd()
            l = ROOT.TLine()
            l.SetLineWidth(2)
            # padmin = min1-0.1*(max1-min1)
            # padmax = max1+0.1*(max1-min1)
            # fraction = (padmax-padmin)/30

            # for i in range(nbins):
            #     x = axor2.GetXaxis().GetBinLowEdge(i) + 1.5;
            #     l.DrawLine(x,padmin-fraction,x,padmin+fraction)
            # l.DrawLine(x+1,padmin-fraction,x+1,padmin+fraction)

            utils.redraw_border()

            for aname in save_names_1D[:-1]:
                _name = utils.rewrite_cut_string(aname, akey, regex=True)
                _name = _name.replace(args.canvas_prefix, atype + '_' + args.canvas_prefix)
                canvas[atype].SaveAs( _name )

        _name_base = utils.rewrite_cut_string(save_names_1D[-1], akey, regex=True)
        for atype in ('eff',):
            _name = _name_base.replace(args.canvas_prefix, atype + '_')
            afile = ROOT.TFile.Open(_name, 'RECREATE')
            afile.cd()

            # save the original histograms, not the ones for visualization
            data1D[atype][akey].SetName(n1dt)
            data1D[atype][akey].Write(n1dt)
            mc1D[atype][akey].SetName(n1mc)
            mc1D[atype][akey].Write(n1mc)
            sf1D[atype][akey].SetName(n1sf)
            sf1D[atype][akey].Write(n1sf)

            if variable in cfg.fit_vars:
                fit_sigmoid_data[akey].Write(n1sigmfunc+"Data")
                fit_sigmoid_mc[akey].Write(n1sigmfunc+"MC")
                fit_sigmoid_ratio[akey].Write(n1sigmfunc+"SF")

    if debug:
        print('[=debug=] 2D Plotting...', flush=True)


def draw_eff_and_sf_2d(proc, channel, joinvars, trig, save_names_2D,
                       tprefix, indir, subtag, mc_name, data_name,
                       intersection_str, debug):
    _name = lambda a,b,c,d : a + b + c + d + '.root'

    name_data = os.path.join(indir, _name( tprefix, data_name, '_Sum', subtag ) )

    file_data = ROOT.TFile.Open(name_data, 'READ')

    name_mc = os.path.join(indir, _name( tprefix, mc_name, '_Sum', subtag ))
    file_mc   = ROOT.TFile.Open(name_mc, 'READ')

    if debug:
        print('[=debug=] Open files:')
        print('[=debug=]  - Data: {}'.format(name_data))
        print('[=debug=]  - MC: {}'.format(name_mc))
        print('[=debug=]  - Args: proc={proc}, channel={channel}'
              .format(proc=proc, channel=channel))
        print('[=debug=]  - Args: joinvars={}, trig_inters={}'
              .format(joinvars, trig))
   
    hnames2D = { 'ref':  utils.get_hnames('Ref2D')(channel, joinvars, trig),
                 'trig': utils.get_hnames('Trig2D')(channel, joinvars, trig) }
   
    keylist_data = utils.get_key_list(file_data, inherits=['TH1'])
    keylist_mc = utils.get_key_list(file_mc, inherits=['TH1'])

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
            histo = utils.get_root_object(k, file_mc)
            stats_cut = 10
            if histo.GetEntries() < stats_cut:
                keys_to_remove.append(k)
            else:
                m = '2D Histogram {} was present in MC but not in data.\n'.format(k)
                m += 'The current statistics cut is {},'.format(stats_cut)
                m += ' but this one had {} events.'.format(histo.GetNbinsX())
                m += 'Check everything is correct, and edit this check if so.'
                raise ValueError(m)
   
    for k in keys_to_remove:
        keylist_mc.remove(k)
    assert(set(keylist_data)==set(keylist_mc))

    hdata2D, hmc2D = ({} for _ in range(2))
    hdata2D['ref'] = utils.get_root_object(hnames2D['ref'], file_data)
    hmc2D['ref'] = utils.get_root_object(hnames2D['ref'], file_mc)   
    hdata2D['trig'], hmc2D['trig'] = ({} for _ in range(2))

    for key in keylist_mc:
        restr2 = utils.rewrite_cut_string(hnames2D['trig'], '')
        if key.startswith(restr2):
            hdata2D['trig'][key] = utils.get_root_object(key, file_data)
            hmc2D['trig'][key] = utils.get_root_object(key, file_mc)

            for h2D in [hdata2D,hmc2D]:
                for ix in range(1,h2D['trig'][key].GetNbinsX()+1):
                    for iy in range(1,h2D['trig'][key].GetNbinsY()+1):
                        if ( h2D['trig'][key].GetBinContent(ix,iy)==0. and
                             h2D['ref'].GetBinContent(ix,iy)==0. ):
                            h2D['ref'].SetBinContent(ix,iy,1)

    # some triggers or their intersection naturally never fire for some channels
    # example: 'IsoMu24' for the etau channel
    if  len(hmc2D['trig']) == 0:
        print('WARNING: Trigger {} never fired for channel {} in MC.'
              .format(trig, channel))
        return
    
    effdata2D, effmc2D = ({} for _ in range(2))
    try:
        for kh, vh in hdata2D['trig'].items():
            effdata2D[kh] = vh.Clone(vh.GetName() + '_' + kh)
            effdata2D[kh].Divide(effdata2D[kh], hdata2D['ref'], 1, 1, "B")

        for kh, vh in hmc2D['trig'].items():
            effmc2D[kh] = vh.Clone(vh.GetName() + '_' + kh)
            effmc2D[kh].Divide(effmc2D[kh], hmc2D['ref'], 1, 1, "B")
            
    except SystemError:
        m = 'There is likely a mismatch in the number of bins.'
        raise RuntimeError(m)

    tmpkey = list(hmc2D['trig'].keys())[0]
    nb2Dx = effmc2D[tmpkey].GetNbinsX()
    nb2Dy = effmc2D[tmpkey].GetNbinsY()
        
    sf2D = {}
    assert len(effdata2D) == len(effmc2D)

    for kh, vh in effdata2D.items():
        sf2D[kh] = vh.Clone(vh.GetName() + '_' + kh)
        sf2D[kh].Sumw2()
        sf2D[kh].Divide(effmc2D[kh])
        
    if debug:
        print('[=debug=] 2D Plotting...', flush=True)

    n2 = ['Data2D', 'MC2D', 'SF2D']
    cnames = [x + '_c' for x in n2]
    ROOT.gStyle.SetPaintTextFormat("4.3f");

    base_name = ( save_names_2D[joinvars][n2[0]][0]
                  .split('.')[0]
                  .replace('EffData_', '')
                  .replace('Canvas2D_', '') )

    hzip = zip(effdata2D.items(), effmc2D.items(), sf2D.items())
    for items in hzip:
        # save 2D histograms
        _name = utils.rewrite_cut_string(base_name, items[0][0], regex=True)
        _name += '.root'

        eff_file_2D = ROOT.TFile.Open(_name, 'RECREATE')
        eff_file_2D.cd()

        for itype, obj in enumerate(items):
            obj[1].SetName(n2[itype])
            obj[1].Write(n2[itype])
  
            # upper and lower 2D uncertainties
            eff_eu = obj[1].Clone(n2[itype] + '_eu')
            eff_ed = obj[1].Clone(n2[itype] + '_ed')
            for i in range(1,obj[1].GetNbinsX()+1):
                for j in range(1,obj[1].GetNbinsY()+1):
                    abin = obj[1].GetBin(i, j)
                    eu2d = obj[1].GetBinErrorLow(abin)
                    ed2d = obj[1].GetBinErrorUp(abin)
                    if obj[1].GetBinContent(abin)==0.:
                        eff_eu.SetBinContent(abin, 0.)
                        eff_ed.SetBinContent(abin, 0.)
                    else:
                        eff_eu.SetBinContent(abin, eu2d)
                        eff_ed.SetBinContent(abin, ed2d)
      
                    if eff_eu.GetBinContent(abin)==0.:
                        eff_eu.SetBinContent(abin, 1.e-10)
                    if eff_ed.GetBinContent(abin)==0.:
                        eff_ed.SetBinContent(abin, 1.e-10)

            canvas = ROOT.TCanvas(cnames[itype]+obj[0], cnames[itype]+obj[0], 600, 600)
            canvas.SetLeftMargin(0.10)
            canvas.SetRightMargin(0.15);
            canvas.cd()
      
            vnames_2D = utils.split_vnames(joinvars)
            vname_x = utils.get_display_variable_name(channel, vnames_2D[0])
            vname_y = utils.get_display_variable_name(channel, vnames_2D[1])
            rx = 2 if 'pt' not in vname_x else 0
            ry = 2 if 'pt' not in vname_y else 0
            veff_new = utils.apply_equal_bin_width(obj[1], roundx=rx, roundy=ry)
            veff_new.GetXaxis().SetTitleOffset(1.0)
            veff_new.GetYaxis().SetTitleOffset(1.3)
            veff_new.GetXaxis().SetTitleSize(0.03)
            veff_new.GetYaxis().SetTitleSize(0.03)
            veff_new.GetXaxis().SetLabelSize(0.03)
            veff_new.GetYaxis().SetLabelSize(0.03)
            veff_new.GetXaxis().SetTitle(vname_x)
            veff_new.GetYaxis().SetTitle(vname_y)
            # useful when adding errors
            # check scripts/draw2DTriggerSF.py
            veff_new.SetBarOffset(0.22)
            veff_new.SetMarkerSize(.75)
            veff_new.SetMarkerColor(ROOT.kOrange+10)
            veff_new.SetMarkerSize(.8)
            ROOT.gStyle.SetPaintTextFormat("4.3f");
            veff_new.Draw('colz text')
      
            # numerator and denominator
            if itype == 0: # data
                htot = hdata2D['ref'].Clone('tot')
                hpass = hdata2D['trig'][obj[0]].Clone('pass')
            elif itype == 1: # MC
                htot = hmc2D['ref'].Clone('tot')
                hpass = hmc2D['trig'][obj[0]].Clone('pass')

            if itype < 2:
                htot.SetName('tot_' + str(itype))
                hpass.SetName('pass_' + str(itype))
                htot = utils.apply_equal_bin_width(htot)
                hpass = utils.apply_equal_bin_width(hpass)
                htot.SetMarkerSize(.65)
                hpass.SetMarkerSize(.65)
                htot.SetMarkerColor(ROOT.kOrange+6)
                hpass.SetMarkerColor(ROOT.kOrange+6)
                hpass.SetBarOffset(-0.10)
                htot.SetBarOffset(-0.22)
                ROOT.gStyle.SetPaintTextFormat("4.3f");
                hpass.Draw("same text")
                htot.Draw("same text")
      
            # up and down errors for the 2D histogram
            eff_eu = utils.apply_equal_bin_width(eff_eu)
            eff_ed = utils.apply_equal_bin_width(eff_ed)
            eff_eu.SetMarkerSize(.6)
            eff_ed.SetMarkerSize(.6)
            eff_eu.SetBarOffset(0.323);
            eff_ed.SetBarOffset(0.10);
            eff_eu.SetMarkerColor(ROOT.kBlack)
            eff_ed.SetMarkerColor(ROOT.kBlack)
            ROOT.gStyle.SetPaintTextFormat("+ 4.3f xxx");
            eff_eu.Draw("same text")
            ROOT.gStyle.SetPaintTextFormat("- 4.3f");
            eff_ed.Draw("same text")
      
            lX, lY, lYstep = 0.8, 0.92, 0.045
            l = ROOT.TLatex()
            l.SetNDC()
            l.SetTextFont(72)
            l.SetTextColor(2)
            textrig = utils.write_trigger_string(trig, intersection_str,
                                                 items_per_line=2)
            l.DrawLatex(lX, lY, n2[itype].replace('2D',''))
      
            paint2d(channel, textrig)
            utils.redraw_border()

            for full in save_names_2D[joinvars][n2[itype]]:
                full = utils.rewrite_cut_string(full, obj[0], regex=True)
                full = full.replace('Canvas2D_', '')
                canvas.SaveAs(full)

def _fit_pp(s):
    return s.replace('>', 'G').replace('<', 'L').replace('=', 'EQ')

def _get_canvas_name(prefix, proc, chn, var, trig, data_name, subtag):
    """
    A 'XXX' placeholder is added for later replacement by all cuts considered for the same channel, variable and trigger combination.
    Without the placeholder one would have to additionally calculate the number of cuts beforehand, which adds complexity with no major benefit.
    """
    add = proc + '_' + chn + '_' + var + '_TRG_' + trig
    n = prefix + data_name + '_' + add
    n += main.placeholder_cuts + subtag
    return n

def run_eff_sf_1d_outputs(outdir, data_name, mc_name,
                          tcomb, channels, variables, subtag):
    outputs = [[] for _ in range(len(main.extensions))]
    processes = [mc_name] #CHANGE !!!! IF STUDYING MCs SEPARATELY
    
    for proc in processes:
        for ch in channels:
            for var in variables:
                canvas_name = _get_canvas_name(args.canvas_prefix,
                                               proc, ch, var,
                                               tcomb,
                                               data_name, subtag)
                thisbase = os.path.join(outdir, ch, var, '')
                utils.create_single_dir(thisbase)

                for ext,out in zip(main.extensions, outputs):
                    _out = os.path.join(thisbase, canvas_name + '.' + ext)
                    out.append(_out)                    

    #join all outputs in the same list
    return sum(outputs, []), main.extensions, processes

def run_eff_sf_2d_outputs(outdir, proc, data_name, cfg, tcomb,
                          channel, subtag, intersection_str, debug):
    """
    This output function is not ready to be used in the luigi framework.
    It returns the outputs corresponding to a single trigger combination.
    This was done for the sake of simplification.
    """
    outputs = {}
    
    # only running when one of the triggers in the intersection
    # matches one of the triggers specified by the user with `cfg.pairs2D`
    splits = tcomb.split(intersection_str)
    run = any({x in cfg.pairs2D.keys() for x in splits})

    if run:
        for onetrig in splits:
            if onetrig in cfg.pairs2D:
                for j in cfg.pairs2D[onetrig]:
                    vname = utils.add_vnames(j[0],j[1])
                    pref2d = args.canvas_prefix.replace('1', '2')
                    cname = _get_canvas_name(pref2d,
                                             proc, channel, vname,
                                             tcomb,
                                             data_name, subtag)

                    outputs[vname] = {}
                    
                    cnames = { 'Data2D': cname.replace(pref2d,'EffData_'),
                               'MC2D': cname.replace(pref2d,'EffMC_'),
                               'SF2D': cname.replace(pref2d,'SF_') }

                    thisbase = os.path.join(outdir, channel, vname, '')

                    utils.create_single_dir(thisbase)

                    for kn,vn in cnames.items():
                        outputs[vname][kn] = []
                        for ext in main.extensions:
                            if ext not in ('C', 'root'): #remove extra outputs to reduce noise
                                outputs[vname][kn].append(os.path.join(thisbase, vn + '.' + ext))
    
    if debug:
        for k,out in outputs.items():
            print(k)
            for o in out:
                print(o)
                print()

    return outputs


def run_eff_sf_1d(indir, outdir, data_name, mc_name, configuration,
                  tcomb, channels, variables, subtag,
                  tprefix, intersection_str, debug):
    
    outs1D, extensions, processes = run_eff_sf_1d_outputs(outdir, data_name, mc_name, tcomb,
                                                          channels, variables, subtag)

    config_module = importlib.import_module(configuration)
    
    triggercomb = {}
    for chn in channels:
        triggercomb[chn] = utils.generate_trigger_combinations(chn, config_module.triggers,
                                                               config_module.exclusive)

    dv = len(args.variables)
    dc = len(args.channels) * dv
    dp = len(processes) * dc

    for ip,proc in enumerate(processes):
        for ic,chn in enumerate(channels):
            if not utils.is_trigger_comb_in_channel(chn, tcomb, config_module.triggers,
                                                    config_module.exclusive):
                continue
            for iv,var in enumerate(variables):
                index = ip*dc + ic*dv + iv
                names1D = [ outs1D[index + dp*x] for x in range(len(extensions)) ]

                if args.debug:
                    for name in names1D:
                        print('[=debug=] {}'.format(name))
                        m = ( 'process={}, channel={}, variable={}'
                              .format(proc, chn, var) )
                        m += ', trigger_combination={}\n'.format(tcomb)
                        print(m)

                draw_eff_and_sf_1d(proc, chn, var, tcomb, names1D, config_module,
                                   tprefix, indir, subtag, mc_name, data_name,
                                   intersection_str, debug)

    splits = tcomb.split(intersection_str)
    for x in splits:
        if x not in main.trig_map:
            mess = 'Trigger {} was not defined in the configuration.'.format(x)
            raise ValueError(mess)
        
    run = any({x in config_module.pairs2D for x in splits})
    if run:
        for ip,proc in enumerate(processes):
            for ic,chn in enumerate(channels):
                if not utils.is_trigger_comb_in_channel(chn, tcomb, config_module.triggers,
                                                        config_module.exclusive):
                    continue
                
                names2D = run_eff_sf_2d_outputs(outdir, proc, data_name, config_module, tcomb, 
                                                chn, subtag, intersection_str, debug)

                for onetrig in splits:
                    if onetrig in config_module.pairs2D:

                        for j in config_module.pairs2D[onetrig]:
                            vname = utils.add_vnames(j[0],j[1])

                            draw_eff_and_sf_2d(proc, chn, vname, tcomb, names2D,
                                               tprefix, indir, subtag, mc_name, data_name,
                                               intersection_str, debug)

parser = argparse.ArgumentParser(description='Draw trigger scale factors')
parser.add_argument('--indir', help='Inputs directory', required=True)
parser.add_argument('--outdir', help='Output directory', required=True, )
parser.add_argument('--tprefix', help='prefix to the names of the produceyd outputs (targets in luigi lingo)', required=True)
parser.add_argument('--canvas_prefix', help='canvas prefix', required=True)
parser.add_argument('--subtag', dest='subtag', required=True, help='subtag')
parser.add_argument('--mc_name', dest='mc_name', required=True, type=str,
                    help='Id for all MC samples')
parser.add_argument('--data_name', dest='data_name', required=True, type=str,
                    help='Id for all data samples',)
parser.add_argument('--triggercomb', dest='triggercomb', required=True,
                    help='Trigger intersection combination.')
parser.add_argument('--channels',   dest='channels',         required=True, nargs='+', type=str,
                    help='Select the channels over which the workflow will be run.' )
parser.add_argument('--variables',        dest='variables',        required=True, nargs='+', type=str,
                    help='Select the variables over which the workflow will be run.' )
parser.add_argument('--intersection_str', dest='intersection_str', required=False, default=main.inters_str,
                    help='String useyd to represent set intersection between triggers.')
parser.add_argument('--configuration', dest='configuration', required=True,
                    help='Name of the configuration module to use.')
parser.add_argument('--debug', action='store_true', help='debug verbosity')
args = utils.parse_args(parser)

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
run_eff_sf_1d(args.indir, args.outdir,
              args.data_name, args.mc_name,
              args.configuration,
              args.triggercomb,
              args.channels, args.variables,
              args.subtag,
              args.tprefix,
              args.intersection_str,
              args.debug)
