# coding: utf-8

_all_ = [ 'draw_2D_sf', 'draw_2D_sf_outputs' ]

import os
import sys
parent_dir = os.path.abspath(__file__ + 3 * '/..')
sys.path.insert(0, parent_dir)

import argparse
import numpy as np
from copy import copy

import ROOT
ROOT.gROOT.SetBatch(True)
from ROOT import (
    TCanvas
    TPad,
    TStyle,
    TFile,
    TEfficiency,
    TGraphAsymmErrors,
    TH1D,
    TH2D,
    TLatex,
    TLine,
    TLegend,
    TString,
    )

import inclusion
from inclusion.utils import utils
from inclusion import config

def set_histo_props(histo, variables):
    histo.GetYaxis().SetNdivisions(6)
    histo.GetXaxis().SetNdivisions(6)
    histo.GetYaxis().SetLabelSize(0.04)
    histo.GetXaxis().SetLabelSize(0.04)
    histo.GetXaxis().SetTitleSize(0.04)
    histo.GetYaxis().SetTitleSize(0.04)
    histo.GetXaxis().SetTitleOffset(0.)
    histo.GetYaxis().SetTitleOffset(1.25)
    histo.GetXaxis().SetTitle(variables[0])
    histo.GetYaxis().SetTitle(variables[1])

def set_histo(histo, variables):
    set_histo_props(histo, variables)
    return histo

def paint_channel_and_trigger(channel, trig):
    lX, lY, lYstep = 0.06, 0.96, 0.03
    l = TLatex()
    l.SetNDC()
    l.SetTextFont(72)
    l.SetTextSize(0.03)
    l.SetTextColor(1)
  
    latexChannel = copy(channel)
    latexChannel.replace('mu','#mu')
    latexChannel.replace('tau','#tau_{h}')
    latexChannel.replace('Tau','#tau_{h}')
    l.DrawLatex( lX, lY, 'Channel: '+latexChannel)
    l.DrawLatex( lX, lY-lYstep, 'Trigger: '+trig)

def check_2D_trigger(args, proc, channel, var, trig, save_names):
    _name = lambda a,b,c,d : a + b + c + d + '.root'
    histo_options = 'colz text'
    name_data = os.path.join(args.indir, _name( args.targetsPrefix, args.data_name,
                                                args.tsuffix, args.subtag ) )
    file_data = TFile( name_data, 'READ')
    
    name_mc = os.path.join(args.indir, _name( args.targetsPrefix, args.mc_name,
                                            args.tsuffix, args.subtag ))
    file_mc   = TFile( name_mc, 'READ');

    if args.debug:
        print('[=debug=] Open files:')
        print('[=debug=]  - Data: {}'.format(name_data))
        print('[=debug=]  - MC: {}'.format(name_mc))
        print('[=debug=]  - Args: proc={proc}, channel={channel}, variable={variable}, trig={trig}'
              .format(proc=proc, channel=channel, variable=variable, trig=trig))
        
    vname = utils.add_vnames( var[0], var[1] )
    eff_names = { 'ref_vs_trig': 'effRefVsTrig_{}_{}_{}'.format(channel, trig, vname),
                 }
    
    eff2D_mc   = { k: utils.getROOTObject(v, file_mc)   for k,v in eff_names.items() }
    eff2D_data = { k: utils.getROOTObject(v, file_data) for k,v in eff_names.items() }

    eff2D = {'Data': eff2D_data, proc: eff2D_mc}
    canvas, histos = ([] for _ in range(len(eff2D)))
    histos_eu, histos_ed = ([] for _ in range(2)) #upper and lower 2D uncertainties
    
    if args.debug:
        print('[=debug=] Plotting...')

    for effname,thiseff2D in eff2D.items():
        saveid = 0 if effname=='Data' else 1
        cname = os.path.basename(save_names[0][saveid]).split('.')[0]
        canvas.append( TCanvas(cname, 'canvas_'+effname, 600, 600) )
        canvas[-1].SetLeftMargin(0.10);
        canvas[-1].SetRightMargin(0.15);
        canvas[-1].cd()
        
        thiseff2D['ref_vs_trig'].Draw('colz')
        ROOT.gPad.Update()
        histos.append( set_histo(thiseff2D['ref_vs_trig'].GetPaintedHistogram(), var) )

        histos_eu.append( histos[-1].Clone(effname+'_eu') )
        histos_ed.append( histos[-1].Clone(effname+'_ed') )
        for i in range(1,histos[-1].GetNbinsX()+1):
            for j in range(1,histos[-1].GetNbinsY()+1):
                abin = histos[-1].GetBin(i, j)
                eu2d = thiseff2D['ref_vs_trig'].GetEfficiencyErrorLow(abin)
                ed2d = thiseff2D['ref_vs_trig'].GetEfficiencyErrorUp(abin)
                if histos[-1].GetBinContent(abin)==0.:
                    histos_eu[-1].SetBinContent(abin, 0.)
                    histos_ed[-1].SetBinContent(abin, 0.)
                else:
                    histos_eu[-1].SetBinContent(abin, eu2d)
                    histos_ed[-1].SetBinContent(abin, ed2d)

                if histos_eu[-1].GetBinContent(abin)==0.:
                    histos_eu[-1].SetBinContent(abin, 1.e-10)
                if histos_ed[-1].GetBinContent(abin)==0.:
                    histos_ed[-1].SetBinContent(abin, 1.e-10)
                    
        histo_pass = thiseff2D['ref_vs_trig'].GetCopyPassedHisto()
        histo_tot = set_histo( thiseff2D['ref_vs_trig'].GetCopyTotalHisto(), var )

        histo_tot.SetMarkerSize(.75)
        histo_tot.SetMarkerColor(ROOT.kOrange+10)
        histo_tot.SetMarkerSize(.8)
        histo_pass.SetMarkerColor(ROOT.kOrange+6)
        histo_pass.SetMarkerSize(.8)
        histos_eu[-1].SetMarkerSize(.6)
        histos_ed[-1].SetMarkerSize(.6)

        ROOT.gStyle.SetPaintTextFormat("4.3f");
        histos[-1].Draw(histo_options)

        histos_eu[-1].SetBarOffset(0.323);
        histos[-1].SetBarOffset(0.20);
        histos_ed[-1].SetBarOffset(0.10);
        histo_pass.SetBarOffset(-0.05);
        histo_tot.SetBarOffset(-0.20);

        ROOT.gStyle.SetPaintTextFormat("+ 4.3f xxx");
        histos_eu[-1].Draw("text same")
        ROOT.gStyle.SetPaintTextFormat("- 4.3f");
        histos_ed[-1].Draw("text same")
        ROOT.gStyle.SetPaintTextFormat("4.3f");
        histo_pass.Draw("text same")
        histo_tot.Draw("text same")
        ROOT.gPad.Update();
  
        lX, lY, lYstep = 0.8, 0.92, 0.045
        l = TLatex()
        l.SetNDC()
        l.SetTextFont(72)
        l.SetTextColor(2)
        l.DrawLatex(lX, lY, effname)

        paint_channel_and_trigger(channel, trig)
        utils.redrawBorder()

    canvas_sf = TCanvas(os.path.basename(save_names[0][2]).split('.')[0],
                        'canvas_sf', 600, 600 )
    canvas_sf.SetLeftMargin(0.10);
    canvas_sf.SetRightMargin(0.15);
    canvas_sf.cd()

    # Data / MC
    histo_sf = histos[0].Clone('sf')
    histo_sf.Divide(histos[1])

    histo_sf.SetAxisRange(-0.5, 2.5, 'Z');
    histo_sf.SetMarkerSize(.8)

    histo_sf_eu = histo_sf.Clone('sf_eu')
    histo_sf_ed = histo_sf.Clone('sf_ed')
    for i in range(1,histo_sf.GetNbinsX()+1):
        for j in range(1,histo_sf.GetNbinsY()+1):
            abin = histo_sf.GetBin(i, j)
            eu_data = histos_eu[0].GetBinContent(abin)
            eu_mc = histos_eu[1].GetBinContent(abin)
            ed_data = histos_ed[0].GetBinContent(abin)
            ed_mc = histos_ed[1].GetBinContent(abin)
            eu = np.sqrt(eu_data*eu_data + eu_mc*eu_mc)
            ed = np.sqrt(ed_data*ed_data + ed_mc*ed_mc)
            if histo_sf.GetBinContent(abin)==0.:
                histo_sf_eu.SetBinContent(abin, 0.)
                histo_sf_ed.SetBinContent(abin, 0.)
            else:
                histo_sf_eu.SetBinContent(abin, eu)
                histo_sf_ed.SetBinContent(abin, ed)

            # trick to display zero errors
            # ('min0' does not work: it also display zeros on empty bins)
            if histo_sf_eu.GetBinContent(abin)==0.:
                histo_sf_eu.SetBinContent(abin, 1.e-10)
            if histo_sf_ed.GetBinContent(abin)==0.:
                histo_sf_ed.SetBinContent(abin, 1.e-10)

                
    ROOT.gStyle.SetPaintTextFormat("4.3f");
    histo_sf.Draw(histo_options)
    
    histo_sf.SetBarOffset(0.0);
    histo_sf_eu.SetBarOffset(0.10);
    histo_sf_ed.SetBarOffset(-0.10);
    histo_sf_eu.SetMarkerSize(0.6);
    histo_sf_ed.SetMarkerSize(0.6);

    ROOT.gStyle.SetPaintTextFormat("+ 4.3f");
    histo_sf_eu.Draw("text same")
    ROOT.gStyle.SetPaintTextFormat("- 4.3f");
    histo_sf_ed.Draw("text same")

    lX, lY, lYstep = 0.6, 0.92, 0.045
    l = TLatex()
    l.SetNDC()
    l.SetTextFont(72)
    l.SetTextColor(2)
    l.DrawLatex(lX, lY, 'Data / {}'.format(proc))

    paint_channel_and_trigger(channel, trig)
    utils.redrawBorder()

    for aname in save_names:
        canvas[0].SaveAs( aname[0] )
        canvas[1].SaveAs( aname[1] )
        canvas_sf.SaveAs( aname[2] )
  
@utils.set_pure_input_namespace
def draw_2D_sf_outputs(args):
    outputs = [[] for _ in range(len(config.extensions))]
    processes = args.mc_processes if args.draw_independent_MCs else [args.mc_name]
  
    for proc in processes:
        for ch in args.channels:
            for trig in args.triggers:
                if trig in config.pairs2D.keys():
                    for variables in config.pairs2D[trig]:
                        add = proc + '_' + ch + '_' + trig + '_' + variables[0] + '_VS_' + variables[1]
                        canvas_data_name = 'EffData_' + args.data_name + '_' + add + args.subtag
                        canvas_mc_name = 'EffMC_' + args.data_name + '_' + add + args.subtag
                        canvas_sf_name = 'SF_' + args.data_name + '_' + add + args.subtag
                        thisbase = os.path.join(args.outdir, ch, '')
                        utils.createSingleDir( thisbase )

                    for ext,out in zip(config.extensions, outputs):
                        out.append( ( os.path.join( thisbase, canvas_data_name + '.' + ext ),
                                      os.path.join( thisbase, canvas_mc_name   + '.' + ext ),
                                      os.path.join( thisbase, canvas_sf_name   + '.' + ext )) )

    #join all outputs in the same list
    return sum(outputs, []), config.extensions
    
@utils.set_pure_input_namespace
def draw_2D_sf(args):
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)

    outputs, extensions = draw_2D_sf_outputs(args)
    processes = args.mc_processes if args.draw_independent_MCs else [args.mc_name]

    # loop through variables, triggers, channels and processes
    dv = 0
    for key in config.pairs2D:
        dv += len(config.pairs2D[key])
    dc = len(args.channels) * dv
    dp = len(processes) * dc
    for ip,proc in enumerate(processes):
        for ic,ch in enumerate(args.channels):
            iv = -1
            for trig in args.triggers:
                 if trig in config.pairs2D.keys():
                    for variables in config.pairs2D[trig]:
                        iv += 1
                        index = ip*dc + ic*dv + iv
                        names = [ outputs[index + dp*x] for x in range(len(extensions)) ]
                       
                        if args.debug:
                            for name in names:
                                print('[=debug=] {}'.format(name))
                            print("process={}, channel={}, variables={}, trigger={}".format(proc, ch, variables, trig))
                            print()
                       
                        check_2D_trigger( args, proc, ch, variables, trig, names )
          
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Draw trigger scale factors')
    parser.add_argument('-i', '--indir', help='Inputs directory', required=True)
    parser.add_argument('-x', '--targetsPrefix', help='prefix to the names of the produced outputs (targets in luigi lingo)', required=True)
    parser.add_argument('-t', '--tag', help='string to diferentiate between different workflow runs', required=True)
    parser.add_argument('-d', '--data', help='dataset to be analyzed/plotted', required=True)
    parser.add_argument('-p', '--mc_processes', help='MC processes to be analyzed: Radions, TT, ...', required=True)
    parser.add_argument('--draw_independent_MCs', action='store_true', help='debug verbosity')
    parser.add_argument('--nocut_dummy_str', dest='tprefix', required=True,
                        help='Dummy string associated to trigger histograms were no cuts are applied.')
    parser.add_argument('--debug', action='store_true', help='debug verbosity')
    args = parser.parse_args()

    draw_2D_sf(args)
