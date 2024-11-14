# coding: utf-8

_all_ = [ 'test_trigger_gains' ]

import os
import sys
parent_dir = os.path.abspath(__file__ + 2 * '/..')
sys.path.insert(0, parent_dir)

import json
import argparse
from inclusion.utils import utils
import numpy as np
from collections import defaultdict as dd
import hist
import mplhep as hep
from hist.intervals import clopper_pearson_interval as clop
import pickle
import uproot
import matplotlib.pyplot as plt


# import bokeh
# from bokeh.plotting import figure, output_file, save
# from bokeh.models import Whisker
# from bokeh.layouts import gridplot
# from bokeh.io import export_svg, export_png

tau = '\u03C4'
mu  = '\u03BC'
pm  = '\u00B1'
ditau = tau+tau

def get_outname(channel, regcuts, ptcuts, met_turnon, bigtau):
    utils.create_single_dir('data')

    name = channel + '_'
    name += '_'.join((*regcuts, 'ptcuts', *[str(x) for x in ptcuts], 'turnon', str(met_turnon)))
    if bigtau:
        name += '_BIGTAU'
    name += '.pkl'

    s = 'data/regions_{}'.format(name)
    return s

def pp(chn):
    if chn == "tautau":
        return ditau
    elif chn == "etau":
        return "e" + tau
    elif chn == "mutau":
        return mu + tau

def rec_dd():
    return dd(rec_dd)

def set_fig(fig, legend=True):
    fig.output_backend = 'svg'
    fig.toolbar.logo = None
    # if legend:
    #     fig.legend.click_policy='hide'
    #     fig.legend.location = 'top_left'
    #     fig.legend.label_text_font_size = '8pt'
    fig.min_border_bottom = 5
    fig.xaxis.visible = True
    fig.title.align = "left"
    fig.title.text_font_size = "20px"
    fig.xaxis.axis_label_text_font_style = "bold"
    fig.yaxis.axis_label_text_font_style = "bold"
    fig.xaxis.axis_label_text_font_size = "13px"
    fig.yaxis.axis_label_text_font_size = "13px"

def main(args):
    vars = [
        # 'dau1_iso', 
        'dau1_pt', 
        # 'dau2_iso',
            'dau2_pt', 'dau1_eta', 'dau2_eta',
            'bjet1_pNet', 'bjet2_pNet', 'bjet1_pt',
            'bjet2_pt', 'bjet1_eta', 'bjet2_eta']
    if args.channels[0] == "etau":
        base_extension = "E"
        channel_text = r"e$\tau_{h}$"
    elif args.channels[0] == "mutau":
        base_extension = "Mu"
        channel_text = r"$\mu \tau_{h}$"
    elif args.channels[0] == "tautau":
        base_extension = "Tau"
        channel_text = r"$\tau_{h} \tau_{h}$"
    with uproot.open("data/regions_preEE_10p20_all.root") as file:
        # from matplotlib.colors import ListedColormap
        petroff6 = ["#5790fc", "#f89c20", "#e42536", "#964a8b", "#9c9ca1", "#7a21dd"]

        channel_all = file['All_baseline_{}_genHH_mass'.format(args.channels[0])]
        channel_base = file['Base{}_baseline_{}_genHH_mass'.format(base_extension, args.channels[0])]
        # channel_base_plus_met = file['Base{}ORMET_baseline_{}_genHH_mass'.format(base_extension, args.channels[0])]
        channel_base_plus_tau = file['Base{}ORTauTauJet_baseline_{}_genHH_mass'.format(base_extension, args.channels[0])]
        # channel_base_plus_Mu50 = file['Base{}ORMu50_baseline_{}_genHH_mass'.format(base_extension, args.channels[0])]
        channel_base_plus_quadJetPNet = file['Base{}OR4JetsPNet_baseline_{}_genHH_mass'.format(base_extension, args.channels[0])]
        # channel_base_plus_met = file['NoBaseMuTTJet_baseline_{}_genHH_mass'.format(args.channels[0])]
        # channel_base_plus_quadJetDeep = file['NoBaseMu4JetsDeepJet_baseline_{}_genHH_mass'.format(args.channels[0])]
        # if args.channels[0] == "tautau":
        #     channel_base_all = file['Base{}ORTauTauJetOR4JetsPNet_baseline_{}_genHH_mass'.format(base_extension, args.channels[0])]
        # if args.channels[0] == "mutau":
        #     channel_base_all = file['Base{}ORMETORTauORMu50OR4JetsPNet_baseline_{}_genHH_mass'.format(base_extension, args.channels[0])]
        # if args.channels[0] == "etau":
        #     channel_base_all = file['Base{}ORMETORTauOR4JetsPNet_baseline_{}_genHH_mass'.format(base_extension, args.channels[0])]
        
        # channel_base_plus_Mu50 = file['NoBaseMuMu50_baseline_{}_genHH_mass'.format(args.channels[0])]
        # channel_base_plus_Ele28 = file['NoBaseEle28_baseline_{}_genHH_mass'.format(args.channels[0])]
        # if args.channels[0] == "tautau":
        #     channel_all = file['NoBaseMuTauMETORTauORTTJet_baseline_{}_genHH_mass'.format(args.channels[0])]
        #     base = "DoubleMediumDeepTau"
        #     last_trig = "TauTauJet"
        # if args.channels[0] == "mutau":
        #     channel_all = file['NoBaseMETORTauORMu50_baseline_{}_genHH_mass'.format(args.channels[0])]
        #     base = "IsoMu + CrossMuTau"
        #     last_trig = "Mu50"
        # if args.channels[0] == "etau":
        #     channel_all = file['NoBaseMETORTauOREle28_baseline_{}_genHH_mass'.format(args.channels[0])]
        #     base = "Ele24 + CrossEleTau"
        #     last_trig = "Ele28HT"
        var_base_tot = channel_base.variances().sum()
        # var_all_tot = channel_all.variances().sum()
        # err_tot_all = np.sqrt(var_all_tot / (channel_base.values().sum() ** 2) + (channel_all.values().sum() ** 2) * var_base_tot / (channel_base.values().sum() ** 4))
        err_base = np.sqrt( (1/channel_all.values() * channel_base.errors()) ** 2 + (channel_base.values()/(channel_all.values() ** 2) * channel_all.errors()) ** 2 )
        # err_base_plus_met = np.sqrt( (1/channel_all.values() * channel_base_plus_met.errors()) ** 2 + (channel_base_plus_met.values()/(channel_all.values() ** 2) * channel_all.errors()) ** 2 )
        err_base_plus_tau = np.sqrt( (1/channel_all.values() * channel_base_plus_tau.errors()) ** 2 + (channel_base_plus_tau.values()/(channel_all.values() ** 2) * channel_all.errors()) ** 2 )
        err_base_plus_quadJetPNet = np.sqrt( (1/channel_all.values() * channel_base_plus_quadJetPNet.errors()) ** 2 + (channel_base_plus_quadJetPNet.values()/(channel_all.values() ** 2) * channel_all.errors()) ** 2 )
        # err_base_plustau_quadJet = np.sqrt( (1/channel_base.values() * channel_base_plus_quadJetPNet.errors()) ** 2 + (channel_base_plus_quadJetPNet.values()/(channel_base.values() ** 2) * channel_base.errors()) ** 2 )
        # err_base_plus_met = np.sqrt( (1/channel_base.values() * channel_base_plus_met.errors()) ** 2 + (channel_base_plus_met.values()/(channel_base.values() ** 2) * channel_base.errors()) ** 2 )
        # err_base_plus_quadJetDeep = np.sqrt( (1/channel_base.values() * channel_base_plus_quadJetDeep.errors()) ** 2 + (channel_base_plus_quadJetDeep.values()/(channel_base.values() ** 2) * channel_base.errors()) ** 2 )
        # err_base_all = np.sqrt( (1/channel_all.values() * channel_base_all.errors()) ** 2 + (channel_base_all.values()/(channel_all.values() ** 2) * channel_all.errors()) ** 2 )
        
        # err_base_plus_Mu50 = np.sqrt( (1/channel_base.values() * channel_base_plus_Mu50.errors()) ** 2 + (channel_base_plus_Mu50.values()/(channel_base.values() ** 2) * channel_base.errors()) ** 2 )
        # err_base_plus_Ele28 = np.sqrt( (1/channel_base.values() * channel_base_plus_Ele28.errors()) ** 2 + (channel_base_plus_Ele28.values()/(channel_base.values() ** 2) * channel_base.errors()) ** 2 )
        # err_all= np.sqrt( (1/channel_base.values() * channel_all.errors()) ** 2 + (channel_all.values()/(channel_base.values() ** 2) * channel_base.errors()) ** 2 )
        ax = channel_base.axes[0].centers()
        edges = channel_base.axes[0].edges()
        xerr_r = ax - edges[:-1]
        xerr_l = edges[1:] - ax
        

        hep.style.use("CMS")
        fig, ax1 = plt.subplots()
        ax2 = ax1.twinx()
        ax2.set_ylabel("Weighted MC events [a.u.]")  # we already handled the x-label with ax1
        print("genHH", ax, channel_all.values())
        ax2.bar(ax, height=channel_all.values(), width=xerr_r * 2, color='lightgrey', alpha=0.35, zorder=0, label="HHbbtautau")

        hep.style.use("CMS")
        # ax1.errorbar(ax - 15, (channel_base.values() / channel_all.values()) * 100, yerr=err_base*100, xerr=[xerr_r - 15, xerr_l + 15], fmt='o', label="Base", color=petroff6[0], zorder=5)
        # ax1.errorbar(ax - 10, (channel_base_plus_met.values() / channel_all.values()) * 100, yerr=err_base_plus_met*100, xerr=[xerr_r -10, xerr_l + 10], fmt='o', label="Base + MET", color=petroff6[1], zorder=5)
        # if args.channels[0] == "tautau":
        #     ax1.errorbar(ax, (channel_base_plus_met.values() / channel_base.values()) * 100, yerr=err_base_plus_met*100, xerr=[xerr_r + 10, xerr_l - 10], fmt='o', label="Base + TauTauJet", color=petroff6[2], zorder=5)

        #     ax1.errorbar(ax + 10, (channel_base_plus_quadJetDeep.values() / channel_base.values()) * 100, yerr=err_base_plus_quadJetDeep*100, xerr=[xerr_r + 10, xerr_l - 10], fmt='o', label="Base + QuadJetDeep", color=petroff6[3], zorder=5)

        #     ax1.errorbar(ax + 15, (channel_base_plus_quadJetPNet.values() / channel_base.values()) * 100, yerr=err_base_plus_quadJetPNet*100, xerr=[xerr_r + 10, xerr_l - 10], fmt='o', label="Base + QuadJetPNet", color=petroff6[4], zorder=5)
        # if args.channels[0] == "mutau":
        #     ax1.errorbar(ax + 5, (channel_base_plus_Mu50.values() / channel_base.values()) * 100, yerr=err_base_plus_Mu50*100, xerr=[xerr_r + 10, xerr_l - 10], fmt='o', label="Base + Mu50", color=petroff6[2], zorder=5)
        #     ax1.errorbar(ax + 15, (channel_base_plus_quadJetPNet.values() / channel_base.values()) * 100, yerr=err_base_plus_quadJetPNet*100, xerr=[xerr_r + 10, xerr_l - 10], fmt='o', label="Base + QuadJetPNet", color=petroff6[4], zorder=5)
        # if args.channels[0] == "etau":
        ax1.errorbar(ax - 5, channel_base.values(), yerr=channel_base.errors(), xerr=[xerr_r -5 , xerr_l +5], fmt='o', label="Base", color=petroff6[1], zorder=5)
        ax1.errorbar(ax + 5, channel_base_plus_quadJetPNet.values(), yerr=channel_base_plus_quadJetPNet.errors(), xerr=[xerr_r + 5, xerr_l - 5], fmt='o', label="Base + QuadJetPNet", color=petroff6[2], zorder=5)
        # ax1.errorbar(ax + 15, (channel_base_all.values() / channel_all.values()) * 100, yerr=err_base_all*100, xerr=[xerr_r + 15, xerr_l - 15], fmt='o', label="Base + TauTauJet\n+ QuadJetPNet", color=petroff6[5], zorder=5)

        # ax1.errorbar(ax + 10, (channel_base_all.values() / channel_all.values()) * 100, yerr=err_base_all*100, xerr=[xerr_r + 10, xerr_l - 10], fmt='o', label="Base + MET + IsoTau\n+ Mu50 + QuadJetPNet", color=petroff6[4], zorder=5)
        # ax1.errorbar(ax + 10, (channel_base_all.values() / channel_all.values()) * 100, yerr=err_base_all*100, xerr=[xerr_r + 10, xerr_l - 10], fmt='o', label="Base + MET + IsoTau\n+ QuadJetPNet", color=petroff6[4], zorder=5)

        ax1.text(0.05, 0.94, 'Channel: '+ channel_text, horizontalalignment='left', verticalalignment='bottom', transform=ax1.transAxes, size=21)
        # ax1.errorbar(ax + 20, (channel_all.values() / channel_base.values()) * 100, yerr=err_all*100, xerr=[xerr_r + 20, xerr_l - 20], fmt='o', label="Base + MET + IsoTau + {}".format( last_trig), color=petroff6[3])
        ax1.set_xlabel(r"m$_{HH}$ [GeV]")
        ax1.set_ylabel("Events/bin [a. u.]")
        lines, labels = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        plt.legend(lines + lines2, labels + labels2, loc=4,  bbox_to_anchor=(1, 0.1))
        # print("Total Gain: {} +/- {}".format(round(channel_all.values().sum()/channel_base.values().sum() * 100, 2), round(err_tot_all * 100, 2)))
        hep.cms.label(lumi="9.8", com="13.6")
        plt.savefig("distribution_{}_quadJet_preEE_10p11.png".format(args.channels[0]))

        for i in vars:
            channel_all = file['All_baseline_{}_{}'.format(args.channels[0], i)]
            channel_base_tautaujet = file['Base{}ORTauTauJetNo4JetsPNet_baseline_{}_{}'.format(base_extension, args.channels[0], i)]
            # channel_base_plus_Mu50 = file['Base{}ORMu50_baseline_{}_{}'.format(base_extension, args.channels[0], i)]
            channel_only_quadJetPNet = file['NoBase{}ORNoTauTauJet4JetsPNet_baseline_{}_{}'.format(base_extension, args.channels[0], i)]
                  
            var_base_tot = channel_base.variances().sum()
            # var_all_tot = channel_all.variances().sum()
            # err_tot_all = np.sqrt(var_all_tot / (channel_base.values().sum() ** 2) + (channel_all.values().sum() ** 2) * var_base_tot / (channel_base.values().sum() ** 4))
            err_base_tautaujet = np.sqrt( (1/channel_all.values() * channel_base_tautaujet.errors()) ** 2 + (channel_base_tautaujet.values()/(channel_all.values() ** 2) * channel_all.errors()) ** 2 )
            # err_base_plus_met = np.sqrt( (1/channel_all.values() * channel_base_plus_met.errors()) ** 2 + (channel_base_plus_met.values()/(channel_all.values() ** 2) * channel_all.errors()) ** 2 )
            # err_base_plus_tau = np.sqrt( (1/channel_all.values() * channel_base_plus_tau.errors()) ** 2 + (channel_base_plus_tau.values()/(channel_all.values() ** 2) * channel_all.errors()) ** 2 )
            err_only_quadJetPNet = np.sqrt( (1/channel_all.values() * channel_only_quadJetPNet.errors()) ** 2 + (channel_only_quadJetPNet.values()/(channel_all.values() ** 2) * channel_all.errors()) ** 2 )

            # err_base_plus_Mu50 = np.sqrt( (1/channel_base.values() * channel_base_plus_Mu50.errors()) ** 2 + (channel_base_plus_Mu50.values()/(channel_base.values() ** 2) * channel_base.errors()) ** 2 )
            # err_base_plus_Ele28 = np.sqrt( (1/channel_base.values() * channel_base_plus_Ele28.errors()) ** 2 + (channel_base_plus_Ele28.values()/(channel_base.values() ** 2) * channel_base.errors()) ** 2 )
            # err_all= np.sqrt( (1/channel_base.values() * channel_all.errors()) ** 2 + (channel_all.values()/(channel_base.values() ** 2) * channel_base.errors()) ** 2 )
            ax = channel_all.axes[0].centers()
            edges = channel_all.axes[0].edges()
            xerr_r = ax - edges[:-1]
            xerr_l = edges[1:] - ax
            print(edges, xerr_r, xerr_l)
        

            hep.style.use("CMS")
            fig, ax1 = plt.subplots()
            ax2 = ax1.twinx()
            ax2.set_ylabel("Weighted MC events [a.u.]")  # we already handled the x-label with ax1
            print(i, channel_all.values())
            ax2.bar(ax, height=channel_all.values(), width=xerr_r * 2, color='lightgrey', alpha=0.35, zorder=0, label="HHbbtautau")

            hep.style.use("CMS")
            # ax1.errorbar(ax - 15, (channel_base.values() / channel_all.values()) * 100, yerr=err_base*100, xerr=[xerr_r - 15, xerr_l + 15], fmt='o', label="Base", color=petroff6[0], zorder=5)
            ax1.errorbar(ax - xerr_r * 0.2, channel_base_tautaujet.values(), yerr=channel_base_tautaujet.errors(), xerr=[xerr_r * 0.8 , xerr_l * 1.2], fmt='o', label="Base + TauTauJet", color=petroff6[3], zorder=5)
            ax1.errorbar(ax + xerr_r * 0.2, channel_only_quadJetPNet.values(), yerr=channel_only_quadJetPNet.errors(), xerr=[xerr_r * 1.2, xerr_l * 0.8], fmt='o', label="QuadJetPNet", color=petroff6[2], zorder=5)
            # ax1.errorbar(ax + 15, (channel_base_all.values() / channel_all.values()) * 100, yerr=err_base_all*100, xerr=[xerr_r + 15, xerr_l - 15], fmt='o', label="Base + TauTauJet\n+ QuadJetPNet", color=petroff6[5], zorder=5)

            # ax1.errorbar(ax + 10, (channel_base_all.values() / channel_all.values()) * 100, yerr=err_base_all*100, xerr=[xerr_r + 10, xerr_l - 10], fmt='o', label="Base + MET + IsoTau\n+ Mu50 + QuadJetPNet", color=petroff6[4], zorder=5)
            # ax1.errorbar(ax + 10, (channel_base_all.values() / channel_all.values()) * 100, yerr=err_base_all*100, xerr=[xerr_r + 10, xerr_l - 10], fmt='o', label="Base + MET + IsoTau\n+ QuadJetPNet", color=petroff6[4], zorder=5)

            ax1.text(0.05, 0.94, 'Channel: '+ channel_text, horizontalalignment='left', verticalalignment='bottom', transform=ax1.transAxes, size=21)
            # ax1.errorbar(ax + 20, (channel_all.values() / channel_base.values()) * 100, yerr=err_all*100, xerr=[xerr_r + 20, xerr_l - 20], fmt='o', label="Base + MET + IsoTau + {}".format( last_trig), color=petroff6[3])
            ax1.set_xlabel("{}".format(i))
            ax1.set_ylabel("Events/bin [a. u.]")
            lines, labels = ax1.get_legend_handles_labels()
            lines2, labels2 = ax2.get_legend_handles_labels()
            plt.legend(lines + lines2, labels + labels2, loc=4,  bbox_to_anchor=(1, 0.1))
            # print("Total Gain: {} +/- {}".format(round(channel_all.values().sum()/channel_base.values().sum() * 100, 2), round(err_tot_all * 100, 2)))
            hep.cms.label(lumi="9.8", com="13.6")
            plt.savefig("distribution_{}_quadJet_preEE_10p20_{}.png".format(args.channels[0], i))
        # print(ratio.axes[0])
    return 
    channels = args.channels
    linear_x = [k for k in range(1,2)]
    edges_x  = [k-0.5 for k in range(1,2)] + [1.5]
    ptcuts = {chn: utils.get_ptcuts(chn, args.year) for chn in args.channels}

    nevents, errors = dd(lambda: dd(dict)), dd(lambda: dd(dict))
    ratios, eratios = dd(lambda: dd(dict)), dd(lambda: dd(dict))

    for adir in main_dir:
        dRstr = str(args.deltaR).replace('.', 'p')
        if len(args.channels) == 1:
            output_name = os.path.join(base_dir, 'trigger_gains_{}_{}_DR{}'.format(args.channels[0],
                                                                                        args.year, dRstr))
        elif len(args.channels) == 2:
            output_name = os.path.join(base_dir, 'trigger_gains_{}_{}_{}_DR{}'.format(*args.channels[:2],
                                                                                           args.year, dRstr))
        elif len(args.channels) == 3:
            output_name = os.path.join(base_dir, 'trigger_gains_all_{}_DR{}'.format(args.year, dRstr))
        if args.bigtau:
            output_name += "_BIGTAU"
        output_name += ".html"
        output_file(output_name)
        print('Saving file {}.'.format(output_name))

        for chn in channels:
            md = adir[chn]
            in_base = os.path.join(base_dir, md)

            nevents[md][chn]['base'], nevents[md][chn]['met'],  nevents[md][chn]['tau'] = [], [], []
            ratios[md][chn]['two'], ratios[md][chn]['met'], ratios[md][chn]['tau'] = [], [], []
            eratios[md][chn]['two'], eratios[md][chn]['met'], eratios[md][chn]['tau'] = [], [], []
            errors[md][chn] = []
              
            outname = get_outname( chn, [str(x) for x in args.region_cuts],
                                    [str(x) for x in ptcuts[chn]], str(args.met_turnon),
                                    args.bigtau)

            print(outname)
            with open(outname, "rb") as f:
                ahistos = pickle.load(f)

                #all regions summed
                sum_base_tot = round(ahistos["Base"]["legacy"]["baseline"].values().sum() +
                                        ahistos["Base"]["tau"]["baseline"].values().sum() +
                                        ahistos["Base"]["met"]["baseline"].values().sum())
                print(ahistos["Base"]["legacy"]["baseline"].values())
                #legacy region
                l1 = lambda x : round(x["legacy"]["baseline"].values().sum(), 2)
                sum_base     = l1(ahistos["Base"])
                sum_vbf      = l1(ahistos["VBF"])
                sum_met      = l1(ahistos["NoBaseMET"])
                sum_only_tau = l1(ahistos["NoBaseNoMETTau"])
                sum_tau      = l1(ahistos["NoBaseTau"])
                sum_basekin  = l1(ahistos["LegacyKin"])
                w2_basekin   = ahistos["METKin"]["legacy"]["baseline"].variances().sum()
                
                #MET region
                l2 = lambda x : round(x["met"]["baseline"].values().sum(), 2)
                sum_metkin = l2(ahistos["METKin"])
                w2_metkin = ahistos["METKin"]["met"]["baseline"].variances().sum()
                
                #Single Tau region
                l3 = lambda x : round(x["tau"]["baseline"].values().sum(), 2)
                sum_taukin = l3(ahistos["TauKin"])
                w2_taukin = ahistos["TauKin"]["tau"]["baseline"].variances().sum()

                #hypothetical VBF region
                sum_vbfkin   = l2(ahistos["VBFKin"]) + l3(ahistos["VBFKin"])

                nevents[md][chn]['base'].append(sum_basekin)
                nevents[md][chn]['met'].append(sum_basekin + sum_metkin)
                nevents[md][chn]['tau'].append(sum_basekin + sum_metkin + sum_taukin)

                rat_met_num = sum_basekin + sum_metkin
                rat_met_all = rat_met_num / sum_base_tot

                rat_tau_num = sum_basekin + sum_taukin
                rat_tau_all = rat_tau_num / sum_base_tot

                rat_all_num = sum_basekin + sum_taukin + sum_metkin
                rat_all = rat_all_num / sum_base_tot
                
                ratios[md][chn]['met'].append(rat_met_all)
                ratios[md][chn]['tau'].append(rat_tau_all)
                ratios[md][chn]['two'].append(rat_all)

                e_metkin = np.sqrt(w2_metkin)
                e_taukin = np.sqrt(w2_taukin)
                e_basekin = np.sqrt(w2_basekin)

                e_tau_num = np.sqrt(w2_taukin + w2_basekin)
                e_met_num = np.sqrt(w2_metkin + w2_basekin)
                e_all_num = np.sqrt(w2_metkin + w2_taukin + w2_basekin)
                print(rat_tau_num, rat_met_num, rat_all_num)

                eratios[md][chn]['tau'].append(rat_tau_all * np.sqrt(e_tau_num**2/rat_tau_num**2 + 1/sum_base_tot))
                eratios[md][chn]['met'].append(rat_met_all * np.sqrt(e_met_num**2/rat_met_num**2 + 1/sum_base_tot))
                eratios[md][chn]['two'].append(rat_all * np.sqrt(e_all_num**2/rat_all_num**2 + 1/sum_base_tot))
                errors[md][chn].append(e_all_num)

    json_name = 'data_' + chn + '_'
    json_name += ('bigtau' if args.bigtau else 'standard') + '.json'
    with open(json_name, 'w', encoding='utf-8') as json_obj:
        json_data = {"vals": {chn: nevents[adir[chn]][chn]['tau'] for chn in channels}}
        json_data.update({"errs": {chn: errors[adir[chn]][chn] for chn in channels}})
        json.dump(json_data, json_obj, ensure_ascii=False, indent=4)

    opt_points = dict(size=8)
    opt_line = dict(width=1.5)
    colors = ('green', 'blue', 'red', 'brown')
    styles = ('solid', 'dashed', 'dotdash')
    legends = {'base': 'Legacy',
               'met': 'MET', 'tau': 'Single Tau',
               'two': 'MET + Single Tau', 'vbf': 'VBF'}
     
    x_str = [str(k) for k in [0]]
    xticks = linear_x[:]
    yticks = [x for x in range(0,110,5)]
    shift_one = {'met': [-0.20, 0., 0.20],  'tau': [-0.20, -0.05, 0.1],
                 'vbf': [-0.10, 0.05, 0.20]}
    shift_both = {'met': [-0.20, 0., 0.20], 'two': [-0.20, -0.05, 0.1]}
    shift_kin = {'met': [-0.09, 0., 0.20],  'tau': [0.03, -0.05, 0.1],
                 'two': [-0.03, 0.05, 0.20], 'vbf': [0.09, 0.1, 0.25]}
     
    for adir in main_dir:
        print(adir)
        p_opt = dict(width=800, height=400, x_axis_label='x', y_axis_label='y')
        p1 = figure(title='Event number (' + pp(channels[0]) + ')', y_axis_type="linear", **p_opt)
        p2 = figure(title='Acceptance Gain (' + pp(channels[0]) + ')', **p_opt) if len(channels)==1 else figure(**p_opt)

        p1.yaxis.axis_label = 'Weighted number of events'
        p2.yaxis.axis_label = 'Trigger acceptance gain (w.r.t. trigger baseline) [%]'
        pics = (p1, p2)
        for p in pics:
            set_fig(p)

        for ichn,chn in enumerate(channels):
            md = adir[chn]
            print(md)
            
            p1.quad(top=nevents[md][chn]["base"], bottom=0,
                    left=edges_x[:-1], right=edges_x[1:],
                    legend_label=legends["base"]+(' ('+pp(chn)+')' if len(channels)>1 else ''),
                    fill_color="dodgerblue", line_color="black")
            p1.quad(top=nevents[md][chn]["met"], bottom=nevents[md][chn]["base"],
                    left=edges_x[:-1], right=edges_x[1:],
                    legend_label=legends["met"]+(' ('+pp(chn)+')' if len(channels)>1 else ''),
                    fill_color="green", line_color="black")
            p1.quad(top=nevents[md][chn]["tau"], bottom=nevents[md][chn]["met"],
                    left=edges_x[:-1], right=edges_x[1:],
                    legend_label=legends["tau"]+(' ('+pp(chn)+')' if len(channels)>1 else ''),
                    fill_color="red", line_color="black")
            print(md)
            for itd,td in enumerate(('met', 'tau', 'two')):
                p2.circle([x+shift_kin[td][ichn] for x in linear_x],
                          [(x-1)*100. for x in ratios[md][chn][td]],
                          color=colors[itd], fill_alpha=1., **opt_points)
                p2.line([x+shift_kin[td][ichn] for x in linear_x],
                        [(x-1)*100. for x in ratios[md][chn][td]],
                        color=colors[itd], line_dash=styles[ichn],
                        legend_label=legends[td]+(' ('+pp(chn)+')' if len(channels)>1 else ''), **opt_line)
                p2.multi_line(
                    [(x+shift_kin[td][ichn],x+shift_kin[td][ichn]) for x in linear_x], 
                    [((y-1)*100-(x*50.),(y-1)*100+(x*50.)) for x,y in zip(eratios[md][chn][td],ratios[md][chn][td])],
                    color=colors[itd], **opt_line)

        p1.legend.location = 'top_right'
        p2.legend.location = 'top_left'
        for p in pics:
            p.xaxis[0].ticker = xticks
            p.xgrid[0].ticker = xticks
            p.xgrid.grid_line_alpha = 0.2
            p.xgrid.grid_line_color = 'black'
            # p.yaxis[0].ticker = yticks
            # p.ygrid[0].ticker = yticks
            p.ygrid.grid_line_alpha = 0.2
            p.ygrid.grid_line_color = 'black'
             
            p.xaxis.axis_label = "m(X) [GeV]"
         
            p.xaxis.major_label_overrides = dict(zip(linear_x,x_str))
     
            p.legend.click_policy='hide'        
     
            p.output_backend = 'svg'
            #export_svg(p, filename='line_graph.svg')
         
        g = gridplot([[p] for p in pics])
        print(g)
        save(g, title=md)
        export_png(g, filename="line_graph_postEE.png")

if __name__ == '__main__':
    desc = "Produce plots of trigger gain VS resonance mass.\n"
    desc += "Uses the output of test_trigger_regions.py."
    desc += "When running on many channels, one should keep in mind each channel has different pT cuts."
    desc += "This might imply moving sub-folders (produced by the previous script) around."
    parser = argparse.ArgumentParser(description=desc, formatter_class=argparse.RawTextHelpFormatter)

    # parser.add_argument('--masses', required=True, nargs='+', type=str,
    #                     help='Resonance mass')
    parser.add_argument('--channels', required=True, nargs='+', type=str, 
                        choices=('etau', 'mutau', 'tautau'),
                        help='Select the channel over which the workflow will be run.' )
    parser.add_argument('--year', required=True, type=str, choices=('2016', '2017', '2018', '2022'),
                        help='Select the year over which the workflow will be run.' )
    parser.add_argument('--deltaR', type=float, default=0.5, help='DeltaR between the two leptons.')
    parser.add_argument('--bigtau', action='store_true',
                        help='Consider a larger single tau region, reducing the ditau one.')
    parser.add_argument('--met_turnon', required=False, type=str,  default=180,
                        help='MET trigger turnon cut [GeV].' )
    parser.add_argument('--region_cuts', required=False, type=float, nargs=2, default=(190, 190),
                        help='High/low regions pT1 and pT2 selection cuts [GeV].' )

    args = utils.parse_args(parser)

    base_dir = '/t3home/fbilandz/TriggerScaleFactors/'
    main_dir = [{"etau":   "Region_Spin2_190_190_PT_33_25_35_DR_{}_TURNON_200_190".format(args.deltaR),
                 "mutau":  "mutau_190_190_DR_0.5_PT_25_21_32_TURNON_180",
                 "tautau": "Region_Spin2_190_190_PT_40_40_DR_{}_TURNON_200_190".format(args.deltaR)},
                ]
    
    main(args)
