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
    fig.title.text_font_size = "15px"
    fig.xaxis.axis_label_text_font_style = "bold"
    fig.yaxis.axis_label_text_font_style = "bold"
    fig.xaxis.axis_label_text_font_size = "13px"
    fig.yaxis.axis_label_text_font_size = "13px"

def main(args):
    with uproot.open("data/regions_preEE_11p4-jet-higher-pt_all.root".format(args.channels[0])) as file:
        # from matplotlib.colors import ListedColormap
        petroff6 = ["#5790fc", "#f89c20", "#e42536", "#964a8b", "#9c9ca1", "#7a21dd"]
        channel_base = file['BaseTau_baseline_{}_genHH_mass'.format(args.channels[0])]
        # channel_inclusive = file['NoBaseTauMETTau_baseline_{}_genHH_mass'.format(args.channels[0])]
        # channel_base_plus_met = file['NoBaseTauMETNoTau_baseline_{}_genHH_mass'.format(args.channels[0])]
        # channel_base_plus_tau = file['NoBaseTauNoMETTau_baseline_{}_genHH_mass'.format(args.channels[0])]

        # err_inclusive = np.sqrt( (1/channel_base.values() * channel_inclusive.errors()) ** 2 + (channel_inclusive.values()/(channel_base.values() ** 2) * channel_base.errors()) ** 2 )
        # err_base_plus_met = np.sqrt( (1/channel_base.values() * channel_base_plus_met.errors()) ** 2 + (channel_base_plus_met.values()/(channel_base.values() ** 2) * channel_base.errors()) ** 2 )
        # err_base_plus_tau = np.sqrt( (1/channel_base.values() * channel_base_plus_tau.errors()) ** 2 + (channel_base_plus_tau.values()/(channel_base.values() ** 2) * channel_base.errors()) ** 2 )
        # sum_base = channel_base.values().sum()
        # err_tot_base = channel_base.variances().sum()
        # sum_inclusive = channel_inclusive.values().sum()
        # err_tot_inclusive = channel_inclusive.variances().sum()
        # sum_met = channel_base_plus_met.values().sum()
        # err_tot_omet = channel_base_plus_met.variances().sum()
        # sum_tau = channel_base_plus_tau.values().sum()
        # err_tot_otau = channel_base_plus_tau.variances().sum()
        # err_tot_met = np.sqrt((err_tot_inclusive + err_tot_omet)/(sum_base**2) + (sum_inclusive + sum_met) ** 2 * err_tot_base / (sum_base ** 4))
        # err_tot_tau = np.sqrt((err_tot_inclusive + err_tot_otau)/(sum_base**2) + (sum_inclusive + sum_tau) ** 2 * err_tot_base /(sum_base ** 4))
        # err_tot_met = round(err_tot_met * 100, 2)
        # err_tot_tau = round(err_tot_tau * 100, 2)
        # sum_overlap = sum_inclusive/(sum_inclusive + sum_met + sum_tau)
        # sum_overlap = round(sum_overlap * 100, 2)
        # err_overlap = np.sqrt((err_tot_inclusive)/((sum_inclusive + sum_met + sum_tau) ** 2) + (sum_inclusive ** 2)/((sum_inclusive + sum_met + sum_tau) ** 4)*(err_tot_inclusive + err_tot_omet + err_tot_otau))
        # err_overlap = round(err_overlap * 100, 2)

        # err_base_plus_ttjet = np.sqrt( (1/channel_base.values() * channel_base_plus_ttjet.errors()) ** 2 + (channel_base_plus_ttjet.values()/(channel_base.values() ** 2) * channel_base.errors()) ** 2 )
        # err_base_plus_Mu50 = np.sqrt( (1/channel_base.values() * channel_base_plus_Mu50.errors()) ** 2 + (channel_base_plus_Mu50.values()/(channel_base.values() ** 2) * channel_base.errors()) ** 2 )
        # err_base_plus_Ele28 = np.sqrt( (1/channel_base.values() * channel_base_plus_Ele28.errors()) ** 2 + (channel_base_plus_Ele28.values()/(channel_base.values() ** 2) * channel_base.errors()) ** 2 )
        # err_all = np.sqrt( (1/channel_base.values() * channel_all.errors()) ** 2 + (channel_all.values()/(channel_base.values() ** 2) * channel_base.errors()) ** 2 )
        ax = channel_base.axes[0].centers()
        edges = channel_base.axes[0].edges()
        xerr_r = ax - edges[:-1]
        xerr_l = edges[1:] - ax
        print(edges, xerr_r, xerr_l)
        print(channel_base.errors()/channel_base.values())
        hep.style.use("CMS")

        plt.figure(0)
        # plt.bar(x=ax, height=(channel_base_plus_met.values()/channel_base.values()) * 100, bottom= (channel_inclusive.values()/channel_base.values()) * 100,width=xerr_r * 2, yerr=err_base_plus_met * 100, color="#fff", edgecolor=petroff6[0], ecolor=petroff6[0], label="MET & !Tau")
        
        # plt.bar(x=ax, height=-(channel_base_plus_tau.values()/channel_base.values()) * 100, bottom=-(channel_inclusive.values()/channel_base.values()) * 100, width=xerr_r * 2, yerr=err_base_plus_tau * 100, label="Tau & !MET", color="#fff", edgecolor=petroff6[1], ecolor=petroff6[1])
        # plt.bar(x=ax, height=(channel_inclusive.values()/channel_base.values()) * 100, width=xerr_r * 2, yerr=err_inclusive * 100, color="#fff", edgecolor=petroff6[2], ecolor=petroff6[2], label="MET & Tau")
        # plt.bar(x=ax, height=-(channel_inclusive.values()/channel_base.values()) * 100, width=xerr_r * 2, yerr=err_inclusive * 100, color="#fff", edgecolor=petroff6[2], ecolor=petroff6[2])
        # plt.axhline(y=0, color='black', linestyle='-')
        
        # plt.text(220, 24, "Total MET Gain: {} +/- {}%".format(round((sum_inclusive + sum_met)/sum_base * 100, 2), err_tot_met))
        # plt.text(250, -30.8, "Total Tau Gain: {} +/- {}%".format(round((sum_inclusive + sum_tau)/sum_base * 100, 2), err_tot_tau))
        # print("Total MET Gain: {} +/- {}%".format(round((sum_inclusive + sum_met)/sum_base * 100, 2), err_tot_met))
        # print("Total Tau Gain: {} +/- {}%".format(round((sum_inclusive + sum_tau)/sum_base * 100, 2), err_tot_tau))
        # print(sum_inclusive, sum_met, sum_tau)
        # print("Total MET/Tau overlap: {} +/- {}%".format(sum_overlap, err_overlap))
        # ticks_locs, ticks_labels = plt.yticks()
        # print(ticks_locs, ticks_labels)
# set labels to absolute values and with integer representation
        # plt.yticks(ticks_locs, [int(abs(tick)) for tick in ticks_locs])

        # plt.errorbar(ax - 10, (channel_base_plus_met.values() / channel_base.values()) * 100, yerr=err_base_plus_met*100, xerr=[xerr_r - 10, xerr_l + 10], fmt='o', label="{} + MET".format(base), color=petroff6[0])
        # plt.errorbar(ax, (channel_base_plus_tau.values() / channel_base.values()) * 100, yerr=err_base_plus_tau*100, xerr=[xerr_r, xerr_l], fmt='o', label="{} + IsoTau".format(base), color=petroff6[1])
        # if args.channels[0] == "tautau":
        #     plt.errorbar(ax + 10, (channel_base_plus_Tau.values() / channel_base.values()) * 100, yerr=err_base_plus_ttjet*100, xerr=[xerr_r + 10, xerr_l - 10], fmt='o', label="{} + TauTauJet".format(base), color=petroff6[2])
        # if args.channels[0] == "mutau":
        #     plt.errorbar(ax + 10, (channel_base_plus_Mu50.values() / channel_base.values()) * 100, yerr=err_base_plus_Mu50*100, xerr=[xerr_r + 10, xerr_l - 10], fmt='o', label="{} + Mu50".format(base), color=petroff6[2])
        # if args.channels[0] == "etau":
        #     plt.errorbar(ax + 10, (channel_base_plus_Ele28.values() / channel_base.values()) * 100, yerr=err_base_plus_Ele28*100, xerr=[xerr_r + 10, xerr_l - 10], fmt='o', label="{} + Ele28".format(base), color=petroff6[2])

        # plt.errorbar(ax + 20, (channel_all.values() / channel_base.values()) * 100, yerr=err_all*100, xerr=[xerr_r + 20, xerr_l - 20], fmt='o', label="{} + MET + IsoTau + {}".format(base, last_trig), color=petroff6[3])
        # plt.xlabel(r"m$_{HH}$ [GeV]")
        # plt.ylabel("Gain [%]")
        # plt.legend()
        # hep.cms.label(lumi="9.8", com="13.6")
        # plt.savefig("orthogonal_gains_{}.png".format(args.channels[0]))

        channel_inclusive = file['NoBaseTauTTJet4JetsPNet_baseline_{}_genHH_mass'.format(args.channels[0])]
        channel_base_plus_ttjet = file['NoBaseTauTauTauJetNo4JetsPNet_baseline_{}_genHH_mass'.format(args.channels[0])]
        channel_base_plus_QuadJetPNet = file['NoBaseTauNoTauTauJet4JetsPNet_baseline_{}_genHH_mass'.format(args.channels[0])]

        err_inclusive = np.sqrt( (1/channel_base.values() * channel_inclusive.errors()) ** 2 + (channel_inclusive.values()/(channel_base.values() ** 2) * channel_base.errors()) ** 2 )
        err_base_plus_ttjet = np.sqrt( (1/channel_base.values() * channel_base_plus_ttjet.errors()) ** 2 + (channel_base_plus_ttjet.values()/(channel_base.values() ** 2) * channel_base.errors()) ** 2 )
        err_base_plus_QuadJetPNet = np.sqrt( (1/channel_base.values() * channel_base_plus_QuadJetPNet.errors()) ** 2 + (channel_base_plus_QuadJetPNet.values()/(channel_base.values() ** 2) * channel_base.errors()) ** 2 )
        sum_base = channel_base.values().sum()
        err_tot_base = channel_base.variances().sum()
        sum_inclusive = channel_inclusive.values().sum()
        err_tot_inclusive = channel_inclusive.variances().sum()
        sum_ttjet = channel_base_plus_ttjet.values().sum()
        err_tot_ottjet = channel_base_plus_ttjet.variances().sum()
        sum_QuadJetPNet = channel_base_plus_QuadJetPNet.values().sum()
        err_tot_oQuadJetPNet = channel_base_plus_QuadJetPNet.variances().sum()
        err_tot_ttjet = np.sqrt((err_tot_inclusive + err_tot_ottjet)/(sum_base**2) + (sum_inclusive + sum_ttjet) ** 2 * err_tot_base / (sum_base ** 4))
        err_tot_QuadJetPNet = np.sqrt((err_tot_inclusive + err_tot_oQuadJetPNet)/(sum_base**2) + (sum_inclusive + sum_QuadJetPNet) ** 2 * err_tot_base /(sum_base ** 4))
        err_tot_ttjet = round(err_tot_ttjet * 100, 2)
        err_tot_QuadJetPNet = round(err_tot_QuadJetPNet * 100, 2)
        sum_overlap = sum_inclusive/(sum_inclusive + sum_ttjet + sum_QuadJetPNet)
        sum_overlap = round(sum_overlap * 100, 2)
        err_overlap = np.sqrt((err_tot_inclusive)/((sum_inclusive + sum_ttjet + sum_QuadJetPNet) ** 2) + (sum_inclusive ** 2)/((sum_inclusive + sum_ttjet + sum_QuadJetPNet) ** 4)*(err_tot_inclusive + err_tot_ottjet + err_tot_oQuadJetPNet))
        err_overlap = round(err_overlap * 100, 2)


        plt.figure(1)
        plt.bar(x=ax, height=(channel_base_plus_ttjet.values()/channel_base.values()) * 100, bottom= (channel_inclusive.values()/channel_base.values()) * 100,width=xerr_r * 2, yerr=err_base_plus_ttjet * 100, color="#fff", edgecolor=petroff6[0], ecolor=petroff6[0], label="TauTauJet & !Tau")
        
        plt.bar(x=ax, height=-(channel_base_plus_QuadJetPNet.values()/channel_base.values()) * 100, bottom=-(channel_inclusive.values()/channel_base.values()) * 100, width=xerr_r * 2, yerr=err_base_plus_QuadJetPNet * 100, label="QuadJetPNet & !TauTauJet", color="#fff", edgecolor=petroff6[1], ecolor=petroff6[1])
        plt.bar(x=ax, height=(channel_inclusive.values()/channel_base.values()) * 100, width=xerr_r * 2, yerr=err_inclusive * 100, color="#fff", edgecolor=petroff6[2], ecolor=petroff6[2], label="TauTauJet & QuadJetPNet")
        plt.bar(x=ax, height=-(channel_inclusive.values()/channel_base.values()) * 100, width=xerr_r * 2, yerr=err_inclusive * 100, color="#fff", edgecolor=petroff6[2], ecolor=petroff6[2])
        plt.axhline(y=0, color='black', linestyle='-')
        
        # plt.text(220, 24, "Total ttjet Gain: {} +/- {}%".format(round((sum_inclusive + sum_ttjet)/sum_base * 100, 2), err_tot_ttjet))
        # plt.text(250, -30.8, "Total QuadJetPNet Gain: {} +/- {}%".format(round((sum_inclusive + sum_QuadJetPNet)/sum_base * 100, 2), err_tot_QuadJetPNet))
        print("Total ttjet Gain: {} +/- {}%".format(round((sum_inclusive + sum_ttjet)/sum_base * 100, 2), err_tot_ttjet))
        print("Total QuadJetPNet Gain: {} +/- {}%".format(round((sum_inclusive + sum_QuadJetPNet)/sum_base * 100, 2), err_tot_QuadJetPNet))
        print(sum_inclusive, sum_ttjet, sum_QuadJetPNet)
        print("Total ttjet/QuadJetPNet overlap: {} +/- {}%".format(sum_overlap, err_overlap))
        ticks_locs, ticks_labels = plt.yticks()
        print(ticks_locs, ticks_labels)
# set labels to absolute values and with integer representation
        plt.yticks(ticks_locs, [int(abs(tick)) for tick in ticks_locs])

        plt.xlabel(r"m$_{HH}$ [GeV]")
        plt.ylabel("Gain [%]")
        plt.legend()
        hep.cms.label(lumi="27.8", com="13.6")
        plt.savefig("orthogonal_gains_{}_QuadJetPNet_preEE_11p4_withOffline_pf75.png".format(args.channels[0]))

        # if args.channels[0] == "tautau":
        #     channel_inclusive = file['NoBaseMETTTJet_baseline_{}_genHH_mass'.format(args.channels[0])]
        #     channel_base_plus_met = file['NoBaseMETNoTTJet_baseline_{}_genHH_mass'.format(args.channels[0])]
        #     channel_base_plus_TTJet = file['NoBaseNoMETTTJet_baseline_{}_genHH_mass'.format(args.channels[0])]
        #     err_inclusive = np.sqrt( (1/channel_base.values() * channel_inclusive.errors()) ** 2 + (channel_inclusive.values()/(channel_base.values() ** 2) * channel_base.errors()) ** 2 )
        #     err_base_plus_met = np.sqrt( (1/channel_base.values() * channel_base_plus_met.errors()) ** 2 + (channel_base_plus_met.values()/(channel_base.values() ** 2) * channel_base.errors()) ** 2 )
        #     err_base_plus_TTJet = np.sqrt( (1/channel_base.values() * channel_base_plus_TTJet.errors()) ** 2 + (channel_base_plus_TTJet.values()/(channel_base.values() ** 2) * channel_base.errors()) ** 2 )
        #     sum_base = channel_base.values().sum()
        #     err_tot_base = channel_base.variances().sum()
        #     sum_inclusive = channel_inclusive.values().sum()
        #     err_tot_inclusive = channel_inclusive.variances().sum()
        #     sum_met = channel_base_plus_met.values().sum()
        #     err_tot_omet = channel_base_plus_met.variances().sum()
        #     sum_TTJet = channel_base_plus_TTJet.values().sum()
        #     err_tot_oTTJet = channel_base_plus_TTJet.variances().sum()
        #     err_tot_met = np.sqrt((err_tot_inclusive + err_tot_omet)/(sum_base**2) + (sum_inclusive + sum_met) ** 2 * err_tot_base / (sum_base ** 4))
        #     err_tot_TTJet = np.sqrt((err_tot_inclusive + err_tot_oTTJet)/(sum_base**2) + (sum_inclusive + sum_TTJet) ** 2 * err_tot_base /(sum_base ** 4))
        #     err_tot_met = round(err_tot_met * 100, 2)
        #     err_tot_TTJet = round(err_tot_TTJet * 100, 2)
        #     sum_overlap = sum_inclusive/(sum_inclusive + sum_met + sum_TTJet)
        #     sum_overlap = round(sum_overlap * 100, 2)
        #     err_overlap = np.sqrt((err_tot_inclusive)/((sum_inclusive + sum_met + sum_TTJet) ** 2) + (sum_inclusive ** 2)/((sum_inclusive + sum_met + sum_TTJet) ** 4)*(err_tot_inclusive + err_tot_omet + err_tot_oTTJet))
        #     err_overlap = round(err_overlap * 100, 2)

        #     plt.figure(1)
        #     plt.bar(x=ax, height=(channel_base_plus_met.values()/channel_base.values()) * 100, bottom= (channel_inclusive.values()/channel_base.values()) * 100,width=xerr_r * 2, yerr=err_base_plus_met * 100, color="#fff", edgecolor=petroff6[0], ecolor=petroff6[0], label="MET & !TTJet")
            
        #     plt.bar(x=ax, height=-(channel_base_plus_TTJet.values()/channel_base.values()) * 100, bottom=-(channel_inclusive.values()/channel_base.values()) * 100, width=xerr_r * 2, yerr=err_base_plus_TTJet * 100, label="TTJet & !MET", color="#fff", edgecolor=petroff6[1], ecolor=petroff6[1])
        #     plt.bar(x=ax, height=(channel_inclusive.values()/channel_base.values()) * 100, width=xerr_r * 2, yerr=err_inclusive * 100, color="#fff", edgecolor=petroff6[2], ecolor=petroff6[2], label="MET & TTJet")
        #     plt.bar(x=ax, height=-(channel_inclusive.values()/channel_base.values()) * 100, width=xerr_r * 2, yerr=err_inclusive * 100, color="#fff", edgecolor=petroff6[2], ecolor=petroff6[2])
        #     plt.axhline(y=0, color='black', linestyle='-')
        #     # plt.text(250, -30.8, "Total Tau Gain: {} +/- {}%".format(round((sum_inclusive + sum_tau)/sum_base * 100, 2), err_tot_tau))
        #     print("Total MET Gain: {} +/- {}%".format(round((sum_inclusive + sum_met)/sum_base * 100, 2), err_tot_met))
        #     print("Total TTJet Gain: {} +/- {}%".format(round((sum_inclusive + sum_TTJet)/sum_base * 100, 2), err_tot_TTJet))
        #     print("Total MET/TTJet overlap: {} +/- {}%".format(sum_overlap, err_overlap))
        #     ticks_locs, ticks_labels = plt.yticks()
        #     print(ticks_locs, ticks_labels)
        #     plt.xlabel(r"m$_{HH}$ [GeV]")
        #     plt.ylabel("Gain [%]")
        #     plt.legend()
        #     hep.cms.label(lumi="9.8", com="13.6")
        #     plt.savefig("orthogonal_gains_{}_ttjet.png".format(args.channels[0]))

        #     channel_inclusive = file['NoBaseTauTTJet_baseline_{}_genHH_mass'.format(args.channels[0])]
        #     channel_base_plus_tau = file['NoBaseTauNoTTJet_baseline_{}_genHH_mass'.format(args.channels[0])]
        #     channel_base_plus_TTJet = file['NoBaseNoTauTTJet_baseline_{}_genHH_mass'.format(args.channels[0])]
        #     err_inclusive = np.sqrt( (1/channel_base.values() * channel_inclusive.errors()) ** 2 + (channel_inclusive.values()/(channel_base.values() ** 2) * channel_base.errors()) ** 2 )
        #     err_base_plus_tau = np.sqrt( (1/channel_base.values() * channel_base_plus_tau.errors()) ** 2 + (channel_base_plus_tau.values()/(channel_base.values() ** 2) * channel_base.errors()) ** 2 )
        #     err_base_plus_TTJet = np.sqrt( (1/channel_base.values() * channel_base_plus_TTJet.errors()) ** 2 + (channel_base_plus_TTJet.values()/(channel_base.values() ** 2) * channel_base.errors()) ** 2 )
        #     sum_base = channel_base.values().sum()
        #     err_tot_base = channel_base.variances().sum()
        #     sum_inclusive = channel_inclusive.values().sum()
        #     err_tot_inclusive = channel_inclusive.variances().sum()
        #     sum_tau = channel_base_plus_tau.values().sum()
        #     err_tot_otau = channel_base_plus_tau.variances().sum()
        #     sum_TTJet = channel_base_plus_TTJet.values().sum()
        #     err_tot_oTTJet = channel_base_plus_TTJet.variances().sum()
        #     err_tot_tau = np.sqrt((err_tot_inclusive + err_tot_otau)/(sum_base**2) + (sum_inclusive + sum_tau) ** 2 * err_tot_base / (sum_base ** 4))
        #     err_tot_TTJet = np.sqrt((err_tot_inclusive + err_tot_oTTJet)/(sum_base**2) + (sum_inclusive + sum_TTJet) ** 2 * err_tot_base /(sum_base ** 4))
        #     err_tot_tau = round(err_tot_tau * 100, 2)
        #     err_tot_TTJet = round(err_tot_TTJet * 100, 2)
        #     sum_overlap = sum_inclusive/(sum_inclusive + sum_tau + sum_TTJet)
        #     sum_overlap = round(sum_overlap * 100, 2)
        #     err_overlap = np.sqrt((err_tot_inclusive)/((sum_inclusive + sum_tau + sum_TTJet) ** 2) + (sum_inclusive ** 2)/((sum_inclusive + sum_tau + sum_TTJet) ** 4)*(err_tot_inclusive + err_tot_otau + err_tot_oTTJet))
        #     err_overlap = round(err_overlap * 100, 2)

        #     plt.figure(2)
        #     plt.bar(x=ax, height=(channel_base_plus_tau.values()/channel_base.values()) * 100, bottom= (channel_inclusive.values()/channel_base.values()) * 100,width=xerr_r * 2, yerr=err_base_plus_tau * 100, color="#fff", edgecolor=petroff6[0], ecolor=petroff6[0], label="Tau & !TauTauJet")
            
        #     plt.bar(x=ax, height=-(channel_base_plus_TTJet.values()/channel_base.values()) * 100, bottom=-(channel_inclusive.values()/channel_base.values()) * 100, width=xerr_r * 2, yerr=err_base_plus_TTJet * 100, label="TauTauJet & !Tau", color="#fff", edgecolor=petroff6[1], ecolor=petroff6[1])
        #     plt.bar(x=ax, height=(channel_inclusive.values()/channel_base.values()) * 100, width=xerr_r * 2, yerr=err_inclusive * 100, color="#fff", edgecolor=petroff6[2], ecolor=petroff6[2], label="Tau & TauTauJet")
        #     plt.bar(x=ax, height=-(channel_inclusive.values()/channel_base.values()) * 100, width=xerr_r * 2, yerr=err_inclusive * 100, color="#fff", edgecolor=petroff6[2], ecolor=petroff6[2])
        #     plt.axhline(y=0, color='black', linestyle='-')
        #     # plt.text(250, -30.8, "Total Tau Gain: {} +/- {}%".format(round((sum_inclusive + sum_tau)/sum_base * 100, 2), err_tot_tau))
        #     print("Total Tau Gain: {} +/- {}%".format(round((sum_inclusive + sum_tau)/sum_base * 100, 2), err_tot_met))
        #     print("Total TTJet Gain: {} +/- {}%".format(round((sum_inclusive + sum_TTJet)/sum_base * 100, 2), err_tot_TTJet))
        #     print("Total Tau/TTJet overlap: {} +/- {}%".format(sum_overlap, err_overlap))
        #     ticks_locs, ticks_labels = plt.yticks()
        #     print(ticks_locs, ticks_labels)
        #     plt.xlabel(r"m$_{HH}$ [GeV]")
        #     plt.ylabel("Gain [%]")
        #     plt.legend()
        #     hep.cms.label(lumi="9.8", com="13.6")
        #     plt.savefig("orthogonal_gains_{}_tau_ttjet.png".format(args.channels[0]))
        # # print(ratio.axes[0])

    #     # elif args.channels[0] == "mutau":
    #         channel_inclusive = file['NoBaseMETMu50_baseline_{}_genHH_mass'.format(args.channels[0])]
    #         channel_base_plus_met = file['NoBaseMETNoMu50_baseline_{}_genHH_mass'.format(args.channels[0])]
    #         channel_base_plus_Mu50 = file['NoBaseNoMETMu50_baseline_{}'.format(args.channels[0])]
    #         err_inclusive = np.sqrt( (1/channel_base.values() * channel_inclusive.errors()) ** 2 + (channel_inclusive.values()/(channel_base.values() ** 2) * channel_base.errors()) ** 2 )
    #         err_base_plus_met = np.sqrt( (1/channel_base.values() * channel_base_plus_met.errors()) ** 2 + (channel_base_plus_met.values()/(channel_base.values() ** 2) * channel_base.errors()) ** 2 )
    #         err_base_plus_Mu50 = np.sqrt( (1/channel_base.values() * channel_base_plus_Mu50.errors()) ** 2 + (channel_base_plus_Mu50.values()/(channel_base.values() ** 2) * channel_base.errors()) ** 2 )
    #         sum_base = channel_base.values().sum()
    #         err_tot_base = channel_base.variances().sum()
    #         sum_inclusive = channel_inclusive.values().sum()
    #         err_tot_inclusive = channel_inclusive.variances().sum()
    #         sum_met = channel_base_plus_met.values().sum()
    #         err_tot_omet = channel_base_plus_met.variances().sum()
    #         sum_Mu50 = channel_base_plus_Mu50.values().sum()
    #         err_tot_oMu50 = channel_base_plus_Mu50.variances().sum()
    #         err_tot_met = np.sqrt((err_tot_inclusive + err_tot_omet)/(sum_base**2) + (sum_inclusive + sum_met) ** 2 * err_tot_base / (sum_base ** 4))
    #         err_tot_Mu50 = np.sqrt((err_tot_inclusive + err_tot_oMu50)/(sum_base**2) + (sum_inclusive + sum_Mu50) ** 2 * err_tot_base /(sum_base ** 4))
    #         sum_overlap = sum_inclusive/(sum_inclusive + sum_met + sum_Mu50)
    #         sum_overlap = round(sum_overlap * 100, 2)
    #         err_overlap = np.sqrt((err_tot_inclusive)/((sum_inclusive + sum_met + sum_Mu50) ** 2) + (sum_inclusive ** 2)/((sum_inclusive + sum_met + sum_Mu50) ** 4)*(err_tot_inclusive + err_tot_omet + err_tot_oMu50))
    #         err_tot_met = round(err_tot_met * 100, 2)
    #         err_tot_Mu50 = round(err_tot_Mu50 * 100, 2)
    #         err_overlap = round(err_overlap * 100, 2)
        
    #         plt.figure(1)
    #         plt.bar(x=ax, height=(channel_base_plus_met.values()/channel_base.values()) * 100, bottom= (channel_inclusive.values()/channel_base.values()) * 100,width=xerr_r * 2, yerr=err_base_plus_met * 100, color="#fff", edgecolor=petroff6[0], ecolor=petroff6[0], label="MET & !Mu50")
            
    #         plt.bar(x=ax, height=-(channel_base_plus_Mu50.values()/channel_base.values()) * 100, bottom=-(channel_inclusive.values()/channel_base.values()) * 100, width=xerr_r * 2, yerr=err_base_plus_Mu50 * 100, label="Mu50 & !MET", color="#fff", edgecolor=petroff6[1], ecolor=petroff6[1])
    #         plt.bar(x=ax, height=(channel_inclusive.values()/channel_base.values()) * 100, width=xerr_r * 2, yerr=err_inclusive * 100, color="#fff", edgecolor=petroff6[2], ecolor=petroff6[2], label="MET & Mu50")
    #         plt.bar(x=ax, height=-(channel_inclusive.values()/channel_base.values()) * 100, width=xerr_r * 2, yerr=err_inclusive * 100, color="#fff", edgecolor=petroff6[2], ecolor=petroff6[2])
    #         plt.axhline(y=0, color='black', linestyle='-')
    #         # plt.text(250, -30.8, "Total Tau Gain: {} +/- {}%".format(round((sum_inclusive + sum_tau)/sum_base * 100, 2), err_tot_tau))
            

    #         print("Total MET Gain: {} +/- {}%".format(round((sum_inclusive + sum_met)/sum_base * 100, 2), err_tot_met))
    #         print("Total Mu50 Gain: {} +/- {}%".format(round((sum_inclusive + sum_Mu50)/sum_base * 100, 2), err_tot_Mu50))
    #         print("Total MET/Mu50 overlap: {} +/- {}%".format(sum_overlap, err_overlap))
    #         ticks_locs, ticks_labels = plt.yticks()
    #         print(ticks_locs, ticks_labels)
    #         plt.xlabel(r"m$_{HH}$ [GeV]")
    #         plt.ylabel("Gain [%]")
    #         plt.legend()
    #         hep.cms.label(lumi="9.8", com="13.6")
    #         plt.savefig("orthogonal_gains_{}_Mu50.png".format(args.channels[0]))

    #         channel_inclusive = file['NoBaseTauMu50_baseline_{}'.format(args.channels[0])]
    #         channel_base_plus_tau = file['NoBaseTauNoMu50_baseline_{}'.format(args.channels[0])]
    #         channel_base_plus_Mu50 = file['NoBaseNoTauMu50_baseline_{}'.format(args.channels[0])]
    #         err_inclusive = np.sqrt( (1/channel_base.values() * channel_inclusive.errors()) ** 2 + (channel_inclusive.values()/(channel_base.values() ** 2) * channel_base.errors()) ** 2 )
    #         err_base_plus_tau = np.sqrt( (1/channel_base.values() * channel_base_plus_tau.errors()) ** 2 + (channel_base_plus_tau.values()/(channel_base.values() ** 2) * channel_base.errors()) ** 2 )
    #         err_base_plus_Mu50 = np.sqrt( (1/channel_base.values() * channel_base_plus_Mu50.errors()) ** 2 + (channel_base_plus_Mu50.values()/(channel_base.values() ** 2) * channel_base.errors()) ** 2 )
    #         sum_base = channel_base.values().sum()
    #         err_tot_base = channel_base.variances().sum()
    #         sum_inclusive = channel_inclusive.values().sum()
    #         err_tot_inclusive = channel_inclusive.variances().sum()
    #         sum_tau = channel_base_plus_tau.values().sum()
    #         err_tot_otau = channel_base_plus_tau.variances().sum()
    #         sum_Mu50 = channel_base_plus_Mu50.values().sum()
    #         err_tot_oMu50 = channel_base_plus_Mu50.variances().sum()
    #         err_tot_tau = np.sqrt((err_tot_inclusive + err_tot_otau)/(sum_base**2) + (sum_inclusive + sum_tau) ** 2 * err_tot_base / (sum_base ** 4))
    #         err_tot_Mu50 = np.sqrt((err_tot_inclusive + err_tot_oMu50)/(sum_base**2) + (sum_inclusive + sum_Mu50) ** 2 * err_tot_base /(sum_base ** 4))
    #         err_tot_tau = round(err_tot_tau * 100, 2)
    #         err_tot_Mu50 = round(err_tot_Mu50 * 100, 2)
    #         sum_overlap = sum_inclusive/(sum_inclusive + sum_tau + sum_Mu50)
    #         sum_overlap = round(sum_overlap * 100, 2)
    #         err_overlap = np.sqrt((err_tot_inclusive)/((sum_inclusive + sum_tau + sum_Mu50) ** 2) + (sum_inclusive ** 2)/((sum_inclusive + sum_tau + sum_Mu50) ** 4)*(err_tot_inclusive + err_tot_otau + err_tot_oMu50))
    #         err_overlap = round(err_overlap * 100, 2)

    #         plt.figure(2)
    #         plt.bar(x=ax, height=(channel_base_plus_tau.values()/channel_base.values()) * 100, bottom= (channel_inclusive.values()/channel_base.values()) * 100,width=xerr_r * 2, yerr=err_base_plus_tau * 100, color="#fff", edgecolor=petroff6[0], ecolor=petroff6[0], label="Tau & !Mu50")
            
    #         plt.bar(x=ax, height=-(channel_base_plus_Mu50.values()/channel_base.values()) * 100, bottom=-(channel_inclusive.values()/channel_base.values()) * 100, width=xerr_r * 2, yerr=err_base_plus_Mu50 * 100, label="Mu50 & !Tau", color="#fff", edgecolor=petroff6[1], ecolor=petroff6[1])
    #         plt.bar(x=ax, height=(channel_inclusive.values()/channel_base.values()) * 100, width=xerr_r * 2, yerr=err_inclusive * 100, color="#fff", edgecolor=petroff6[2], ecolor=petroff6[2], label="Tau & Mu50")
    #         plt.bar(x=ax, height=-(channel_inclusive.values()/channel_base.values()) * 100, width=xerr_r * 2, yerr=err_inclusive * 100, color="#fff", edgecolor=petroff6[2], ecolor=petroff6[2])
    #         plt.axhline(y=0, color='black', linestyle='-')
    #         # plt.text(250, -30.8, "Total Tau Gain: {} +/- {}%".format(round((sum_inclusive + sum_tau)/sum_base * 100, 2), err_tot_tau))
    #         print("Total Tau Gain: {} +/- {}%".format(round((sum_inclusive + sum_tau)/sum_base * 100, 2), err_tot_met))
    #         print("Total Mu50 Gain: {} +/- {}%".format(round((sum_inclusive + sum_Mu50)/sum_base * 100, 2), err_tot_Mu50))
    #         print("Total Tau/Mu50 overlap: {} +/- {}%".format(sum_overlap, err_overlap))
    #         ticks_locs, ticks_labels = plt.yticks()
    #         print(ticks_locs, ticks_labels)
    #         plt.xlabel(r"m$_{HH}$ [GeV]")
    #         plt.ylabel("Gain [%]")
    #         plt.legend()
    #         hep.cms.label(lumi="9.8", com="13.6")
    #         plt.savefig("orthogonal_gains_{}_tau_Mu50.png".format(args.channels[0]))
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
    shift_one = {'met': [-0.15, 0., 0.15],  'tau': [-0.20, -0.05, 0.1],
                 'vbf': [-0.10, 0.05, 0.20]}
    shift_both = {'met': [-0.15, 0., 0.15], 'two': [-0.20, -0.05, 0.1]}
    shift_kin = {'met': [-0.09, 0., 0.15],  'tau': [0.03, -0.05, 0.1],
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
        export_png(g, filename="line_graph.png")

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
