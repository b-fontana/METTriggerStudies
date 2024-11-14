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
    vars = ['dau1_pt', 
            'dau2_pt', 'dau1_eta', 'dau2_eta',
            'dau1_tauIdVSjet', 'dau2_tauIdVSjet',
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
    genHH_mass_baseline = uproot.open("data/regions_postBPix_15p1-no-trigger_all.root")['tautau_genHH_mass']
    genHH_mass_base_tau = uproot.open("data/regions_2023_postBPix_15p35-base-tau_all.root")['tautau_genHH_mass']
    genHH_mass_ditaujet = uproot.open("data/regions_2023_postBPix_15p36-ditaujet_all.root")['tautau_genHH_mass']
    genHH_mass_quadjet = uproot.open("data/regions_2023_postBPix_15p34_all.root")['tautau_genHH_mass']
    genHH_mass_both = uproot.open("data/regions_2023_postBPix_15p33_all.root")['tautau_genHH_mass']
    # ht = uproot.open("data/regions_postBPix_15p11-quadjet-confirm_all.root")['tautau_SoftActivityJetHT']
    petroff6 = ["#5790fc", "#f89c20", "#e42536", "#964a8b", "#9c9ca1", "#7a21dd"]

    ditaujet_gain = genHH_mass_ditaujet.values() - genHH_mass_base_tau.values()
    quadjet_gain = genHH_mass_quadjet.values() - genHH_mass_base_tau.values()
    overlap_gain = genHH_mass_ditaujet.values() + genHH_mass_quadjet.values() - genHH_mass_base_tau.values() - genHH_mass_both.values()

    hep.style.use("CMS")
    
    ditaujet_gain_err = np.sqrt(ditaujet_gain)/genHH_mass_base_tau.values()
    quadjet_gain_err = np.sqrt(quadjet_gain)/genHH_mass_base_tau.values()
    overlap_gain_err = np.sqrt(overlap_gain)/genHH_mass_base_tau.values()

    ax = genHH_mass_base_tau.axes[0].centers()
    edges = genHH_mass_base_tau.axes[0].edges()
    xerr_r = ax - edges[:-1]
    xerr_l = edges[1:] - ax

    plt.figure(1)
    plt.bar(x=ax, height=(ditaujet_gain/genHH_mass_base_tau.values()) * 100,width=xerr_r * 2, yerr=ditaujet_gain_err * 100, color="#fff", edgecolor=petroff6[0], ecolor=petroff6[0], label="DiTauJet & !DiTau")
    
    plt.bar(x=ax, height=-(quadjet_gain/genHH_mass_base_tau.values()) * 100, width=xerr_r * 2, yerr=quadjet_gain_err * 100, label="QuadJet & !DiTauJet", color="#fff", edgecolor=petroff6[1], ecolor=petroff6[1])
    plt.bar(x=ax, height=(overlap_gain/genHH_mass_base_tau.values()) * 100, width=xerr_r * 2, yerr=overlap_gain_err * 100, color="#fff", edgecolor=petroff6[2], ecolor=petroff6[2], label="DiTauJet & QuadJet")
    plt.bar(x=ax, height=-(overlap_gain/genHH_mass_base_tau.values()) * 100, width=xerr_r * 2, yerr=overlap_gain_err * 100, color="#fff", edgecolor=petroff6[2], ecolor=petroff6[2])
    plt.axhline(y=0, color='black', linestyle='-')
    ticks_locs, ticks_labels = plt.yticks()
    plt.yticks(ticks_locs, [int(abs(tick)) for tick in ticks_locs])

    plt.xlabel(r"m$_{HH}$ [GeV]")
    plt.ylabel("Gain [%]")
    plt.legend()
    hep.cms.label(lumi="9.96", com="13.6")
    plt.savefig("orthogonal_gains_postBPix.png".format(args.channels[0]))

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
    parser.add_argument('--year', required=True, type=str, choices=('2016', '2017', '2018', '2022', '2023'),
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