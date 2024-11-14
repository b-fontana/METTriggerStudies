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

tau = '\u03C4'
mu  = '\u03BC'
pm  = '\u00B1'
ditau = tau+tau

def main(args):
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

    err_ditaujet = genHH_mass_ditaujet.errors()/ genHH_mass_base_tau.values()
    err_quadjet = genHH_mass_quadjet.errors()/ genHH_mass_base_tau.values()
    err_both = genHH_mass_both.errors()/ genHH_mass_base_tau.values()

    ax = genHH_mass_baseline.axes[0].centers()    
    edges = genHH_mass_baseline.axes[0].edges()
    xerr_r = ax - edges[:-1]
    xerr_l = edges[1:] - ax
    hep.style.use("CMS")
    # print("CMS style")
    fig, ax1 = plt.subplots()
    # ax1.errorbar(ax[0:10], ht.values()[0:10], yerr=ht.errors()[0:10], xerr=[xerr_r[0:10], xerr_l[0:10]], fmt='o', label="SoftActivityJetHT", color=petroff6[2], zorder=5)
    # hep.cms.label(lumi="9.96", com="13.6")
    # plt.legend()
    # ax1.set_xlabel(r"HT [GeV]")
    # ax1.set_ylabel("N [a. u.]")
    # plt.savefig("ht_postBPix_pf75.png".format(args.channels[0]))
    
    # return 
    ax1.set_ylim(-4, 30)
    # print("subplots")
    ax2 = ax1.twinx()
    print("generated additional axis")
    ax2.set_ylabel("Weighted MC events [a.u.]")  # we already handled the x-label with ax1
    ax2.bar(ax, height=genHH_mass_baseline.values(), width=xerr_r * 2, color=petroff6[4], alpha=0.3, zorder=0, label="HHbbtautau")
    # print(channel_base_60.values(), channel_base_70.values(), channel_base_80.values())
    # hep.style.use("CMS")
    print(genHH_mass_ditaujet.values()/ genHH_mass_base_tau.values())
    ax1.errorbar(ax - 10, (genHH_mass_ditaujet.values() / genHH_mass_base_tau.values() - 1) * 100, yerr=err_ditaujet*100, xerr=[xerr_r - 10, xerr_l + 10], fmt='o', label="DiTau+Jet", color=petroff6[0], zorder=5)
    ax1.errorbar(ax, (genHH_mass_quadjet.values() / genHH_mass_base_tau.values() - 1) * 100, yerr=err_quadjet*100, xerr=[xerr_r, xerr_l], fmt='o', label="QuadJet parking", color=petroff6[2], zorder=5)
    ax1.errorbar(ax + 10, (genHH_mass_both.values() / genHH_mass_base_tau.values() - 1) * 100, yerr=err_both*100, xerr=[xerr_r + 10, xerr_l - 10], fmt='o', label="QuadJet + DiTau+Jet", color=petroff6[4], zorder=5)
    ax1.set_xlabel(r"m$_{HH}$ [GeV]")
    ax1.set_ylabel("Gain [%]")
    # print("here")
    lines, labels = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    plt.legend(lines + lines2, labels + labels2)
    # # print("Total Gain: {} +/- {}".format(round(channel_all.values().sum()/channel_base.values().sum() * 100, 2), round(err_tot_all * 100, 2)))
    hep.cms.label(lumi="9.7", com="13.6")
    plt.savefig("gains_ditaujet_postBPix.png".format(args.channels[0]))

    return

if __name__ == '__main__':
    desc = "Produce plots of trigger gain VS resonance mass.\n"
    desc += "Uses the output of setup_regions.py."
    desc += "When running on many channels, one should keep in mind each channel has different pT cuts."
    desc += "This might imply moving sub-folders (produced by the previous script) around."
    parser = argparse.ArgumentParser(description=desc, formatter_class=argparse.RawTextHelpFormatter)

    # parser.add_argument('--masses', required=True, nargs='+', type=str,
    #                     help='Resonance mass')
    parser.add_argument('--channels', required=True, nargs='+', type=str, 
                        choices=('etau', 'mutau', 'tautau'), default="tautau",
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
