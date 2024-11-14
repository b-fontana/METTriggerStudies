# Coding: utf-8

_all_ = [ 'drawCuts' ]

# script to be applied on the output of the KLUB skimming
# stored in this repo due to the old CMSSW release we are forced to use in KLUB

import os
import argparse
import glob
import uproot as up
import hist
from hist.intervals import clopper_pearson_interval as clop

import matplotlib; import matplotlib.pyplot as plt
import matplotlib.colors as colors
import mplhep as hep
plt.style.use(hep.style.ROOT)

def compute_acceptance_times_efficiency(basepath, masses, spin, year):
    nums_leg, nums_met, nums_tau, denoms = {}, {}, {}, {}
    signal = {'0': 'Rad', '2': 'Grav'}
    for mx in masses:
        path = os.path.join(basepath, signal[spin] + str(mx), "output_*.root")
        flist = glob.glob(path)

        n1_leg, n2_leg, n3_leg, n4_leg = 0, 0, 0, 0
        n1_met, n2_met, n3_met, n4_met = 0, 0, 0, 0
        n1_tau, n2_tau, n3_tau, n4_tau = 0, 0, 0, 0
        d1, d2, d3, d4 = 0, 0, 0, 0
        for afile in flist:
            with up.open(afile) as upfile:
                histo = upfile['h_AccEff'].to_hist()
                n1_leg += histo['MuTau_legacy']
                n1_met += histo['MuTau_met']
                n1_tau += histo['MuTau_singletau']
                d1 += histo['MuTau']

                n2_leg += histo['ETau_legacy']
                n2_met += histo['ETau_met']
                n2_tau += histo['ETau_singletau']
                d2 += histo['ETau']
                
                n3_leg += histo['TauTau_legacy']
                n3_met += histo['TauTau_met']
                n3_tau += histo['TauTau_singletau']
                d3 += histo['TauTau']
                
                n4_leg += histo['MuMu_legacy']
                n4_met += histo['MuMu_met']
                n4_tau += histo['MuMu_singletau']
                d4 += histo['MuMu']

        nums_leg[mx]   = (n1_leg, n2_leg, n3_leg, n4_leg)
        nums_met[mx]   = (n1_leg + n1_met, n2_leg + n2_met, n3_leg + n3_met, n4_leg + n4_met)
        nums_tau[mx]   = (n1_leg + n1_met + n1_tau, n2_leg + n2_met + n2_tau,
                          n3_leg + n3_met + n3_tau, n4_leg + n4_met + n4_tau)
        denoms[mx] = (d1, d2, d3, d4)

    return nums_leg, nums_met, nums_tau, denoms

def plot_acceptance_times_efficiency(vals, spin, compare, lowm):
    masses = [float(k) for k in vals[0].keys()]
    chn_unicodes = {"ETau":   r'$e\tau$',
                    "MuTau":  r'$\mu\tau$',
                    "TauTau": r'$\tau\tau$',
                    "MuMu":   r'$\mu\mu$'}
    spin_map = {'0': "Radion", '2': "BulkGraviton"}
    cov = 0.68
    
    nums_leg, nums_met, nums_tau, denoms = vals

    if compare:
        channels = ("LepTau", "", "TauTau", "MuMu")
    else:
        channels = ("MuTau", "ETau", "TauTau", "MuMu")
    assert len(nums_leg[list(nums_leg.keys())[0]]) == len(channels) # all numerators
    assert len(nums_met[list(nums_met.keys())[0]]) == len(channels) # all numerators
    assert len(nums_tau[list(nums_tau.keys())[0]]) == len(channels) # all numerators
    assert len(denoms[list(denoms.keys())[0]]) == len(channels) # all denominators
    
    for ichn, chn in enumerate(channels):
        if chn == "":
            continue

        fig, ax = plt.subplots(1)
        plt.subplots_adjust(wspace=0, hspace=0)
        ax.set_xlabel("m(X) [GeV]", fontsize=20)
        ax.set_ylabel(r"Acceptance $\times$ Efficiency", fontsize=20)

        if chn == "LeptonTau":
            yvals_leg = [0 if d[0]==0 else (n[0]+n[1])/(d[0]+d[1]) for n,d in zip(nums_leg.values(),denoms.values())]
            yvals_met = [0 if d[0]==0 else (n[0]+n[1])/(d[0]+d[1]) for n,d in zip(nums_met.values(),denoms.values())]
            yvals_tau = [0 if d[0]==0 else (n[0]+n[1])/(d[0]+d[1]) for n,d in zip(nums_tau.values(),denoms.values())]

            err_leg = [clop(num=n[ichn],denom=d[ichn],coverage=cov) for n,d in zip(nums_leg.values(),denoms.values())]
            err_met = [clop(num=n[ichn],denom=d[ichn],coverage=cov) for n,d in zip(nums_met.values(),denoms.values())]
            err_tau = [clop(num=n[ichn],denom=d[ichn],coverage=cov) for n,d in zip(nums_tau.values(),denoms.values())]

        else:
            yvals_leg = [0 if d[ichn]==0 else n[ichn]/d[ichn] for n,d in zip(nums_leg.values(),denoms.values())]
            yvals_met = [0 if d[ichn]==0 else n[ichn]/d[ichn] for n,d in zip(nums_met.values(),denoms.values())]
            yvals_tau = [0 if d[ichn]==0 else n[ichn]/d[ichn] for n,d in zip(nums_tau.values(),denoms.values())]

            err_leg = [clop(num=n[ichn],denom=d[ichn],coverage=cov) for n,d in zip(nums_leg.values(),denoms.values())]
            err_met = [clop(num=n[ichn],denom=d[ichn],coverage=cov) for n,d in zip(nums_met.values(),denoms.values())]
            err_tau = [clop(num=n[ichn],denom=d[ichn],coverage=cov) for n,d in zip(nums_tau.values(),denoms.values())]
        
        edo_leg, eup_leg = [x for x,_ in err_leg], [x for _,x in err_leg]
        edo_rel_leg, eup_rel_leg = [y-x for y,x in zip(yvals_leg,edo_leg)], [x-y for y,x in zip(yvals_leg,eup_leg)]

        edo_met, eup_met = [x for x,_ in err_met], [x for _,x in err_met]
        edo_rel_met, eup_rel_met = [y-x for y,x in zip(yvals_met,edo_met)], [x-y for y,x in zip(yvals_met,eup_met)]

        edo_tau, eup_tau = [x for x,_ in err_tau], [x for _,x in err_tau]
        edo_rel_tau, eup_rel_tau = [y-x for y,x in zip(yvals_tau,edo_tau)], [x-y for y,x in zip(yvals_tau,eup_tau)]

        gap = 2 if lowm else 6
        ax.errorbar(masses, yvals_leg,
                    yerr=(edo_rel_leg,eup_rel_leg), fmt='--o', color="red", label="Legacy")
        ax.errorbar([m+gap for m in masses], yvals_met,
                    yerr=(edo_rel_met,eup_rel_met), fmt='--^', color="orange", label="Legacy + MET")
        ax.errorbar([m+gap*2 for m in masses], yvals_tau,
                    yerr=(edo_rel_tau,eup_rel_tau), fmt='--s', color="green", label=r"Legacy + MET + Single-$\tau$")
        
        # ax.plot(masses, yvals_leg, '--o', color="red", label="Legacy")
        # ax.plot([m+10 for m in masses], yvals_met, '--^', color="orange", label="Legacy + MET")
        # ax.plot([m+20 for m in masses], yvals_tau, '--s', color="green", label=r"Legacy + MET + Single-$\tau$")
        
        if compare:
            ax.legend(loc="lower right")
            ax.set_xlim(200, 1700)
            # ax.set_ylim(0., 0.149)
        elif lowm:
            ax.legend(loc="upper left")
        else:
            ax.legend(loc="upper right")
        
        hep.cms.text(' Preliminary', fontsize=22, ax=ax)
        if compare and chn == "LepTau":
            lumitext = spin_map[spin] + ', ' + chn_unicodes["ETau"] + ' + ' + chn_unicodes["MuTau"]
        else:
            lumitext = spin_map[spin] + ', ' + chn_unicodes[chn]
        hep.cms.lumitext(lumitext  + ' | ' + r"138 $fb^{-1}$ (13 TeV)", fontsize=21, ax=ax)
        
        output = "AccEff_" + chn
        if compare:
            output += "_comp"
        if lowm:
            output += "_lowm"
        for ext in ('.png', '.pdf'):
            fig.savefig(output + ext)
        plt.close()
    
if __name__ == '__main__':
    com = "python3 tests/acceptanceTimesEfficiency.py --year 2018"
    parser = argparse.ArgumentParser(description="Draw acceptance vs. efficiency plots. Example: '{}'".format(com),
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--basepath', default="/data_CMS/cms/alves/HHresonant_SKIMS/SKIMS_UL18_AccEff/",
                        type=str, help='base path where signal input files are stored.')
    parser.add_argument('--masses', nargs='+', type=str, help='Resonance masses',
                        default=(250, 260, 270, 280, 300, 320, 350, 400, 450, 500, 550, 600, 650,
                                 700, 750, 800, 850, 900, 1000, 1250, 1500, 1750, 2000, 2500, 3000))
    parser.add_argument('--year', required=True, choices=("2016", "2016APV", "2017", "2018"),
                        type=str, help='Data taking period.')
    parser.add_argument('--spin', required=True, choices=('0', '2'), help='Signal spin hypothesis.')
    parser.add_argument('--lowm', action='store_true',
                        help="Run on low masses only.")
    parser.add_argument('--compare_atlas', action='store_true',
                        help="Format the plot to faciliate a comparison with the ATLAS acceptance times efficiency non-resonant plot.")
    
    FLAGS = parser.parse_args()

    if FLAGS.compare_atlas:
        FLAGS.masses = (250, 260, 270, 280, 300, 320, 350, 400, 450, 500, 550, 600, 650,
                        700, 750, 800, 850, 900, 1000, 1250, 1500)
    if FLAGS.lowm:
        FLAGS.masses = (250, 260, 270, 280, 300, 320, 350, 400, 450, 500, 550, 600, 650)
        
    values = compute_acceptance_times_efficiency(FLAGS.basepath, FLAGS.masses, FLAGS.spin, FLAGS.year)
    plot_acceptance_times_efficiency(values, FLAGS.spin, FLAGS.compare_atlas, FLAGS.lowm)
