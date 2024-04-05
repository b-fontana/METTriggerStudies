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
    vals = {}
    signal = {'0': 'GluGluToRadionToHHTo2B2Tau',
              '2': 'GluGluToGravitonToHHTo2B2Tau'}
    for mx in masses:
        path = os.path.join(basepath, "{}_M-{}_".format(signal[spin], mx), "output_*.root")
        flist = glob.glob(path)

        n1, n2, d = 0, 0, 0
        for afile in flist:
            with up.open(afile) as upfile:
                histo = upfile['h_effSummary'].to_hist()
                n1 += histo['METfilter']
                n2 += histo['Trigger']
                d += histo['all']
        vals[mx] = (n1, n2, d)

    return vals

def plot_acceptance_times_efficiency(vals, spin):
    colors = iter(("red", "dodgerblue"))
    masses = [float(k) for k in vals.keys()]
    chn_unicodes = {"etau":   r'$bb\: e\tau$',
                    "mutau":  r'$bb\: \mu\tau$',
                    "tautau": r'$bb\: \tau\tau$',
                    "mumu":   r'$bb\: \mu\mu$'}
    spin_map = {'0': "Radion", '2': "BulkGraviton"}
    
    fig, ax = plt.subplots(1)
    plt.subplots_adjust(wspace=0, hspace=0)
    ax.set_xlabel("m(X) [GeV]", fontsize=20)
    ax.set_ylabel(r"Acceptance $\times$ Efficiency", fontsize=20)

    channels = ("etau", "mutau")
    assert len(vals[list(vals.keys())[0]]) == len(channels) + 1 # all numerators + the denominator
    
    for ichn, chn in enumerate(channels):
        yvals = [0 if v[-1]==0 else v[ichn]/v[-1] for v in vals.values()]
        err = [clop(num=v[ichn],denom=v[-1],coverage=0.95) for v in vals.values()]
        edo, eup = [x for x,_ in err], [x for _,x in err]
        edo_rel = [y-x for y,x in zip(yvals,edo)]
        eup_rel = [x-y for y,x in zip(yvals,eup)]

        ax.errorbar(masses, yvals,
                    yerr=(edo_rel,eup_rel), fmt='-o', color=next(colors), label=chn_unicodes[chn])

    ax.legend(loc="upper right")

    hep.cms.text(' Simulation Preliminary', fontsize=22, ax=ax)
    hep.cms.lumitext(spin_map[spin] + ' | ' + r"138 $fb^{-1}$ (13 TeV)", fontsize=21, ax=ax)

    output = "acceptanceTimesEfficiency"
    for ext in ('.png', '.pdf'):
        fig.savefig(output + ext)
    
if __name__ == '__main__':
    com = "python3 tests/acceptanceTimesEfficiency.py --year 2018"
    parser = argparse.ArgumentParser(description="Draw acceptance vs. efficiency plots. Example: '{}'".format(com),
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--basepath', default="/data_CMS/cms/alves/HHresonant_SKIMS/SKIMS_UL18_OpenCADI_Sig/",
                        type=str, help='base path where signal input files are stored.')
    parser.add_argument('--masses', nargs='+', type=str, help='Resonance masses',
                        default=(250, 260, 270, 280, 300, 320, 350, 400, 450, 500, 550, 600, 650,
                                 700, 750, 800, 850, 900, 1000, 1250, 1500, 1750, 2000, 2500, 3000))
    parser.add_argument('--year', required=True, choices=("2016", "2016APV", "2017", "2018"),
                        type=str, help='Data taking period.')
    parser.add_argument('--spin', required=True, choices=('0', '2'), help='Signal spin hypothesis.')
    
    FLAGS = parser.parse_args()

    values = compute_acceptance_times_efficiency(FLAGS.basepath, FLAGS.masses, FLAGS.spin, FLAGS.year)
    plot_acceptance_times_efficiency(values, FLAGS.spin)
