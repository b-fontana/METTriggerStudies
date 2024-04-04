# Coding: utf-8

_all_ = [ 'drawCuts' ]

# script to be applied on the output of the KLUB skimming
# stored in this repo due to the old CMSSW release we are forced to use in KLUB

import os
import argparse
import glob
import uproot as up
import hist

import matplotlib; import matplotlib.pyplot as plt
import matplotlib.colors as colors
import mplhep as hep
plt.style.use(hep.style.ROOT)

def compute_acceptance_times_efficiency(basepath, masses, spin, year):
    eff = {}
    signal = {'0': 'GluGluToRadionToHHTo2B2Tau',
              '2': 'GluGluToGravitonToHHTo2B2Tau'}
    for mx in masses:
        path = os.path.join(basepath, "{}_M-{}_".format(signal[spin], mx), "output_*.root")
        flist = glob.glob(path)

        num, den = 0, 0
        for afile in flist:
            with up.open(afile) as upfile:
                histo = upfile['h_effSummary'].to_hist()
                num += histo['Trigger']
                den += histo['all']
        eff[mx] = 0. if den == 0 else num / den

    return eff

def plot_acceptance_times_efficiency(vals, spin):
    colors = iter(("red", "dodgerblue"))
    
    fig, ax = plt.subplots(1)
    plt.subplots_adjust(wspace=0, hspace=0)
    ax.set_xlabel("mX", fontsize=20)
    ax.set_ylabel("Acceptance x efficiency", fontsize=20)

    masses = [float(k) for k in vals.keys()]
    yvals = [k for k in vals.values()]
    err = [x/6 for x in yvals]
    ax.errorbar(masses, yvals,
                yerr=(err, err), fmt='-o', color=next(colors), label='test')

    ax.legend(loc="upper right")

    chn_unicodes = {"etau":   r'$bb\: e\tau$',
                    "mutau":  r'$bb\: \mu\tau$',
                    "tautau": r'$bb\: \tau\tau$',
                    "mumu":   r'$bb\: \mu\mu$'}
    spin_map = {'0': "Radion", '2': "BulkGraviton"}
     
    hep.cms.text(' Preliminary', fontsize=22, ax=ax)
    hep.cms.lumitext(chn_unicodes["mutau"] + " (baseline) | " + spin_map[spin],
                     fontsize=21, ax=ax)

    output = "acceptanceTimesEfficiency"
    for ext in ('.png', '.pdf'):
        fig.savefig(output + ext)
    
if __name__ == '__main__':
    com = "python3 tests/acceptanceTimesEfficiency.py --year 2018"
    parser = argparse.ArgumentParser(description="Draw acceptance vs. efficiency plots. Example: '{}'".format(com),
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--basepath', default="/data_CMS/cms/alves/HHresonant_SKIMS/SKIMS_UL18_OpenCADI_Sig/",
                        type=str, help='base path where signal input files are stored.')
    parser.add_argument('--masses', nargs='+', type=str,
                        default=(260,300,500,1000), help='Resonance masses')
    parser.add_argument('--year', required=True, choices=("2016", "2016APV", "2017", "2018"),
                        type=str, help='Data taking period.')
    parser.add_argument('--spin', required=True, choices=('0', '2'), help='Signal spin hypothesis.')
    
    FLAGS = parser.parse_args()

    values = compute_acceptance_times_efficiency(FLAGS.basepath, FLAGS.masses, FLAGS.spin, FLAGS.year)
    plot_acceptance_times_efficiency(values, FLAGS.spin)
