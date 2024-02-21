# Coding: utf-8

_all_ = [ 'compare_ratios' ]

import os
import numpy as np
import glob
import argparse
import uproot as up

import matplotlib
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.style.ROOT)

def build_path(base, channel, variable):
    path = os.path.join(base, channel, variable)
    return os.path.join(path, "eff_Data_Mu_MC_TT_DY_WJets_" + channel + "_" + variable + "_TRG_METNoMu120_CUTS_*.root")


def sigmoid(x, params):
    """
    Sigmoid function to mimick the TF1 object.
    Uproot does not yet support TF1 reading.
    """
    return params[2] / (1 + np.exp(-params[0] * (x - params[1])))

def compare_ratios(tag, channels, variable, var_units):
    base = os.path.join("/data_CMS/cms/alves/TriggerScaleFactors/", tag, "Outputs")

    #bpaths = build_path(base, channels[0], variable)
    _paths = (build_path(base, channels[1], variable),
              "180_mumu_fit.root", "160_mumu_fit.root", "150_mumu_fit.root")

    paths = []
    for bpath in _paths:
        tmp = glob.glob(bpath)
        if len(tmp) != 1:
            print(tmp)
            raise RuntimeError('[ERROR] Path {} must have lenght 1.'.format(tmp))
        paths.append(tmp[0])

    mu = '\u03BC'
    tau = '\u03C4'
    dd = {"mumu": mu+mu,
          "mutau": mu+tau}

    colors = ("blue", "green", "red", "brown")
    labels = ("full", r"partial: $[180:\infty[$", r"partial: $[160:\infty[$", r"partial: $[150:\infty[$")
    met_cuts = [180., 160., 150.]
    var_map = dict(metnomu_et=r"MET (no-$\mu$)")

    ratios = []
    
    fig, (ax1, ax2) = plt.subplots(2, sharex=True, gridspec_kw={'height_ratios': [3., 1.]})
    plt.subplots_adjust(wspace=0, hspace=0)

    for ipath, bpath in enumerate(paths):
        graph = up.open(bpath + ":SF1D")
        
        fit_data = up.open(bpath + ":SigmoidFuncData")
        fit_mc = up.open(bpath + ":SigmoidFuncMC")

        fit_data_pars   = fit_data.member('fFormula').member('fClingParameters')[:]
        fit_data_xrange = (fit_data.member('fXmin'), fit_data.member('fXmax'))
        fit_mc_pars     = fit_mc.member('fFormula').member('fClingParameters')[:]
        fit_mc_xrange   = (fit_mc.member('fXmin'), fit_mc.member('fXmax'))
        assert fit_data_xrange == fit_mc_xrange

        # fit_data_xvals = np.linspace(*fit_data_xrange, num=500)
        # fit_mc_xvals   = np.linspace(*fit_mc_xrange, num=500)
        fit_xvals = np.linspace(0, 350, num=5000)

        fit_data_yvals = sigmoid(fit_xvals, fit_data_pars)
        fit_mc_yvals   = sigmoid(fit_xvals, fit_mc_pars)
        breakpoint()
        ratios.append(fit_data_yvals / fit_mc_yvals)
        
        ax1.plot(fit_data_xvals, ratios[-1], '--', color=colors[ipath])
        ax1.errorbar(graph.values(axis="x"), graph.values(axis="y"),
                     xerr=(graph.errors(axis="x", which="low"),graph.errors(axis="x", which="high")),
                     yerr=(graph.errors(axis="y", which="low"),graph.errors(axis="y", which="high")),
                     fmt='o', color=colors[ipath], label=labels[ipath])

        if ipath != 0:
            ax1.axvline(x=met_cuts[ipath-1], color="grey", linestyle="--")
            ax1.axhline(y=1., color="grey", linestyle="--")
            ax2.axvline(x=met_cuts[ipath-1], color="grey", linestyle="--")
            
    ax1.set_ylabel("Data / MC", fontsize=20)
    ax1.legend(loc="lower right")
    ax1.set_ylim(0.64, 1.015)
    for yval in (-0.05, 0., 0.05, 0.1):
        ax2.axhline(y=yval, color="grey", linestyle="--")
        ax2.set_ylim(-0.075, 0.145)

    # comparison of ratio relative using the first partial fit as reference
    for ipath, bpath in enumerate(paths[1:]):
        ax2.plot(fit_data_xvals, (ratios[1]/ratios[ipath+1])-1., '--', color=colors[ipath+1])

    ax2.set_ylabel(r"$(SF_{{180}}/SF) - 1$", fontsize=20)
    ax2.set_xlabel(var_map[variable] + " [" + var_units + "]", fontsize=21)
    
    hep.cms.text(' Preliminary', fontsize=22, ax=ax1)
    hep.cms.lumitext(r"59.7 $fb^{-1}$ (13 TeV)", fontsize=21, ax=ax1)

    output = os.path.join("/eos/home-b/bfontana/www/TriggerScaleFactors/CompareRatios/",
                          os.path.basename(__file__[:-3]))
    for ext in ('.png', '.pdf'):
        fig.savefig(output + ext)
        print('Plot saved under {}'.format(output))

    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Compare efficiency ratios obtained with two different methods.')
    parser.add_argument('--tag', help='Tag used to produce the graphs. Same used by the inclusion/run.py command.')

    FLAGS = parser.parse_args()

    compare_ratios(FLAGS.tag,
                   channels=("mutau", "mumu"),
                   variable="metnomu_et", var_units="GeV")
