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
    #return os.path.join(path, "eff_Data_Mu_MC_DY_TT_WJets_" + channel + "_" + variable + "_TRG_METNoMu120_CUTS_*.root")
    return os.path.join(path, "eff_Data_Mu_MC_TT_DY_" + channel + "_" + variable + "_TRG_METNoMu120_CUTS_*.root")

def sigmoid(x, params):
    """
    Sigmoid function to mimick the TF1 object.
    Uproot does not yet support TF1 reading.
    """
    return params[2] / (1 + np.exp(-params[0] * (x - params[1])))

def compare_ratios(tag, channels, variable, var_units):
    base = os.path.join("/data_CMS/cms/alves/TriggerScaleFactors/", tag, "Outputs")

    bpath1 = build_path(base, channels[0], variable)
    path1 = glob.glob(bpath1)
    if len(path1) != 1:
        print(bpath1)
        print(path1)
        raise RuntimeError('s[ERROR] Path must have lenght 1.')
    bpath2 = build_path(base, channels[1], variable)
    path2 = glob.glob(bpath2)
    if len(path2) != 1:
        print(bpath2)
        print(path2)
        raise RuntimeError('[ERROR] Path must have lenght 1.')
    assert len(path2) == 1

    mu = '\u03BC'
    tau = '\u03C4'
    dd = {"mumu": mu+mu,
          "mutau": mu+tau}
    
    g1 = up.open(path1[0] + ":SF1D")
    g2 = up.open(path2[0] + ":SF1D")
    fit1_data = up.open(path1[0] + ":SigmoidFuncData")
    fit1_mc   = up.open(path1[0] + ":SigmoidFuncMC")
    fit2_data = up.open(path2[0] + ":SigmoidFuncData")
    fit2_mc   = up.open(path2[0] + ":SigmoidFuncMC")

    fit1_data_pars   = fit1_data.member('fFormula').member('fClingParameters')[:]
    fit1_data_xrange = (fit1_data.member('fXmin'), fit1_data.member('fXmax'))
    fit1_mc_pars     = fit1_mc.member('fFormula').member('fClingParameters')[:]
    fit1_mc_xrange   = (fit1_mc.member('fXmin'), fit1_mc.member('fXmax'))
    assert fit1_data_xrange == fit1_mc_xrange

    fit2_data_pars   = fit2_data.member('fFormula').member('fClingParameters')[:]
    fit2_data_xrange = (fit2_data.member('fXmin'), fit2_data.member('fXmax'))
    fit2_mc_pars     = fit2_mc.member('fFormula').member('fClingParameters')[:]
    fit2_mc_xrange   = (fit2_mc.member('fXmin'), fit2_mc.member('fXmax'))
    assert fit2_data_xrange == fit2_mc_xrange

    fit1_data_xvals = np.linspace(*fit1_data_xrange, num=500)
    fit1_data_yvals = sigmoid(fit1_data_xvals, fit1_data_pars)
    fit1_mc_xvals   = np.linspace(*fit1_mc_xrange, num=500)
    fit1_mc_yvals   = sigmoid(fit1_mc_xvals, fit1_mc_pars)

    fit2_data_xvals = np.linspace(*fit2_data_xrange, num=500)
    fit2_data_yvals = sigmoid(fit2_data_xvals, fit2_data_pars)
    fit2_mc_xvals   = np.linspace(*fit2_mc_xrange, num=500)
    fit2_mc_yvals   = sigmoid(fit2_mc_xvals, fit2_mc_pars)

    ratio1 = fit1_data_yvals / fit1_mc_yvals
    ratio2 = fit2_data_yvals / fit2_mc_yvals
    met_cut = 180.
    
    fig, (ax1, ax2) = plt.subplots(2, sharex=True, gridspec_kw={'height_ratios': [3., 1.]})
    plt.subplots_adjust(wspace=0, hspace=0)
    ax1.plot(fit1_data_xvals, ratio1, '--', color="blue")
    ax1.errorbar(g1.values(axis="x"), g1.values(axis="y"),
                 xerr=(g1.errors(axis="x", which="low"),g1.errors(axis="x", which="high")),
                 yerr=(g1.errors(axis="y", which="low"),g1.errors(axis="y", which="high")),
                 fmt='o', label=dd[channels[0]], color="blue")

    ax1.plot(fit2_data_xvals, ratio2, '--', color="orange")
    ax1.errorbar(g2.values(axis="x"), g2.values(axis="y"),
                 xerr=(g2.errors(axis="x", which="low"),g2.errors(axis="x", which="high")),
                 yerr=(g2.errors(axis="y", which="low"),g2.errors(axis="y", which="high")),
                 fmt='o', label=dd[channels[1]], color="orange")

    ax1.axvline(x=met_cut, color="grey", linestyle="--")
    ax1.axhline(y=1., color="grey", linestyle="--")
    ax1.set_ylabel("Efficiency Ratio", fontsize=20)
    ax1.legend(loc="best")

    ax2.axvline(x=met_cut, color="grey", linestyle="--")
    ax2.axhline(y=1., color="grey", linestyle="--")
    ax2.set_ylim(0.905, 1.095)
    ax2.plot(fit1_data_xvals, ratio1 / ratio2, '--', color="green")
    ax2.set_ylabel("{} / {}".format(dd[channels[0]], dd[channels[1]]), fontsize=20)
    ax2.set_xlabel(variable + " [" + var_units + "]", fontsize=21)
    
    hep.cms.text(' Preliminary', fontsize=22, ax=ax1)
    hep.cms.lumitext(r"59.7 $fb^{-1}$ (13 TeV)", fontsize=21, ax=ax1)

    output = __file__[:-3] + '.png'
    fig.savefig(output)
    print('Plot save under {}.'.format(output))

    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Compare efficiency ratios obtained with two different methods.')
    parser.add_argument('--tag', help='Tag used to produce the graphs. Same used by the inclusion/run.py command.')

    FLAGS = parser.parse_args()

    compare_ratios(FLAGS.tag,
                   channels=("mutau", "mumu"),
                   variable="metnomu_et", var_units="GeV")
