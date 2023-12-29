# Coding: utf-8

_all_ = [ 'compare_ratios' ]

import os
import glob
import argparse
import uproot as up

import matplotlib
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.style.ROOT)

def build_path(base, channel, variable):
    path = os.path.join(base, channel, variable)
    return os.path.join(path, "eff_Data_Mu_MC_TT_DY_" + channel + "_" + variable + "_TRG_METNoMu120_CUTS_*.root")

def compare_ratios(tag, channels, variable, var_units):
    base = os.path.join("/data_CMS/cms/alves/TriggerScaleFactors/", tag, "Outputs")

    path1 = glob.glob(build_path(base, channels[0], variable))
    assert len(path1) == 1
    path2 = glob.glob(build_path(base, channels[1], variable))
    assert len(path2) == 1
    
    g1 = up.open(path1[0] + ":SF1D")
    g2 = up.open(path2[0] + ":SF1D")
    # fit1 = up.open(path1[0] + ":SigmoidFuncSF")
    # fit2 = up.open(path2[0] + ":SigmoidFuncSF")
    # breakpoint()
    plt.figure()
    plt.errorbar(g1.values(axis="x"), g1.values(axis="y"),
                 xerr=(g1.errors(axis="x", which="low"),g1.errors(axis="x", which="high")),
                 yerr=(g1.errors(axis="y", which="low"),g1.errors(axis="y", which="high")),
                 fmt='o', label=channels[0])
    plt.errorbar(g2.values(axis="x"), g2.values(axis="y"),
                 xerr=(g2.errors(axis="x", which="low"),g2.errors(axis="x", which="high")),
                 yerr=(g2.errors(axis="y", which="low"),g2.errors(axis="y", which="high")),
                 fmt='o', label=channels[1])
    plt.axvline(x=180., color="black", linestyle="--")
    plt.ylabel("Efficiency Ratio", fontsize=20)
    plt.xlabel(variable + " [" + var_units + "]", fontsize=21)
    plt.legend(loc="best")

    hep.cms.text(' Preliminary', fontsize=22)
    hep.cms.lumitext(r"59.7 $fb^{-1}$ (13 TeV)", fontsize=21)

    output = __file__[:-3] + '.png'
    plt.savefig(output)
    print('Plot save under {}.'.format(output))

    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Compare efficiency ratios obtained with two different methods.')
    parser.add_argument('--tag', help='Tag used to produce the graphs. Same used by the inclusion/run.py command.')

    FLAGS = parser.parse_args()

    compare_ratios(FLAGS.tag,
                   channels=("mutau", "mumu"),
                   variable="metnomu_et", var_units="GeV")
