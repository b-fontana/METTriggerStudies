# Coding: utf-8

_all_ = [ 'drawCuts' ]

import os
import argparse
import glob
import uproot as up
import hist

import matplotlib; import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.style.ROOT)

def sel_category(batch):
    """Applies analysis category cuts to an awkward batch."""
    pass

def sel_cuts(batch, channel):
    """Applies analysis selection cuts to an awkward batch."""
    chn_map = {"etau": 1, "mutau": 0, "tautau": 2, "mumu": 3}
    batch = batch[batch.pairType == chn_map[channel]]
    
    # When one only has 0 or 1 bjet th HH mass is not well defined,
    # and a value of -1 is assigned. One thus has to remove the cut below
    # when considering events with less than 2 b-jets.
    batch = batch[batch.HHKin_mass > 1]

    # third lepton veto
    batch = batch[batch.nleps == 0]

    # require at least two b jet candidates
    batch = batch[batch.nbjetscand > 1]

    # opposite sign leptons
    batch = batch[batch.isOS == 1]

    # Loose / Medium / Tight
    iso_allowed = { 'dau1_ele': 1., 'dau1_mu': 0.15, 'dau2_mu': 0.15,
                    'dau1_tau': 5., 'dau2_tau': 5. }

    # lepton id and isolation
    if channel == "etau":
        batch = batch[(batch.dau1_eleMVAiso == 1.) & (batch.dau2_deepTauVsJet >= 5.)]
    elif channel == "mutau":
        batch = batch[(batch.dau1_iso < 0.15) & (batch.dau2_deepTauVsJet >= 5.)]
    elif channel == "tautau":
        batch = batch[(batch.dau1_deepTauVsJet >= 5.) & (batch.dau2_deepTauVsJet >= 5.)]
    elif channel = "mumu":
        batch = batch[(batch.dau1_iso < 0.15) & (batch.dau2_iso < 0.15)]

    return batch
    
def createHisto(x, y, inputs, channel, other_vars=None):
    avars = (x, y) if other_vars is None else (x, y, *other_vars)
    for inp in inputs:
        assert inp[-5:] == ".root"

    xmax, ymax = 240, 240
    histogram = hist.Hist(
        hist.axis.Regular(50, 15, xmax, name=x),
        hist.axis.Regular(50, 15, ymax, name=y),
        storage=hist.storage.Double()
    )

    for ibatch, batch in enumerate(up.iterate(inputs, step_size="200MB", library='ak',
                                              filter_name=avars)):
        print("Batch {}".format(ibatch))
        batch = sel_cuts(batch, channel)
        batch = sel_category(batch)
        histogram.fill(batch.dau1_pt, batch.dau2_pt)

    return histogram

def drawCuts(inputs, dtype, sample, channel):
    xvar, yvar = "dau1_pt", "dau2_pt"
    other_vars = ('HHKin_mass', 'pairType', 'isOS', 'dau1_eleMVAiso', 'dau1_iso', 'dau2_iso',
                  'dau1_deepTauVsJet', 'dau2_deepTauVsJet', 'nleps', 'nbjetscand')

    histogram = createHisto(x=xvar, y=yvar, channel=channel,
                            other_vars=other_vars, inputs=inputs)

    wsize, hsize = 16, 16
    
    fig = plt.figure(figsize=(wsize, hsize),)
    ax = plt.subplot(111)
    ax.title.set_size(100)

    ax.set_xlabel(r"$p^1_T$ [GeV]")
    ax.set_ylabel(r"$p^2_T$ [GeV]")
    
    cbar = hep.hist2dplot(histogram.values(), histogram.axes[0].edges, histogram.axes[1].edges,
                          flow=None)
    cbar.cbar.ax.set_ylabel("No. Events", rotation=90, labelpad=0.5, loc='top')

    hep.cms.text('Preliminary', fontsize=wsize*2.5)
    hep.cms.lumitext(sample, fontsize=wsize*1.5) # r"138 $fb^{-1}$ (13 TeV)"

    # lines
    xmin, xmax = histogram.axes[0].edges[0], histogram.axes[0].edges[-1]
    ymin, ymax = histogram.axes[1].edges[0], histogram.axes[1].edges[-1]
    plt.plot([42., 42., xmax], [ymax, 42., 42.], c='lightgreen', linewidth=10, label="DiTau")
    plt.plot([xmin, 39., 39.], [191., 191., ymax], c='red', linewidth=10, label="SingleTau")
    plt.plot([xmax, 191., 191.], [39., 39., ymin], c='red', linewidth=10)
    plt.plot([189., 189., 39., 39., xmin], [ymin, 39., 39., 189., 189.], c='blue', linewidth=10, label="MET")

    rect1 = matplotlib.patches.Rectangle((188,211), 49, 26, color='white')
    ax.add_patch(rect1)
    
    plt.legend(loc="upper right", facecolor="white", edgecolor="white", framealpha=1)
    
    for ext in ('.pdf',):
        savename = "plot"
        plt.savefig(savename + ext, dpi=600)
        print('Stored in {}'.format(savename + ext))
    plt.close()
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="draw cuts on top of signal or MC distributions",
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('--dtype', choices=('signal', 'mc'), default='signal',
                        type=str, help='Data type')
    parser.add_argument('--signal', choices=('Radion', 'BulkGraviton'), default='Radion',
                        type=str, help='Signal particle type')
    parser.add_argument('--mass', default='1000',
                        type=str, help='Signal particle type')
    parser.add_argument('--channel', required=True, choices=("etau", "mutau", "tautau", "mumu"),
                        type=str, help='Signal particle type')
    FLAGS = parser.parse_args()

    path = {'signal': "/data_CMS/cms/alves/HHresonant_SKIMS/SKIMS_UL18_validateMETNoSF_Sig_V2/",
            'mc'    : "/data_CMS/cms/alves/HHresonant_SKIMS/SKIMS_UL18_validateMETNoSF_MC_V2/"}
    name = "GluGluTo{}ToHHTo2B2Tau_M-{}_".format(FLAGS.signal, FLAGS.mass)
    
    infiles = glob.glob(os.path.join(path[FLAGS.dtype], name, "output_*.root"))

    drawCuts(infiles, dtype=FLAGS.dtype, channel=FLAGS.channel,
             sample=FLAGS.signal + " " + FLAGS.mass + " GeV")
