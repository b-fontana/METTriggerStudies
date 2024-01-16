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

def sel_category(batch, category, year):
    """Applies analysis category cuts to an awkward batch."""
    deepJetWP = {'2016'    : (0.048, 0.249),
                 '2016APV' : (0.051, 0.260),
                 '2017'    : (0.0532, 0.3040),
                 '2018'    : (0.0490, 0.2783)}[year]

    if category == "sboosted":
        batch = batch[(batch.isBoosted == 1) &
                      (batch.bjet1_bID_deepFlavor > deepJetWP[0]) & (batch.bjet2_bID_deepFlavor > deepJetWP[0])]
    elif category == "s1b1jresolved":
        batch = batch[(batch.isBoosted != 1) &
                      (((batch.bjet1_bID_deepFlavor > deepJetWP[1]) & (batch.bjet2_bID_deepFlavor < deepJetWP[1])) |
                       ((batch.bjet1_bID_deepFlavor < deepJetWP[1]) & (batch.bjet2_bID_deepFlavor > deepJetWP[1])))]
    elif category == "s2b0jresolved":
        batch = batch[(batch.isBoosted != 1) &
                      (batch.bjet1_bID_deepFlavor > deepJetWP[1]) & (batch.bjet2_bID_deepFlavor > deepJetWP[1])]
        
    return batch

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
    elif channel == "mumu":
        batch = batch[(batch.dau1_iso < 0.15) & (batch.dau2_iso < 0.15)]

    return batch
    
def createHisto(x, y, inputs, channel, category, year, other_vars=None):
    avars = (x, y) if other_vars is None else (x, y, *other_vars)
    for inp in inputs:
        assert inp[-5:] == ".root"

    if channel == "etau":
        xmax, ymax = 135, 135
        bins = (40, 15, xmax)
    elif channel == "mutau":
        xmax, ymax = 135, 135
        bins = (40, 15, xmax)
    if channel == "tautau":
        xmax, ymax = 240, 240
        bins = (40, 15, xmax)
    
    histogram = hist.Hist(
        hist.axis.Regular(*bins, name=x),
        hist.axis.Regular(*bins, name=y),
        storage=hist.storage.Double()
    )

    for ibatch, batch in enumerate(up.iterate(inputs, step_size="200MB", library='ak',
                                              filter_name=avars)):
        print("Batch {}".format(ibatch))
        batch = sel_cuts(batch, channel)
        batch = sel_category(batch, category, year)
        histogram.fill(batch.dau1_pt, batch.dau2_pt)

    return histogram

def drawCuts(inputs, dtype, sample, channel, category, year):
    xvar, yvar = "dau1_pt", "dau2_pt"
    other_vars = ('HHKin_mass', 'pairType', 'isOS', 'dau1_eleMVAiso', 'dau1_iso', 'dau2_iso',
                  'dau1_deepTauVsJet', 'dau2_deepTauVsJet', 'nleps', 'nbjetscand',
                  'bjet1_bID_deepFlavor', 'bjet2_bID_deepFlavor', 'isBoosted')

    histogram = createHisto(x=xvar, y=yvar, channel=channel, category=category, year=year,
                            other_vars=other_vars, inputs=inputs)

    wsize, hsize = 16, 16
    
    fig = plt.figure(figsize=(wsize, hsize),)
    ax = plt.subplot(111)
    ax.title.set_size(100)

    ax.set_xlabel(r"$p_T(1)$ [GeV]")
    ax.set_ylabel(r"$p_T(2)$ [GeV]")
    
    cbar = hep.hist2dplot(histogram.values(), histogram.axes[0].edges, histogram.axes[1].edges,
                          flow=None)
    cbar.cbar.ax.set_ylabel("No. Events", rotation=90, labelpad=0.5, loc='top')

    hep.cms.text('Preliminary', fontsize=wsize*2.5)

    chn_unicodes = {"etau":   r'$e\tau$',
                    "mutau":  r'$\mu\tau$',
                    "tautau": r'$\tau\tau$'}
    cat_map = {'baseline': "baseline", 'sboosted': "boosted",
               's1b1jresolved': "res 1b", 's2b0jresolved': "res 2b"}
    hep.cms.lumitext(chn_unicodes[channel] + " (" + cat_map[category] + ") | " + sample + " (" + year + ")",
                     fontsize=wsize*1.5) # r"138 $fb^{-1}$ (13 TeV)"

    # lines
    xmin, xmax = histogram.axes[0].edges[0], histogram.axes[0].edges[-1]
    ymin, ymax = histogram.axes[1].edges[0], histogram.axes[1].edges[-1]

    if channel == "etau":
        if year == "2016":
            plt.plot([26., 26.], [ymin, ymax],
                     c='lightgreen', linewidth=10, label=r"single-e + e$\tau$")
            plt.plot([24., 24.], [ymin, ymax],
                     c='blue', linewidth=10, label=r"MET")
        elif year == "2017" or year == "2018":
            plt.plot([26., 26., 34., 34.], [ymax, 26., 26., ymin],
                     c='lightgreen', linewidth=10, label=r"single-e + e$\tau$")
            plt.plot([24., 24., 32., 32.], [ymax, 24., 24., ymin],
                     c='blue', linewidth=10, label=r"MET")
        
    elif channel == "mutau":
        if year == "2016":
            plt.plot([21., 21., 26., 26.], [ymax, 26., 26., ymin],
                     c='lightgreen', linewidth=10, label=r"single-$\mu$ + $\mu\tau$")
            plt.plot([19., 19., 24., 24.], [ymax, 24., 24., ymin],
                     c='blue', linewidth=10, label=r"MET")
        elif year == "2017":
            plt.plot([22., 22., 29., 29.], [ymax, 33., 33., ymin],
                     c='lightgreen', linewidth=10, label=r"single-$\mu$ + $\mu\tau$")
            plt.plot([20., 20., 27., 27.], [ymax, 31., 31., ymin],
                     c='blue', linewidth=10, label=r"MET")
        elif year == "2018":
            plt.plot([22., 22., 26., 26.], [ymax, 33., 33., ymin],
                     c='lightgreen', linewidth=10, label=r"single-$\mu$ + $\mu\tau$")
            plt.plot([20., 20., 24., 24.], [ymax, 31., 31., ymin],
                     c='blue', linewidth=10, label=r"MET")

    elif channel == "tautau":
        plt.plot([42., 42., xmax], [ymax, 42., 42.], c='lightgreen', linewidth=10, label=r"$\tau\tau$")
        plt.plot([xmin, 39., 39.], [191., 191., ymax], c='red', linewidth=10, label=r"single-$\tau$")
        plt.plot([xmax, 191., 191.], [39., 39., ymin], c='red', linewidth=10)
        plt.plot([189., 189., 39., 39., xmin], [ymin, 39., 39., 189., 189.], c='blue', linewidth=10, label="MET")

    if channel == "etau":
        rect = matplotlib.patches.Rectangle((103,119), 30, 14, color='white')
    elif channel == "mutau":
        rect = matplotlib.patches.Rectangle((103,119), 30, 14, color='white')
    elif channel == "tautau":
        rect = matplotlib.patches.Rectangle((194,203), 42, 34, color='white')
    ax.add_patch(rect)

    plt.legend(loc="upper right", facecolor="white", edgecolor="white", framealpha=1, title="Triggers")
    
    for ext in ('.pdf',):
        savename = '_'.join(("drawCuts", sample.replace(' ', '-'), channel, category, year))
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
    parser.add_argument('--category', default="baseline",
                        choices=("baseline", "sboosted", "s1b1jresolved", "s2b0jresolved"),
                        type=str, help='Analysis category')
    parser.add_argument('--year', required=True, choices=("2016", "2016APV", "2017", "2018"),
                        type=str, help='Signal particle type')
    FLAGS = parser.parse_args()

    path = {'signal': "/data_CMS/cms/alves/HHresonant_SKIMS/SKIMS_UL18_validateMETNoSF_Sig_V2/",
            'mc'    : "/data_CMS/cms/alves/HHresonant_SKIMS/SKIMS_UL18_validateMETNoSF_MC_V2/"}
    name = "GluGluTo{}ToHHTo2B2Tau_M-{}_".format(FLAGS.signal, FLAGS.mass)
    
    infiles = glob.glob(os.path.join(path[FLAGS.dtype], name, "output_*.root"))

    drawCuts(infiles, dtype=FLAGS.dtype, channel=FLAGS.channel, year=FLAGS.year, category=FLAGS.category,
             sample=FLAGS.signal + " " + FLAGS.mass + " GeV")
