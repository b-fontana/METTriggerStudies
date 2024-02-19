# Coding: utf-8

_all_ = [ 'drawCuts' ]

import os
import argparse
import glob
import uproot as up
import hist
import pickle

import matplotlib; import matplotlib.pyplot as plt
import matplotlib.colors as colors
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
    
    # When one only has 0 or 1 bjets the HH mass is not well defined,
    # and a value of -1 is assigned. One thus has to remove the cut below
    # when considering events with less than 2 b-jets.
    # batch = batch[batch.HHKin_mass > 1]

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
    
def getHisto(x, y, inputs, xbins, ybins,
             channel, category, year, dtype, savename, save=False, other_vars=None):
    avars = (x, y) if other_vars is None else (x, y, *other_vars)
    for inp in inputs:
        assert inp[-5:] == ".root"

    if save and os.path.isfile(savename):
        with open(savename, 'rb') as f:
            histogram = pickle.load(f)
        return histogram

    histogram = hist.Hist(
        hist.axis.Regular(*xbins, name=x),
        hist.axis.Regular(*ybins, name=y),
        storage=hist.storage.Double()
    )

    for ibatch, batch in enumerate(up.iterate(inputs, step_size="200MB", library='ak',
                                              filter_name=avars)):
        print("Batch {}".format(ibatch))
        batch = sel_cuts(batch, channel)
        batch = sel_category(batch, category, year)
        histogram.fill(getattr(batch, x), getattr(batch, y))

    with open(savename, 'wb') as f:
        pickle.dump(histogram, f)
        
    return histogram

class DrawCuts():
    def __init__(self, inputs, sample, channel, category, year, dtype, save):
        self.inputs = inputs
        self.sample = sample
        self.channel = channel
        self.category = category
        self.year = year
        self.dtype = dtype
        self.save = save

        self.other_vars = ('HHKin_mass', 'pairType', 'isOS', 'dau1_eleMVAiso', 'dau1_iso', 'dau2_iso',
                           'dau1_deepTauVsJet', 'dau2_deepTauVsJet', 'nleps', 'nbjetscand',
                           'bjet1_bID_deepFlavor', 'bjet2_bID_deepFlavor', 'isBoosted')

        self.outname = lambda mode: os.path.join('pickles',
                                                 '_'.join((mode, sample, channel, dtype, category, year)) + ".pkl")

        self.rect_opt = dict(facecolor='white', edgecolor="black", linewidth=2)
        self.leg_opt = dict(loc="upper right" if channel!="mumu" else "upper center",
                            facecolor="white", edgecolor="white", framealpha=1)

    def triggers(self):
        mode = "trigger"
        
        fig, ax = self.set_figure(16, 16)
        
        bins = self.get_bins(mode=mode, dtype=self.dtype)
        histogram = getHisto(x="dau1_pt", y="dau2_pt", xbins=bins[0], ybins=bins[1],
                             channel=self.channel, category=self.category, year=self.year,
                             other_vars=self.other_vars, inputs=self.inputs, dtype=self.dtype,
                             savename=self.outname(mode), save=self.save)

        xlabel, ylabel = self.set_axis_labels(mode=mode, channel=self.channel)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)

        values = histogram.values()
        cbar = hep.hist2dplot(values, histogram.axes[0].edges, histogram.axes[1].edges,
                              flow=None, norm=colors.LogNorm(vmin=1, vmax=values.max()))
        cbar.cbar.ax.set_ylabel("No. Events", rotation=90, labelpad=0.5, loc='top')
    
        self.set_hep_labels(self.channel)
        
        self.set_lines(histogram, mode=mode)

        if self.channel == "etau":
            rect = matplotlib.patches.Rectangle((76,202), 24, 35, **self.rect_opt)
        elif self.channel == "mutau":
            rect = matplotlib.patches.Rectangle((76,202), 24, 35, **self.rect_opt)
        elif self.channel == "tautau":
            rect = matplotlib.patches.Rectangle((194,203), 42, 34, **self.rect_opt)
        elif self.channel == "mumu":
            rect = matplotlib.patches.Rectangle((76,202), 24, 35, **self.rect_opt)

        ax.add_patch(rect)

        plt.legend(title="Triggers", **self.leg_opt)
        self.savefig(mode=mode)

    def mass(self):
        mode = "mass"
        
        fig, ax = self.set_figure(16, 16)
        
        bins = self.get_bins(mode=mode, dtype=self.dtype)
        histogram = getHisto(x="tauH_mass", y="bH_mass", xbins=bins[0], ybins=bins[1],
                             channel=self.channel, category=self.category, year=self.year,
                             other_vars=self.other_vars, inputs=self.inputs, dtype=self.dtype,
                             savename=self.outname(mode), save=self.save)

        xlabel, ylabel = self.set_axis_labels(mode=mode, channel=self.channel)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)

        values = histogram.values()
        cbar = hep.hist2dplot(values, histogram.axes[0].edges, histogram.axes[1].edges,
                              flow=None, norm=colors.LogNorm(vmin=1, vmax=values.max()))
        cbar.cbar.ax.set_ylabel("No. Events", rotation=90, labelpad=0.5, loc='top')
    
        self.set_hep_labels(self.channel)
     
        self.set_lines(histogram, mode=mode)

        rect = matplotlib.patches.Rectangle((140,332), 57, 13, **self.rect_opt)
        ax.add_patch(rect)

        plt.legend(loc="upper right", facecolor="white", edgecolor="white", framealpha=1)
        self.savefig(mode)
        
    def get_bins(self, mode, dtype):
        if mode == "trigger":
            if self.channel == "etau":
                nbinsx = 25 if dtype == "signal" else 35
                nbinsy = 25 if dtype == "signal" else 35
                xbins = (nbinsx, 5, 101)
                ybins = (nbinsy, 5, 240)
            elif self.channel == "mutau":
                nbinsx = 25 if dtype == "signal" else 35
                nbinsy = 25 if dtype == "signal" else 35
                xbins = (nbinsx, 10, 101)
                ybins = (nbinsy, 10, 240)
            elif self.channel == "tautau":
                nbinsx = 25 if dtype == "signal" else 40
                nbinsy = 25 if dtype == "signal" else 40
                xbins = (nbinsx, 12, 240)
                ybins = (nbinsy, 12, 240)
            elif self.channel == "mumu":
                nbinsx = 25 if dtype == "signal" else 35
                nbinsy = 25 if dtype == "signal" else 35
                xbins = (nbinsx, 10, 170)
                ybins = (nbinsy, 10, 130)

        elif mode == "mass":
            nbinsx = 50 if dtype == "signal" else 100
            nbinsy = 50 if dtype == "signal" else 100
            if self.channel == "etau":
                xbins = (nbinsx, 5, 200)
                ybins = (nbinsy, 5, 350)
            elif self.channel == "mutau":
                xbins = (nbinsx, 10, 200)
                ybins = (nbinsy, 10, 350)
            elif self.channel == "tautau":
                xbins = (nbinsx, 15, 200)
                ybins = (nbinsy, 15, 350)
            elif self.channel == "mumu":
                xbins = (nbinsx, 10, 200)
                ybins = (nbinsy, 10, 350)

        return xbins, ybins

    def savefig(self, mode):
        for ext in ('.png', '.pdf',):
            smpl = self.sample.replace(' ', '-').replace('+', '-')
            savename = os.path.join("/eos/home-b/bfontana/www/DrawCuts/", mode)
            savename = os.path.join(savename,
                                    '_'.join(("draw_" + mode, smpl,
                                              self.channel, self.category, self.year)))
            plt.savefig(savename + ext, dpi=600)
            print('Stored in {}'.format(savename + ext))
        plt.close('all')

    def set_figure(self, wsize, hsize):
        fig = plt.figure(figsize=(wsize, hsize),)
        ax = plt.subplot(111)
        ax.title.set_size(100)
        return fig, ax

    def set_axis_labels(self, mode, channel):
        """Set X and Y label."""
        if mode == "trigger":
            if channel == "etau":
                xlabel = r"$p_T(e)$ [GeV]"
                ylabel = r"$p_T(\tau)$ [GeV]"
            elif channel == "mutau":
                xlabel = r"$p_T(\mu)$ [GeV]"
                ylabel = r"$p_T(\tau)$ [GeV]"
            elif channel == "tautau":
                xlabel = r"$p_T(\tau_1)$ [GeV]"
                ylabel = r"$p_T(\tau_2)$ [GeV]"
            elif channel == "mumu":
                xlabel = r"$p_T(\mu)$ [GeV]"
                ylabel = r"$p_T(\mu)$ [GeV]"
        elif mode == "mass":
            xlabel = r"$m_{{\tau\tau}}^{{vis}}$  [GeV]"
            ylabel = r"$m_{{bb}}^{{vis}}$  [GeV]"
        return xlabel, ylabel

    def set_hep_labels(self, channel):
        hep.cms.text('Preliminary', fontsize=40)
        chn_unicodes = {"etau":   r'$bb\: e\tau$',
                        "mutau":  r'$bb\: \mu\tau$',
                        "tautau": r'$bb\: \tau\tau$',
                        "mumu": r'$bb\: \mu\mu$'}
        cat_map = {'baseline': "baseline", 'baseline_boosted': "baseline boosted",
                   'boostedL_pnet': "boosted", 'res1b': "res 1b", 'res2b': "res 2b"}
        sample_header = self.sample.replace('TT', r"$t\bar{t}$")
        hep.cms.lumitext((chn_unicodes[self.channel] + " (" + cat_map[self.category] + ") | " +
                          sample_header + " (" + self.year + ")"),
                         fontsize=24) # r"138 $fb^{-1}$ (13 TeV)"

    def set_lines(self, h, mode):
        xmin, xmax = h.axes[0].edges[0], h.axes[0].edges[-1]
        ymin, ymax = h.axes[1].edges[0], h.axes[1].edges[-1]

        if mode == "trigger":
            if self.channel == "etau":
                if self.year == "2016":
                    plt.plot([26.5, 26.5], [ymin, ymax],
                             c='black', linewidth=10, label=r"single-e + e$\tau$")
                    plt.plot([xmin, 25.5, 25.5], [191., 191., ymax],
                             c='black', linewidth=10, label=r"single-$\tau$")
                    plt.plot([xmin, 25.5, 25.5], [189., 189., ymin],
                             c='black', linewidth=10, label=r"MET")
                elif self.year == "2017" or self.year == "2018":
                    plt.plot([24.5, 24.5, 32.5, 32.5], [ymax, 36., 36., ymin],
                             c='black', linewidth=10, label=r"single-e + e$\tau$")
                    plt.plot([xmin, 23.5, 23.5], [191., 191., ymax],
                             c='red', linewidth=10, label=r"single-$\tau$")
                    plt.plot([xmin, 23.5, 23.5, 31.5, 31.5], [189., 189., 34., 34., ymin],
                             c='deepskyblue', linewidth=10, label="MET")
                
            elif self.channel == "mutau":
                if self.year == "2016":
                    plt.plot([19.5, 19.5, 24.5, 24.5], [ymax, 24.5, 24.5, ymin],
                             c='black', linewidth=10, label=r"single-$\mu$ + $\mu\tau$")
                    plt.plot([xmin, 18.5, 18.5], [191., 191., ymax],
                             c='red', linewidth=10, label=r"single-$\tau$")
                    plt.plot([xmin, 18.5, 18.5, 23.5, 23.5], [189., 189., 31., 31., ymin],
                             c='deepskyblue', linewidth=10, label="MET")
                elif self.year == "2017":
                    plt.plot([20.5, 20.5, 27.5, 27.5], [ymax, 33., 33., ymin],
                             c='black', linewidth=10, label=r"single-$\mu$ + $\mu\tau$")
                    plt.plot([xmin, 19.5, 19.5], [191., 191., ymax],
                             c='red', linewidth=10, label=r"single-$\tau$")
                    plt.plot([xmin, 19.5, 19.5, 26.5, 26.5], [189., 189., 31., 31., ymin],
                             c='deepskyblue', linewidth=10, label="MET")
                elif self.year == "2018":
                    plt.plot([20.5, 20.5, 24.5, 24.5], [ymax, 33., 33., ymin],
                             c='black', linewidth=10, label=r"single-$\mu$ + $\mu\tau$")
                    plt.plot([xmin, 19.5, 19.5], [191., 191., ymax],
                             c='red', linewidth=10, label=r"single-$\tau$")
                    plt.plot([xmin, 19.5, 19.5, 23.5, 23.5], [189., 189., 31., 31., ymin],
                             c='deepskyblue', linewidth=10, label="MET")
         
            elif self.channel == "tautau":
                plt.plot([42., 42., xmax], [ymax, 42., 42.], c='black', linewidth=10, label=r"$\tau\tau$")
                plt.plot([xmin, 39., 39.], [191., 191., ymax], c='red', linewidth=10, label=r"single-$\tau$")
                plt.plot([xmax, 191., 191.], [39., 39., ymin], c='red', linewidth=10)
                plt.plot([189., 189., 39., 39., xmin], [ymin, 39., 39., 189., 189.], c='deepskyblue', linewidth=10, label="MET")
                
            elif self.channel == "mumu":
                if self.year == "2016":
                    plt.plot([24.5, 24.5], [ymax, ymin],
                             c='black', linewidth=10, label=r"single-$\mu$")
                elif self.year == "2017":
                    plt.plot([27.5, 27.5], [ymax, ymin],
                             c='black', linewidth=10, label=r"single-$\mu$")
                elif self.year == "2018":
                    plt.plot([24.5, 24.5], [ymax, ymin],
                             c='black', linewidth=10, label=r"single-$\mu$")

        elif mode == "mass":
            arrow_h = dict(color='red', width=3, head_width=10, head_length=5, shape="full")
            arrow_v = dict(color='red', width=2.6, head_width=7, head_length=6, shape="full")
            line_d = dict(c='red', linewidth=10)
            if self.channel == "mumu":
                if self.dtype == "mc":
                    # plt.plot([20., 130., 130., 20., 20.], [50., 50., 270., 270., 50.],
                    #          label=r"Mass window cut", **line_d)
                    # plt.arrow(20., (ymax-ymin)/2, -4, 0, **arrow_h)
                    # plt.arrow(130., (ymax-ymin)/2, 10, 0, **arrow_h)
                    # plt.arrow(75., 50, 0, -12, **arrow_v)
                    # plt.arrow(75., 270, 0, 12, **arrow_v)

                    dy_l, dy_r = 86., 95.  # DY left and right cuts
                    tt_d, tt_u = 50., 150. # ttbar down and up cuts
                    
                    plt.plot([xmin, dy_l, dy_l], [tt_u, tt_u, ymax], **line_d, label=r"Mass window cut")
                    plt.plot([dy_r, dy_r, xmax], [ymax, tt_u, tt_u], **line_d)
                    plt.plot([xmin, dy_l, dy_l], [tt_d, tt_d, ymin], **line_d)
                    plt.plot([dy_r, dy_r, xmax], [ymin, tt_d, tt_d], **line_d)

                    arrow_length_h, arrow_length_v = 8, 10
                    # top left corner
                    plt.arrow(dy_l/2, tt_u, 0, arrow_length_v, **arrow_v)
                    plt.arrow(dy_l, (ymax-tt_u)/2. + tt_u, -arrow_length_h, 0, **arrow_h)

                    # top right cornerd
                    plt.arrow(dy_r, (ymax-tt_u)/2. + tt_u, arrow_length_h, 0, **arrow_h)
                    plt.arrow((xmax-dy_r)/2.+dy_r, tt_u, 0, arrow_length_v, **arrow_v)

                    # bottom left corner
                    plt.arrow(dy_l, tt_d/2, -arrow_length_h, 0, **arrow_h)
                    plt.arrow(dy_l/2., tt_d, 0, -arrow_length_v, **arrow_v)

                    # bottom right corner
                    plt.arrow(dy_r, tt_d/2, arrow_length_h, 0, **arrow_h)
                    plt.arrow((xmax-dy_r)/2. + dy_r, tt_d, 0, -arrow_length_v, **arrow_v)

                elif self.dtype == "dy":
                    arrow_opt = dict(color='red', width=2, head_width=8, head_length=4, shape="full")
                    plt.plot([86., 86.], [ymin, ymax], label=r"Mass window cut", **line_d)
                    plt.plot([95., 95.], [ymin, ymax], **line_d)

                    plt.arrow(86., (ymax-ymin)/5., 2, 0, **arrow_opt)
                    plt.arrow(95., 2*(ymax-ymin)/5., -2, 0, **arrow_opt)
                    plt.arrow(86., 3*(ymax-ymin)/5., 2, 0, **arrow_opt)
                    plt.arrow(95., 4*(ymax-ymin)/5., -2, 0, **arrow_opt)
                    
                elif self.dtype == "tt":
                    dy_l, dy_r = 86., 95.  # DY left and right cuts
                    tt_d, tt_u = 50., 150. # ttbar down and up cuts
                    
                    plt.plot([xmin, dy_l, dy_l], [tt_u, tt_u, ymax], **line_d, label=r"Mass window cut")
                    plt.plot([dy_r, dy_r, xmax], [ymax, tt_u, tt_u], **line_d)
                    plt.plot([xmin, dy_l, dy_l], [tt_d, tt_d, ymin], **line_d)
                    plt.plot([dy_r, dy_r, xmax], [ymin, tt_d, tt_d], **line_d)

                    arrow_length_h, arrow_length_v = 8, 10
                    # top left corner
                    plt.arrow(dy_l/2, tt_u, 0, arrow_length_v, **arrow_v)
                    plt.arrow(dy_l, (ymax-tt_u)/2. + tt_u, -arrow_length_h, 0, **arrow_h)

                    # top right corner
                    plt.arrow(dy_r, (ymax-tt_u)/2. + tt_u, arrow_length_h, 0, **arrow_h)
                    plt.arrow((xmax-dy_r)/2.+dy_r, tt_u, 0, arrow_length_v, **arrow_v)

                    # bottom left corner
                    plt.arrow(dy_l, tt_d/2, -arrow_length_h, 0, **arrow_h)
                    plt.arrow(dy_l/2., tt_d, 0, -arrow_length_v, **arrow_v)

                    # bottom right corner
                    plt.arrow(dy_r, tt_d/2, arrow_length_h, 0, **arrow_h)
                    plt.arrow((xmax-dy_r)/2. + dy_r, tt_d, 0, -arrow_length_v, **arrow_v)
   
            else:
                plt.plot([20., 130., 130., 20., 20.], [50., 50., 270., 270., 50.],
                         label=r"Mass window cut", **line_d)
            

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="draw cuts on top of signal or MC distributions",
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('--dtype', choices=('signal', 'mc', 'dy', 'tt'), default='signal',
                        type=str, help='Data type')
    parser.add_argument('--skim_tag', required=True, type=str, help='Tag of input skims.')
    parser.add_argument('--signal', choices=('Radion', 'BulkGraviton'), default='Radion',
                        type=str, help='Signal particle type')
    parser.add_argument('--mass', default='1000',
                        type=str, help='Signal particle type')
    parser.add_argument('--channel', required=True, choices=("etau", "mutau", "tautau", "mumu"),
                        type=str, help='Signal particle type')
    parser.add_argument('--category', default="baseline",
                        choices=("baseline", "baseline_boosted", "boostedL_pnet", "res1b", "res2b"),
                        type=str, help='Analysis category')
    parser.add_argument('--year', required=True, choices=("2016", "2016APV", "2017", "2018"),
                        type=str, help='Signal particle type')
    parser.add_argument('--mode', default="trigger", choices=("trigger", "mass"),
                        type=str, help='Signal particle type')
    parser.add_argument('--save', action="store_false",
                        help='Wether to save the histograms or to use the ones produced.')
    FLAGS = parser.parse_args()

    base = "/data_CMS/cms/alves/HHresonant_SKIMS/"
    if FLAGS.dtype == "signal":
        base = os.path.join(base, "SKIMS_UL18_{}_Sig/".format(FLAGS.skim_tag))
        name = os.path.join(base, "GluGluTo{}ToHHTo2B2Tau_M-{}_".format(FLAGS.signal, FLAGS.mass))
        infiles = glob.glob(os.path.join(name, "output_*.root"))
    elif FLAGS.dtype == "dy":
        base = os.path.join(base, "SKIMS_UL18_{}_MC/".format(FLAGS.skim_tag))
        names = ["DYJetsToLL_M-50_TuneCP5_13TeV-amc"]
        infiles = [os.path.join(base, name, "hadded.root") for name in names]
    elif FLAGS.dtype == "tt":
        base = os.path.join(base, "SKIMS_UL18_{}_MC/".format(FLAGS.skim_tag))
        names = ["TTTo2L2Nu", "TTToHadronic", "TTToSemiLeptonic"]
        infiles = [os.path.join(base, name, "hadded.root") for name in names]
    elif FLAGS.dtype == "mc":
        base = os.path.join(base, "SKIMS_UL18_{}_MC/".format(FLAGS.skim_tag))
        names = ["DYJetsToLL_M-50_TuneCP5_13TeV-amc",
                 "TTTo2L2Nu", "TTToHadronic", "TTToSemiLeptonic"]
        infiles = [os.path.join(base, name, "hadded.root") for name in names]

    if FLAGS.dtype == "signal":
        sample = FLAGS.signal + " " + FLAGS.mass + " GeV"
    elif FLAGS.dtype == "dy":
        sample = "DY"
    elif FLAGS.dtype == "tt":
        sample = "TT"
    elif FLAGS.dtype == "mc":
        sample = "TT+DY"

    draw = DrawCuts(infiles, channel=FLAGS.channel, year=FLAGS.year, category=FLAGS.category,
                    sample=sample, dtype=FLAGS.dtype, save=FLAGS.save)
    if FLAGS.mode == "trigger":
        draw.triggers()
    elif FLAGS.mode == "mass":
        draw.mass()
