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

def createHisto(x, y, inputs, avars=None):
    avars = (x, y) if avars is None else (x, y, *avars)
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
        histogram.fill(batch.dau1_pt, batch.dau2_pt)

    return histogram

def drawCuts(inputs, dtype, sample):
    xvar, yvar = "dau1_pt", "dau2_pt"
    histogram = createHisto(x=xvar, y=yvar, inputs=inputs)

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
    FLAGS = parser.parse_args()

    path = {'signal': "/data_CMS/cms/alves/HHresonant_SKIMS/SKIMS_UL18_validateMETNoSF_Sig_V2/",
            'mc'    : "/data_CMS/cms/alves/HHresonant_SKIMS/SKIMS_UL18_validateMETNoSF_MC_V2/"}
    name = "GluGluTo{}ToHHTo2B2Tau_M-{}_".format(FLAGS.signal, FLAGS.mass)
    
    infiles = glob.glob(os.path.join(path[FLAGS.dtype], name, "output_*.root"))

    drawCuts(infiles, dtype=FLAGS.dtype,
             sample=FLAGS.signal + " " + FLAGS.mass + " GeV")
