import os
import argparse
import json
import numpy as np

import matplotlib
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.style.ROOT)

def compare_gains(chn):    
    data = {}
    for mode in ('standard', 'bigtau'):
        with open('data_' + mode + '.json') as json_data:
            data[mode] = json.load(json_data)

    masses = [300, 400, 500, 600, 700, 750, 800, 850, 900, 1000, 1250, 1500, 1750, 2000, 2500, 3000]
    colors = iter(("red", "dodgerblue"))
     
    fig, (ax1, ax2) = plt.subplots(2, sharex=True, gridspec_kw={'height_ratios': [3., 1.]})
    plt.subplots_adjust(wspace=0, hspace=0)

    ax1.set_ylabel("Number of weighted events", fontsize=20)
    for mode in ('standard', 'bigtau'):
        err = [e/2. for e in data[mode]['errs'][chn]]
        ax1.errorbar(masses, data[mode]['vals'][chn],
                     yerr=(err, err), fmt='-o', color=next(colors), label=mode)
    ax1.legend(loc="upper right")
    
    ax2.set_ylim(-0.02, 0.14)
    ax2.set_ylabel(r"(Standard - Bigtau) - 1", fontsize=20)
    ax2.set_xlabel(r"$M(X)\: [GeV]$", fontsize=21)

    erratio = [1/2 * np.sqrt(ea**2/a**2 + eb**2/b**2) * a/b
               for a,b,ea,eb in zip(data['standard']['vals'][chn],data['bigtau']['vals'][chn],
                                    data['standard']['errs'][chn],data['bigtau']['errs'][chn])]
    ax2.errorbar(masses, [x/y-1. for x, y in zip(data['standard']['vals'][chn],data['bigtau']['vals'][chn])],
                 yerr=[erratio,erratio], fmt='-o', color='black')
     
    hep.cms.text(' Preliminary', fontsize=22, ax=ax1)
    hep.cms.lumitext(r"59.7 $fb^{-1}$ (13 TeV)", fontsize=21, ax=ax1)
     
    output = os.path.join("/eos/home-b/bfontana/www/TriggerScaleFactors/CompareRatios/",
                          chn + "_" + os.path.basename(__file__[:-3]))
    for ext in ('.png', '.pdf'):
        fig.savefig(output + ext)
        print('Plot saved under {}'.format(output + ext))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Compare two selection regions with and without bigtau flag.')
    parser.add_argument('--channel', required=True, choices=('etau', 'mutau', 'tautau'), help='Which analysis channel to consider.')

    FLAGS = parser.parse_args()

    compare_gains(chn=FLAGS.channel)
                  
