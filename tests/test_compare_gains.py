import os
import argparse
import json
import numpy as np

import matplotlib
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.style.ROOT)

def compare_gains(chn, spin):
    data = {}
    for mode in ('standard', 'bigtau'):
        with open('data_' + chn + '_' + mode + '.json') as json_data:
            data[mode] = json.load(json_data)
        
    chn_unicodes = {"etau":   r'$bb\: e\tau$',
                    "mutau":  r'$bb\: \mu\tau$',
                    "tautau": r'$bb\: \tau\tau$',
                    "mumu":   r'$bb\: \mu\mu$'}
    spin_map = {'0': "Radion", '2': "BulkGraviton"}
    masses = [300, 400, 500, 600, 700, 750, 800, 850, 900, 1000,
              1250, 1500, 1750, 2000, 2500, 3000]
    colors = iter(("red", "dodgerblue"))
     
    fig, (ax1, ax2) = plt.subplots(2, sharex=True, gridspec_kw={'height_ratios': [3., 1.]})
    plt.subplots_adjust(wspace=0, hspace=0)

    ax1.set_ylabel("Number of weighted events", fontsize=20)
    for mode in ('standard', 'bigtau'):
        err = [e/2. for e in data[mode]['errs'][chn]]
        ax1.errorbar(masses, data[mode]['vals'][chn],
                     yerr=(err, err), fmt='-o', color=next(colors), label=mode)
    ax1.legend(loc="upper right")
    
    ax2.set_ylim(-2, 15)
    yticks = [0., 5., 10.]
    ax2.set_yticks(yticks)
    line_opt = dict(color="grey", linestyle="--")
    for yval in yticks:
        ax2.axhline(y=yval, **line_opt)
    ax2.set_ylabel(r"Ratio - 1 [%]", fontsize=20)
    ax2.set_xlabel(r"$m(X)\:\:[GeV]$", fontsize=21)

    erratio = [1/2 * np.sqrt(ea**2/a**2 + eb**2/b**2) * a/b * 100
               for a,b,ea,eb in zip(data['standard']['vals'][chn],data['bigtau']['vals'][chn],
                                    data['standard']['errs'][chn],data['bigtau']['errs'][chn])]
    ax2.errorbar(masses, [(x/y-1.)*100 for x, y in zip(data['standard']['vals'][chn],data['bigtau']['vals'][chn])],
                 yerr=[erratio,erratio], fmt='-o', color='black')
     
    hep.cms.text(' Preliminary', fontsize=22, ax=ax1)
    hep.cms.lumitext(chn_unicodes[chn] + " (baseline) | " + spin_map[spin],
                     fontsize=21, ax=ax1)
     
    output = os.path.join("/eos/home-b/bfontana/www/TriggerScaleFactors/CompareRatios/",
                          chn + "_" + os.path.basename(__file__[:-3]))
    for ext in ('.png', '.pdf'):
        fig.savefig(output + ext)
        print('Plot saved under {}'.format(output + ext))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Compare two selection regions with and without bigtau flag.')
    parser.add_argument('--channel', required=True, choices=('etau', 'mutau', 'tautau'), help='Which analysis channel to consider.')
    parser.add_argument('--spin', required=True, choices=('0', '2'),
                        help='Signal spin hypothesis.')

    FLAGS = parser.parse_args()

    compare_gains(chn=FLAGS.channel, spin=FLAGS.spin)
                  
