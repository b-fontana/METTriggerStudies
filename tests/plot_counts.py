from pathlib import Path
import csv
import bokeh
from bokeh.plotting import figure, save
from bokeh.layouts import gridplot
#from bokeh.io import export_svg

base_dir = '/eos/user/b/bfontana/www/TriggerScaleFactors/'
#main_dir = 
main_dir = ['TriggerStudy_MET200_SingleTau190',
            'TriggerStudy_MET200_SingleTau190_CUT_entries_ditau_deltaR_GT_0_5',
            'TriggerStudy_MET200_SingleTau190_CUT_entries_ditau_deltaR_ST_0_5']

table_dir = {'met': 'counts_total_baseline_met',
             'tau': 'counts_total_baseline_tau',
             'met_tau': 'counts_total_baseline_met_tau'}
channels = ('etau', 'mutau', 'tautau')
#alphas_cuts = {200: 1, 180: 0.5, 220: 0.5}

x = [400, 500, 600, 700, 800, 900, 1000, 1250, 1500, 1750, 2000, 3000]
linear_x = [k for k in range(1,len(x)+1)]

y = {}
for md in main_dir:
    y[md] = {}
    for chn in channels:
        y[md][chn] = {}
        for ktd, vtd in table_dir.items():
            y[md][chn][ktd] = []
            d = Path(base_dir) / md / chn / vtd
            with open(Path(d) / 'table.csv') as f:
                reader = csv.reader(f, delimiter=',', quotechar='|')
                next(reader, None) #ignore header line
                for sample, frac in reader:
                    y[md][chn][ktd].append(float(frac))

for md in main_dir:
    for chn in channels:
        y[md][chn]['add_met_tau'] = []
        for elem1, elem2 in zip(y[md][chn]['met'], y[md][chn]['met_tau']):
            y[md][chn]['add_met_tau'].append(elem1+elem2)

opt_points = dict(size=8)
opt_line = dict(width=1.5)
colors = ('green', 'blue', 'red')
styles = ('solid', 'dashed')
legends = {'met': ' (MET)',
           'tau': ' (Tau)',
           'add_met_tau': ' (MET + Tau)'}

x_str = [str(k) for k in x]
xticks = linear_x[:]
yticks = [x for x in range(0,100,10)]

for md in main_dir:
    p1 = figure(title='Inclusion of MET and Single Tau triggers',
                plot_width=800, plot_height=400,
                x_axis_label='x', y_axis_label='y')
    p2 = figure(title='Acceptance gain of MET + SingleTau triggers',
                plot_width=800, plot_height=400,
                x_axis_label='x', y_axis_label='y')
    for ichn,chn in enumerate(channels):
        for itd,td in enumerate(('met', 'tau')):
            p1.circle(linear_x, y[md][chn][td], color=colors[ichn], fill_alpha=1., **opt_points)
            p1.line(linear_x, y[md][chn][td], color=colors[ichn], line_dash=styles[itd],
                    legend_label=chn+legends[td], **opt_line)
        for itd,td in enumerate(('met', 'add_met_tau')):
            p2.circle(linear_x, y[md][chn][td], color=colors[ichn], fill_alpha=1., **opt_points)
            p2.line(linear_x, y[md][chn][td], color=colors[ichn], line_dash=styles[itd],
                    legend_label=chn+legends[td], **opt_line)
            
    for p in (p1, p2):
        p.xaxis[0].ticker = xticks
        p.xgrid[0].ticker = xticks
        p.xgrid.grid_line_alpha = 0.2
        p.xgrid.grid_line_color = 'black'
        p.yaxis[0].ticker = yticks
        p.ygrid[0].ticker = yticks
        p.ygrid.grid_line_alpha = 0.2
        p.ygrid.grid_line_color = 'black'
         
        p.legend.location = 'top_left'
        p.xaxis.axis_label = 'mHH [GeV]'
        p.yaxis.axis_label = 'Trigger acceptance gain (w.r.t. trigger baseline) [%]'
     
        p.xaxis.major_label_overrides = dict(zip(linear_x,x_str))
        
        p.output_backend = 'svg'
        #export_svg(p, filename='line_graph.svg')
     
    g = gridplot([[p1], [p2]])
    save(g, title=md, filename='plots/'+md+'.html')
