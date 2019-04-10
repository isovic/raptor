#! /usr/bin/python

import os;
import sys;
import math;

import numpy as np;
from scipy.stats.stats import pearsonr
from scipy.stats.stats import spearmanr
from scipy.optimize import curve_fit

USE_MATPLOTLIB = True;
# try:
	# import matplotlib;
	# matplotlib.use('Agg')
import matplotlib.pyplot as plt;
# plt.use('agg')
from matplotlib.font_manager import FontProperties;

try:
	import seaborn as sns;
	USE_SEABORN = True;
except Exception, e:
	USE_SEABORN = False;
	pass;

# except Exception, e:
# 	USE_MATPLOTLIB = False;

HIGH_DPI_PLOT = False;
# HIGH_DPI_PLOT = True;



def LineFunction(x, b):
    return (1*x + b);

def PlotMedianLine(ax, min_x, max_x, l_median, color='r'):
	x0 = 0;		y0 = 1*x0 + l_median;
	x1 = max_x;	y1 = 1*x1 + l_median;
	ax.plot([x0, x1], [y0, y1], color, lw=1);

def PlotLines(ax, min_x, max_x, l_median, threshold, color='purple'):
	threshold_l = threshold * 2.0 / (math.sqrt(2.0));
	l_min = l_median - threshold_l;
	l_max = l_median + threshold_l;

	x0 = 0;		y0 = 1*x0 + l_min;
	x1 = max_x;	y1 = 1*x1 + l_min;
	ax.plot([x0, x1], [y0, y1], color, lw=1);

	x0 = 0;		y0 = 1*x0 + l_max;
	x1 = max_x;	y1 = 1*x1 + l_max;
	ax.plot([x0, x1], [y0, y1], color, lw=1);



def load_csv(csv_path):
	fp = open(csv_path, 'r');
	lines = fp.readlines();
	fp.close()

	x = [];
	y = [];
	c = [];
	i = 0;
	for line in lines:
		split_line = line.split('\t');
		if (i == 0):
			if (len(split_line) == 6):
				query_header = split_line[0];
				query_id = int(split_line[1]);
				query_length = int(split_line[2]);
				l1_used = True if (int(split_line[3]) == 1) else 0;
				l_median = float(split_line[4]);
				l_maximum_allowed = float(split_line[5]);
			elif (len(split_line) == 7):
				query_header = split_line[0];
				query_id = int(split_line[1]);
				query_length = int(split_line[2]);
				ref_header = split_line[3];
				# ref_id = int(split_line[4]);
				# ref_length = int(split_line[5]);

				l1_used = True if (int(split_line[4]) == 1) else 0;
				l_median = float(split_line[5]);
				l_maximum_allowed = float(split_line[6]);
			else:
				sys.stderr.write('ERROR: Heading line does not contain a valid number of parameters!\n');
				exit(1);

			i += 1;
			continue;
		x.append(float(split_line[0].strip()));
		y.append(float(split_line[1].strip()));
		if (len(split_line) > 2):
			c.append(int(split_line[2].strip()));
		else:
			c.append(0);
		i += 1;
	return [x, y, c, query_header, query_id, query_length, l1_used, l_median, l_maximum_allowed];

def plot_data(fig, ax, subplot_coords, x, y, c, query_length, ymin, ymax, l_median, threshold_L1_under_max, plot_mode, plot_title, out_png_path=''):
	if USE_MATPLOTLIB == True:
		ax.grid();

		ax.text(0.5, 1.02, plot_title,
			horizontalalignment='center',
			fontsize=12,
			transform = ax.transAxes)

		plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0), useOffset=False)
		plt.xlim(0, query_length);
		plt.ylim(ymin, ymax);
		plt.xticks(np.arange(0, query_length, query_length/5))

		all_colors = 'bgrcmyk';

		colors1 = [all_colors[val%len(all_colors)] for val in c[0::2]];
		colors2 = [all_colors[val%len(all_colors)] for val in c[1::2]];

		i = 0;
		while (i < len(x)):
			ax.plot(x[i:(i+2)], y[i:(i+2)], color=all_colors[(c[i]) % len(all_colors)]);
			i += 2;

		ax.scatter(x[0::2], y[0::2], s=10, edgecolor=colors1, facecolor=colors1, lw = 0.2)
		ax.scatter(x[1::2], y[1::2], s=10, edgecolor=colors2, facecolor=colors2, lw = 0.2)

		if (plot_mode == 1):
			try:
				PlotMedianLine(ax, min(x), max(x), l_median, 'r');
				PlotLines(ax, min(x), max(x), l_median, threshold_L1_under_max, 'purple');
				PlotMedianLine(ax, 0, query_length, l_median, 'r');
				PlotLines(ax, 0, query_length, l_median, threshold_L1_under_max, 'purple');
			except Exception, e:
				sys.stderr.write(str(e) + '\n');

		if (subplot_coords == 223 or subplot_coords == 224):
			ax.set_xlabel('Query coordinates');
		if (subplot_coords == 221 or subplot_coords == 223):
			ax.set_ylabel('Reference coordinates');

if __name__ == "__main__":
	if (len(sys.argv) < 2):
		sys.stdout.write('Plots intermediate results from the LCSk-L1 step of the GraphMap algorithm.\n')
		sys.stdout.write('\n')
		sys.stdout.write('Usage:\n')
		sys.stdout.write('  %s <seed_hits.csv> [<seed_hits_2.csv> <seed_hits_3.csv> ...]' % sys.argv[0])
		sys.stdout.write('\n')
		exit(1)

	for seed_hit_path in sys.argv[1:]:
		if (USE_SEABORN == True):
			sns.set_style("darkgrid");
			sns.set_style("white")

		fig = plt.figure()
		ax1 = fig.add_subplot(111)

		sys.stderr.write('Reading: %s\n' % seed_hit_path);
		[x, y, c, query_header, query_id, query_length, l1_used, l_median, l_maximum_allowed] = load_csv(seed_hit_path);
		if (len(x) == 0 or len(y) == 0): continue;
		ymin = min(y);
		ymax = max(y);
		plot_data(fig, ax1, 221, x, y, c, query_length, ymin, ymax, l_median, l_maximum_allowed, 0, query_header if len(query_header) < 30 else query_header[0:30], '');

		out_png_path = '%s.png' % (seed_hit_path);
		if (out_png_path != ''):
			sys.stderr.write('Writing image to file: %s\n\n' % out_png_path);
			if (HIGH_DPI_PLOT == False):
				plt.savefig(out_png_path, bbox_inches='tight');
			else:
				plt.savefig(out_png_path, bbox_inches='tight', dpi=1000);

	plt.show();
