#!/usr/bin/env python 
import csv
import math
import matplotlib.pyplot as pp
import itertools

# NOTE: This script must be run from the reports directory where the results are

# This class opens a file and iterates over the lines in it, skipping commented-out
# lines
class CommentedFile:
    def __init__(self, f, commentstring="#"):
        self.f = f
        self.commentstring = commentstring
    def next(self):
        line = self.f.next()
        # If the current line is commented, read another line.
        while line.startswith(self.commentstring):
            line = self.f.next()
        return line
    def __iter__(self):
        return self

# Function for creating a plot for some fixed Reynolds number
def make_off_lattice_plot(csv_in, png_out, y_title, kernel, reynold, lattice_name):

        def get_angles(axis_dir):
            y = int(round(180*math.asin(axis_dir[0]) / math.pi))
            x = int(round(180*math.acos(axis_dir[2] / math.cos(y*math.pi/180) )/math.pi ))
            return [x,y]

        def get_summary(data_point):
            axis = [float(x) for x in data_point[4][1:-1].split(',')]
            [x,y] = get_angles(axis)
            return 'X%dY%d' % (x,y)

	raw_data=csv.reader(CommentedFile(open(csv_in,'r')))

        # Filter out the non-data rows
	filtered_data=[x for x in raw_data if len(x)>2 and x[1] != '']

        mapped_kernel = {'EntropicAnsumali' : 'enta',
                         'EntropicChik' : 'entc',
                         'LBGK' : 'lbgk',
                         'MRT' : 'mrt'}

        # Filter out the rows for other kernels
        filtered_data = [x for x in filtered_data if 
          mapped_kernel[x[2]] == kernel.lower() 
          and (float(x[0]) / float(reynold)) < 1.05 
          and (float(x[0]) / float(reynold)) > 0.95
          and get_summary(x) != 'X0Y0']

        # Define a labelling function based on attributes of the results
	def label(datum):
		return "%s %s"%(datum[2],datum[3])

        # Group the data by label
	groupeddata={k:list(v) for k,v in itertools.groupby(sorted(filtered_data,key=label),label)}

        pp.figure(figsize=(16,12))

        order = {}

        for x in range(15, 60, 15):
            for y in range(0, x+15, 15):
                order ['X%iY%i' % (x,y)] = len(order) + 1

        # Add data to the figure, with a log-log plot
	for k in groupeddata:
		cdat=sorted(groupeddata[k],key=get_summary)
                pp.plot([order[get_summary(x)] for x in cdat],[x[1] for x in cdat],'+-',label=k)
        pp.xticks([order[x] for x in order], [x for x in order])
        pp.legend(loc='upper left')
        pp.title('%s over axis orientation for %s' % (y_title, lattice_name))
        pp.xlabel('Axis orientation')
        pp.ylabel(y_title)
        pp.savefig(png_out)

reynolds = [1]

# Plot all graphs for each lattice type
for lattice_n in [15, 19, 27]:
  for kernel in ['Lbgk', 'Mrt', 'EntA', 'EntC']:
    for field in ['velocity']:
      for norm in ['max', 'ave']:
        for reynold in reynolds:
          y_title = '%s normalised residual %s' % ('Average' if norm=='ave' else 'Maximum', field)
          make_off_lattice_plot('%i_lattice/data_files/%s_%s_error_by_re.csv' % (lattice_n, norm, field),
            '%i_lattice_%s_off_lattice_%s_%s_comparison.eps' % (lattice_n, kernel, norm, field),
            y_title,
            kernel,
            reynold,
            'D3Q%i' % lattice_n)
