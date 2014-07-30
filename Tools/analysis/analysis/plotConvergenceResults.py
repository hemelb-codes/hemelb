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
def make_convergence_plot(csv_in, png_out, y_title, kernel, reynold, lattice_name):
	raw_data=csv.reader(CommentedFile(open(csv_in,'r')))

        # Filter out the non-data rows
	filtered_data=[x for x in raw_data if len(x)>2 and x[1] != '']

        mapped_kernel = {'EntropicAnsumali' : 'enta',
                         'EntropicChik' : 'entc',
                         'LBGK' : 'lbgk',
                         'MRT' : 'mrt'}
 
        # Filter out the rows for other kernels
        filtered_data = [x for x in filtered_data if 
          mapped_kernel[x[3]] == kernel.lower() 
          and (float(x[1]) / reynold) < 1.05 
          and (float(x[1]) / reynold) > 0.95
          and lattice_name == x[4]]

        # Define a labelling function based on attributes of the results
	def label(datum):
		return "%s %s"%(datum[3],datum[5])

        # Group the data by label
	groupeddata={k:list(v) for k,v in itertools.groupby(sorted(filtered_data,key=label),label)}

        pp.figure(figsize=(16,12))

        # Add data to the figure, with a log-log plot
	for k in groupeddata:
		cdat=groupeddata[k]
                pp.plot([ int(round(8e-3/float(x[0]))) for x in cdat],[x[2] for x in cdat],'+-',label=k)
        pp.legend(loc='upper left')
        pp.title('%s over resolution for %s, %s' % (y_title, lattice_name, kernel))
        pp.xlabel('Resolution')
        pp.ylabel(y_title)
        pp.savefig(png_out)

reynolds = [10,100]

# Plot all graphs for each lattice type
for lattice_n in [15, 19, 27]:
  for kernel in ['Lbgk', 'Mrt', 'EntA', 'EntC']:
    for field in ['velocity']:
      for norm in ['max', 'ave']:
        for reynold in reynolds:
          y_title = '%s normalised residual %s' % ('Average' if norm=='ave' else 'Maximum', field)
          make_convergence_plot('convergence/data_files/%s_%s_error_by_resolution.csv' % (norm, field),
            '%i_lattice_%s_%s_%s_convergence.eps' % (lattice_n, kernel, norm, field),
            y_title,
            kernel,
            reynold,
            'D3Q%i' % lattice_n)
