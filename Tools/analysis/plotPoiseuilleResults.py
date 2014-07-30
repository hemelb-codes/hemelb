#!/usr/bin/env python 
import csv
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
def make_Re_plot_for_kernel(csv_in, png_out, y_title, kernel, lattice_name):
	raw_data=csv.reader(CommentedFile(open(csv_in,'r')))

        # Filter out the non-data rows
	filtered_data=[x for x in raw_data if len(x)>2 and x[1] != '' and x[4] == str((0.,0.,1.))]

        mapped_kernel = {'EntropicAnsumali' : 'enta',
                         'EntropicChik' : 'entc',
                         'LBGK' : 'lbgk',
                         'MRT' : 'mrt'}
 
        # Filter out the rows for other kernels
        filtered_data = [x for x in filtered_data if mapped_kernel[x[2]] == kernel.lower()]

        # Define a labelling function based on attributes of the results
	def label(datum):
		return "%s %s"%(datum[2],datum[3])

        # Group the data by label
	groupeddata={k:list(v) for k,v in itertools.groupby(sorted(filtered_data,key=label),label)}

        pp.figure(figsize=(16,12))

        # Add data to the figure, with a log-log plot
	for k in groupeddata:
		cdat=groupeddata[k]
                pp.loglog([x[0] for x in cdat],[x[1] for x in cdat],'+-',label=k)
        pp.legend(loc='upper left')
        pp.title('%s over Reynolds number for %s' % (y_title, lattice_name))
        pp.xlabel('Reynolds number')
        pp.ylabel(y_title)
        pp.savefig(png_out)

# Function for creating a plot for some fixed Reynolds number
def make_Re_plot_for_bc(csv_in, png_out, y_title, bc, lattice_name):
        raw_data=csv.reader(CommentedFile(open(csv_in,'r')))

        # Filter out the non-data rows
        filtered_data=[x for x in raw_data if len(x)>2 and x[1] != '' and x[4] == str((0.,0.,1.))]

        mapped_bc = {'FINTERPOLATION' : 'fin',
                         'JUNKYANG' : 'jyg',
                         'SIMPLEBOUNCEBACK' : 'sbb',
                         'GZS' : 'gzs'}

        # Filter out the rows for other kernels
        filtered_data = [x for x in filtered_data if mapped_bc[x[3]] == bc.lower()]

        # Define a labelling function based on attributes of the results
        def label(datum):
                return "%s %s"%(datum[2],datum[3])

        # Group the data by label
        groupeddata={k:list(v) for k,v in itertools.groupby(sorted(filtered_data,key=label),label)}

        pp.figure(figsize=(16,12))

        # Add data to the figure, with a log-log plot
        for k in groupeddata:
                cdat=groupeddata[k]
                pp.loglog([x[0] for x in cdat],[x[1] for x in cdat],'+-',label=k)
        pp.legend(loc='upper left')
        pp.title('%s over Reynolds number for %s' % (y_title, lattice_name))
        pp.xlabel('Reynolds number')
        pp.ylabel(y_title)
        pp.savefig(png_out)


# Plot all graphs for each lattice type
for lattice_n in [15, 19, 27]:
  for kernel in ['Lbgk', 'Mrt', 'EntA', 'EntC']:
    for field in ['pressure', 'velocity']:
      for norm in ['max', 'ave']:
        y_title = '%s normalised residual %s' % ('Average' if norm=='ave' else 'Maximum', field)
        make_Re_plot_for_kernel('%i_lattice/data_files/%s_%s_error_by_re.csv' % (lattice_n, norm, field),
          '%i_lattice_%s_%s_%s_comparison.eps' % (lattice_n, kernel, norm, field),
          y_title,
          kernel,
          'D3Q%i' % lattice_n)
  for bc in ['Sbb', 'Jyg', 'Gzs', 'Fin']:
    for field in ['pressure', 'velocity']:
      for norm in ['max', 'ave']:
        y_title = '%s normalised residual %s' % ('Average' if norm=='ave' else 'Maximum', field)
        make_Re_plot_for_bc('%i_lattice/data_files/%s_%s_error_by_re.csv' % (lattice_n, norm, field),
          '%i_lattice_%s_%s_%s_comparison.eps' % (lattice_n, bc, norm, field),
          y_title,
          bc,
          'D3Q%i' % lattice_n)
