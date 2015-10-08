#!/usr/bin/env python 
import csv
import matplotlib.pyplot as pp
import itertools
import re
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
def make_velocity_plot(csv_in, png_out):
	raw_data=csv.reader(CommentedFile(open(csv_in,'r')))

        # Filter out the non-data rows
	filtered_data=[x for x in raw_data if len(x)>2 and x[1] != '']
        formatted_data = [re.sub('[\[\,\]]*','',x) for x in filtered_data[0]]

        pp.figure(figsize=(16,12))

        # Plot.
        xs = [float(x) for x in formatted_data[0].split()]
        y_theor = [float(x) for x in formatted_data[1].split()]
        y_sim = [float(x) for x in formatted_data[2].split()]

        # Set x axis to have corresponding lattice types
        pp.scatter(xs, y_theor, label='Theoretical velocity profile', c='r', s=60, marker='+')
        pp.scatter(xs, y_sim, label='Simulated velocity profile', c='b', s=60, marker='x')

        pp.xticks([0.000, 0.001,0.002,0.003,0.004])
        pp.title('Theoretical and simulated velocity profiles')
        pp.xlabel('Radius (m)')
        pp.ylabel('Velocity (m/s)')
        pp.legend(loc='upper right')
        pp.savefig(png_out)

def make_pressure_plot(csv_in, png_out):
        raw_data=csv.reader(CommentedFile(open(csv_in,'r')))

        # Filter out the non-data rows
        filtered_data=[x for x in raw_data if len(x)>2 and x[1] != '']
        formatted_data = [re.sub('[\[\,\]]*','',x) for x in filtered_data[0]]

        pp.figure(figsize=(16,12))

        # Plot.
        xs = [float(x) for x in formatted_data[0].split()]
        y_theor = [float(x) for x in formatted_data[1].split()]
        y_sim = [float(x) for x in formatted_data[2].split()]

        # Set x axis to have corresponding lattice types
        pp.scatter(xs, y_theor, label='Theoretical pressure profile', c='r', s=60, marker='+')
        pp.scatter(xs, y_sim, label='Simulated pressure profile', c='b', s=60, marker='x')

        pp.autoscale

        locs,labels = pp.yticks()
        pp.yticks(locs, map(lambda x: "%.5f" % x, locs))

        pp.title('Theoretical and simulated pressure profiles along cylinder')
        pp.xlabel('Z (m)')
        pp.ylabel('Pressure (Pa)')
        pp.legend(loc='upper right')
        pp.savefig(png_out)

make_velocity_plot('poiseuille_compare/data_files/velocity_comparison.csv', 'sample_velocity_profile.eps')
make_pressure_plot('poiseuille_compare/data_files/pressure_comparison.csv', 'sample_pressure_profile.eps')
