#!/usr/bin/python

# libraries
import matplotlib
matplotlib.use('TkAgg')

import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import numpy as np
import pandas as pd
import seaborn as sns

#######################################################################
#    Plot Functions
#######################################################################

# Simple histogram test
# Not good: does not normalize correctly,
# as it only takes the values in the range into account
def plot_simple_hist( list_filename, max_x, num_bins ):
	data = pd.read_csv( list_filename, header=None )
	plt.hist( data.iloc[:,0], bins=num_bins, normed=1, range=(0,max_x), histtype='step', cumulative=True, linewidth=1 )

# Plot a proper histogram, normalized by all values
def plot_histogram( list_filename, max_x, num_bins ):
	data = pd.read_csv( list_filename, header=None )
	hist, bins = np.histogram( data.iloc[:,0], bins = np.linspace( 0.0, max_x, num_bins + 1 ) )
	cum = np.cumsum(hist) / float(len(data.iloc[:,0]))

	# the step plot does not draw the last line segement.
	# duplicate it, so that we can see it in the plot.
	cum = np.append( cum, cum[-1] )
	plt.step(bins[:], cum, where='post')

# Plot the exact curve instead of a histogram approximation
def plot_exact_hist( list_filename ):
	# Load list of distances, sort them
	data = pd.read_csv( list_filename, header=None )
	data = np.sort( data.iloc[:,0] )

	# For each distance, add a step on the y-axis,
	# normalized by how many values there are.
	steps = np.array(range(len(data))) / float(len(data))

	# We want to avoid the ugly jump in the beginning.
	# Set all steps in the beginning to the height at the last 0 step.
	last_zero_idx=0
	for i in range(len(data)):
		if data[i] == 0:
			last_zero_idx=i
		else:
			break
	for i in range(last_zero_idx):
		steps[i] = steps[last_zero_idx]

	# Plot the data
	plt.plot(data, steps)

def plot_dataset(datasets, distance, constraint, blacklist, method):
	if distance == "branch_distances":
		max_x = 1.0
	elif distance == "edge_distances":
		max_x = 10.0
	else:
		raise error

	# Plot the graph.
	for dataset in datasets:
		filename = "lists/weighted_placements_" + distance + "_" + constraint + "_" + blacklist + "_" + dataset + ".csv"

		if method == "simple":
			plot_simple_hist( filename, max_x, 10 )
		elif method == "histogram":
			plot_histogram( filename, max_x, 10 )
		elif method == "exact":
			plot_exact_hist( filename )
		else:
			raise error

	# Set texts.
	plt.legend(datasets, loc='lower right')
	if distance == "branch_distances":
		plt.xlabel("Distance (Branch Length Units)")
	elif distance == "edge_distances":
		plt.xlabel("Distance (Number of Branches)")
	else:
		raise error
	plt.ylabel("Cumulative Frequency")

	# Set 10 ticks on x and y axis, and set the desired axis ranges.
	plt.locator_params( nbins=10, axis='x' )
	plt.locator_params( nbins=10, axis='y' )
	plt.xlim( 0.0, max_x )
	plt.ylim( 0.0, 1.0 )

	# Set y tick labels every 20% and grid lines every 10%.
	# Rename labels to be percent.
	ax = plt.gca()
	ax.set_yticks( np.linspace( 0.0, 1.0, 6 ), minor=False)
	ax.set_yticks( np.linspace( 0.1, 0.9, 5 ), minor=True)
	ax.yaxis.grid(True, which='major')
	ax.yaxis.grid(True, which='minor')
	ax.set_yticklabels([ '{:3.0f}%'.format(x*100) for x in ax.get_yticks() ])
	#ax.set_xticklabels([ '< {:1.0f}'.format(x) for x in ax.get_xticks() ], rotation=45)

	plt.tight_layout()
	#plt.figure().tight_layout()

	# Save as PDF with not much white border around it, and as PNG with more border.
	plt.savefig("figures/weighted_placements_" + distance + "_" + constraint + "_" + blacklist + "_" + method + ".pdf", format='pdf', bbox_inches='tight')
	plt.savefig("figures/weighted_placements_" + distance + "_" + constraint + "_" + blacklist + "_" + method + ".png", format='png')
	#plt.show()

#######################################################################
#    Dataset Information
#######################################################################

datasets    = [ "General", "Archaea", "Bacteria", "Eukaryota" ]
distances   = [ "edge_distances", "branch_distances" ]
constraints = [ "unconstr", "constr" ]
blacklists  = [ "no-blacklist", "blacklist" ]
methods     = [ "exact", "histogram" ]

#######################################################################
#    Plotting Loop
#######################################################################

sns.set(font_scale=1.6)

plotnum = 0
for distance in distances:
	#for constraint in constraints:
		#for blacklist in blacklists:
			for method in methods:
				constraint = "unconstr"
				blacklist = "no-blacklist"

				plotnum += 1
				plt.figure( plotnum )

				plot_dataset(datasets, distance, constraint, blacklist, method)

#plt.show()
