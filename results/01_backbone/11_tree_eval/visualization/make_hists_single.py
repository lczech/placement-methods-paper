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

#def plot_dataset():
	

#######################################################################
#    Dataset Information
#######################################################################

datasets   = [ "General", "Archaea", "Bacteria", "Eukaryota" ]
distance   = "edge_distances"
constraint = "unconstr"
blacklist  = "no-blacklist"
method     = "exact"

if distance == "branch_distances":
	max_x = 1.0
else:
	max_x = 10.0

#######################################################################
#    Plotting Loop
#######################################################################

sns.set(font_scale=1.6)

for dataset in datasets:
	filename = "lists/weighted_placements_" + distance + "_" + constraint + "_" + blacklist + "_" + dataset + ".csv"
	
	if method == "simple":
		plot_simple_hist( filename, max_x, 10 )
	elif method == "histogram":
		plot_histogram( filename, max_x, 10 )
	elif method == "exact":
		plot_exact_hist( filename )

plt.legend(datasets, loc='lower right')

plt.xlabel("Distance (Number of Edges)")
plt.ylabel("Cumulative Frequency")

plt.locator_params( nbins=10, axis='x' )
plt.locator_params( nbins=10, axis='y' )
plt.xlim( 0.0, max_x )
plt.ylim( 0.0, 1.0 )

ax = plt.gca()
ax.set_yticklabels([ '{:3.0f}%'.format(x*100) for x in ax.get_yticks() ])
#ax.set_xticklabels([ '< {:1.0f}'.format(x) for x in ax.get_xticks() ], rotation=45)

plt.tight_layout()
#plt.savefig("weighted_placements_edge_distances_unconstr_hist_10b.pdf", format='pdf')
plt.savefig("weighted_placements_" + distance + "_" + constraint + "_" + blacklist + "_" + method + ".png", format='png')
plt.show()

