#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 0.1.0
# Date: 2017-08-17
# Project: GPSeq - bin-size.
# Description: identify smallest optimal bin-size.
# 
# ------------------------------------------------------------------------------



# DEPENDENCIES =================================================================

import argparse
import numpy as np
import pandas as pd
import progressbar
import sys

# PARAMETERS ===================================================================


# Add script description
parser = argparse.ArgumentParser(description = '''
Identify smallest bin-size that provides a restriction probability distribution
comparable with similar bin-sizes. Regions in the bed file are assigned to a bin
if their middle point is inside it. Left-boundary is inclusive.
''')

# Add mandatory arguments
parser.add_argument('bedFile', type = str, nargs = 1,
	help = """GPSeq condition bedfile.""")

# Add arguments with default value
parser.add_argument('--min-bin-size', type = int, nargs = 1,
	metavar = 'binSize', help = """Default to 100e3 bp (100 kbp).""",
	default = [1e5])
parser.add_argument('--max-bin-size', type = int, nargs = 1,
	metavar = 'binSize', help = """Default to 10e6 bp (10 Mbp).""",
	default = [1e7])
parser.add_argument('--bin-size-step', type = int, nargs = 1,
	metavar = 'binSize', help = """Step between bin sizes.
	Default to 100e3 bp (100 kbp).""",
	default = [1e5])
parser.add_argument('--bin-step', type = float, nargs = 1,
	metavar = 'binSize', help = """Percentage of bin size to use for the bin
	step, as moving window. Set to 1 for non-overlapping bins.
	Default to 0.1 (10%% bin size).""",
	default = [0.1])

# Parse arguments
args = parser.parse_args()

# Assign to in-script variables
bedFile = args.bedFile[0]
binSize_min = int(args.min_bin_size[0])
binSize_max = int(args.max_bin_size[0])
binSize_step = int(args.bin_size_step[0])
binStep_perc = args.bin_step[0]

# Additional checks
if 0 == binSize_min:
	sys.exit("!ERROR! Bin size min must be a positive non-zero integer.")
if 0 == binSize_min:
	sys.exit("!ERROR! Bin size max must be a positive non-zero integer.")
if 0 == binSize_min:
	sys.exit("!ERROR! Bin size step must be a positive non-zero integer.")
if 0 == binStep_perc:
	sys.exit("!ERROR! Bin step must be a positive non-zero integer.")

# Show settings
print("""
  binSize : %d - %d ; by %d
  binStep : %.1f%%
 bed file : %s
""" % (binSize_min, binSize_max, binSize_step, binStep_perc * 100, bedFile))

# FUNCTIONS ====================================================================

def get_chr_size(bed):
	'''Get chromosome size from bed file.

	Args:
		bed (pd.DataFrame): bed file content.

	Returns:
		dict: keys are chromosome labels, values are the size.
	'''

	# Empty dictionary
	d = {}

	# Set default values
	for chromosome in set(bed['chr'].tolist()):
		d[chromosome] = 0

	# Update dictionary
	for i in range(bed.shape[0]):
		d[bed['chr'][i]] = max(d[bed['chr'][i]], bed['end'][i])

	# Output
	return(d)

def sum_to_bins(bins, bed, chr_list, chr_sizes):
	'''Assign rows from bed to bins based on middle position. Also, it sorts
	both bed and bins structure to avoid useless comparisons.

	Args:
		bins (pd.DataFrame): bed file with regions of interest.
		bed (pd.DataFrame): bed file with regions to be assigned to ROIs.

	Returns:
		pd.DataFrame: bed file with added rois column.
	'''

	# Add region middle point
	bed['middle'] = bed['start'] + (bed['end'] - bed['start']) / 2.
	midi = np.where('middle' == bed.columns)[0]

	# Setup progress bar
	bar = progressbar.ProgressBar(max_value = bins.shape[0])
	bc = 0

	binned = bins.copy()
	binned['score'] = binned['score'].astype('float')

	# # Per chromosome
	# for c in chr_list:
	# 	# Subset dataframes
	# 	cbed = bed.iloc[np.where(c == bed['chr'])[0],]
	# 	cbins = bins.iloc[np.where(c == bins['chr'])[0],]

	# 	# Sort dataframes
	# 	cbed = cbed.iloc[np.argsort(np.array(cbed['start']).astype('int'))]
	# 	cbins = cbins.iloc[np.argsort(np.array(cbins['start']).astype('int'))]

	# 	for i in cbed.index:
	# 		bar.update(bc)
	# 		bc += 1

	# 		for j in cbins.index:
	# 			if int(cbins.loc[j, 'end']) > cbed.loc[i, 'middle']:
	# 				cbins = cbins.loc[j:,]
	# 				break

	# 		selected = []
	# 		for j in cbins.index:
	# 			if int(cbins.loc[j, 'start']) <= cbed.loc[i, 'middle']:
	# 				selected.append(j)
	# 			else:
	# 				break

	# 		binned.loc[selected, 'score'] = binned.loc[selected, 'score'] + cbed.loc[i, 'score']

	# return(cbins)


	# Per chromosome
	for c in chr_list:
		# Subset dataframes
		cbed = bed.iloc[np.where(c == bed['chr'])[0],]
		cbins = bins.iloc[np.where(c == bins['chr'])[0],]

		# Sort dataframes
		cbed = cbed.iloc[np.argsort(np.array(cbed['start']).astype('int'))]
		cbins = cbins.iloc[np.argsort(np.array(cbins['start']).astype('int'))]

		# Assign to bins
		for i in range(cbins.shape[0]):
			# Update progress status
			bar.update(bc)
			bc += 1

			#print("%d - %d" % (i, cbed.shape[0]))
			bin_start = float(cbins.iloc[i, 1])
			bin_end = float(cbins.iloc[i, 2])

			# Identify shrinking position
			shrink_from = False
			for j in range(cbed.shape[0]):
				if cbed.iloc[j, midi][0] >= bin_start:
					shrink_from = j
					break;

			# If the chromosome bed is over, skip to next chromosome
			if False is shrink_from:
				continue
			else:
				cbed = cbed.iloc[shrink_from:,]

			# Identify contained rows
			for j in range(cbed.shape[0]):
				if cbed.iloc[j, midi][0] <= bin_end:
					binned.iloc[i, 4] = binned.iloc[i, 4] + cbed.iloc[j, 4]
				else:
					break

	# Output
	return(binned)
			




def bin_chr(chr_lab, chr_len, size, step, last_bin):
	'''Generate bins covering a chromosome.

	Args:
		chr_lab (string): chromosome label.
		chr_len (int): chromosome length.
		size (int): bin size.
		step (int): bin step. Use step == size for not overlapping bins.
		last_bin (bool): whether to add extra final bin.

	Returns:
		pd.Dataframe: a chr-start-end-name table.
	'''

	# Calculate bin borders
	starts = np.arange(0, int(chr_len) - size, step)
	ends = starts + size - 1
	scores = np.zeros(starts.shape)

	if last_bin:
		# Add last bin
		starts = starts.tolist()
		starts.append(starts[-1] + size)
		starts = np.array(starts)
		ends = ends.tolist()
		ends.append(ends[-1] + size)
		ends = np.array(ends)
		scores = scores.tolist()
		scores.append(0)
		scores = np.array(scores)

	# Generate bin names
	names = ['bin_' + str(i + 1) for i in range(len(starts))]

	# Produce output bedfile
	out = pd.DataFrame(
		data = np.transpose([np.tile(chr_lab, starts.shape[0]),
			starts, ends, names, scores]),
		index = np.arange(starts.shape[0]),
		columns = ['chr', 'start', 'end', 'name', 'score']
	)

	# Output
	return(out)

def bin_genome(chr_list, chr_sizes, bin_size, bin_step, last_bin):
	'''Bin the whole genome.

	Args:
		chr_list (list): sorted chromosome list.
		chr_sizes (dict): (chr, length in nt) dictionary.
		bin_sizes (int): bin size in nt.
		bin_step (int): bin step in nt.
		last_bin (bool): whether to add extra final bin.

	Returns:
		pd.DataFrame: binned genome.
	'''

	# Empty binned list
	bins = []

	# Per chromosome
	for c in chr_list:
		bins.append(bin_chr(c, chr_sizes[c], bin_size, bin_step, True))

	# Concatenate
	bins = pd.concat(bins)
	bins.index = range(bins.shape[0])

	# Output
	return(bins)

# RUN ==========================================================================

# Read bed file
bed = pd.read_csv(bedFile, "\t", index_col = False, skiprows = [0],
	names = ['chr', 'start', 'end', 'name', 'score'])

# Identify bin-sizes
bin_sizes = range(binSize_min, binSize_max + 1, binSize_step)

# Identify chromosome max
chr_list = np.unique(bed['chr'], return_index = True)
chr_list = [bed['chr'][i] for i in sorted(chr_list[1])]
chr_sizes = get_chr_size(bed)

#for bin_size_id in range(len(bin_sizes)):
for bin_size_id in [0]:

	# Calculate bin characteristics
	bin_size = bin_sizes[bin_size_id]
	bin_step = int(bin_size * binStep_perc)
	print(""" · Building bins...
    Size : %d
    Step : %d""" % (bin_size, bin_step))
	bins = bin_genome(chr_list, chr_sizes, bin_size, bin_step, True)
	print("   > Generated %d bins." % (bins.shape[0],))

	# Assign bed rows to bins
	print(" · Assigning to bins...")
	binned = sum_to_bins(bins, bed, chr_list, chr_sizes)
	binned.to_csv('test.csv')



	

# END ==========================================================================

################################################################################
