#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 0.1.0
# Date: 2017-08-17
# Project: GPSeq - centrality estimation
# Description: estimate genomic region nuclear centrality.
# 
# ------------------------------------------------------------------------------



# DEPENDENCIES =================================================================

import argparse
import numpy as np
import pandas as pd
from StringIO import StringIO
import sys

# PARAMETERS ===================================================================


# Add script description
parser = argparse.ArgumentParser(description = '''
Estimate the nuclear centrality of a genomic region from GPSeq sequencing data.
''')

# Add mandatory arguments
parser.add_argument('bedFile', type = str, nargs = '+',
	help = """At least two (2) GPSeq condition bedfiles, space-separated and in
	increasing order of restriction conditions intensity.""", default = [])

# Add arguments with default value
parser.add_argument('--bin-size', type = int, nargs = 1,
	metavar = 'binSize', help = """Default to chromosome-wide bins.""",
	default = [-1])
parser.add_argument('--bin-step', type = int, nargs = 1,
	metavar = 'binSize', help = """Default to binStep.""",
	default = [-1])

# Parse arguments
args = parser.parse_args()

# Assign to in-script variables
bedFiles = args.bedFile
binSize = args.bin_size[0]
binStep = args.bin_step[0]
chrWide = False

# Additional checks
if 0 == len(bedFiles):
	sys.exit("!ERROR! At least 2 bed files are required.")
if 0 == binSize:
	sys.exit("!ERROR! Bin-size must be a positive non-zero integer.")
if 0 == binStep:
	sys.exit("!ERROR! Bin step must be a positive non-zero integer.")
if -1 == binSize:
	chrWide = True
if not chrWide and -1 == binStep:
	binStep = binSize

# Show settings
print("""
binSize : %s
binStep : %s

bed files:
 %s
""" % (
	"chr-wide" if chrWide else binSize,
	"-" if chrWide else binStep,
	"\n ".join([b for b in bedFiles])
))

# FUNCTIONS ====================================================================

def read_bed(fname):
	'''Read BED file.

	Args:
		fname (string): path to bed file.

	Returns:
		pd.DataFrame: bed file content.
	'''

	# Empty content receiver
	content = ""

	# Line by line
	with open(fname, 'r') as f:
		for line in f:
			# Skip header if present
			if not line.startswith('track'):
				# Strip and split
				content += line

	# Convert to DataFrame
	bed = pd.read_csv(StringIO(content), delimiter = '\t', header = None)
	bed.columns = ['chr', 'start', 'end', 'name', 'score']

	# Output
	return(bed)

def get_chr_size(bed, chr_list):
	'''Get max end position per chromosome.

	Args:
		bed (np.ndarray): bed content.
		chr_list (list): list of chromosome.

	Returns:
		dict: (chr, size) couples.
	'''

	# Empty dictionary receiver
	d = {}

	# Per chromosome
	for chrlab in chr_list:
		d[chrlab] = bed['start'][np.where(chrlab == bed['chr'])[0]].max()

	return(d)

def assign_to_bin(bin_starts, bin_size, bed):
	'''Assign bed scorse to bins.

	Args:
		bin_starts (list): list of bin start positions.
		bed (pd.DataFrame): bed content.

	Returns:
		list: score values grouped per bin.
	'''

	scores = []

	bed_start = 0

	for start in bin_starts:
		print(start)
		# Identify starting bed line
		bed_start = np.where(start <= bed['start'][bed_start:])[0][0]

		# Compare
		selected = np.where(start + bin_size > bed['end'][bed_start:])[0]
		scores.append(bed['score'][bed_start:][selected])

	return(scores)



# RUN ==========================================================================

# Read bed files ---------------------------------------------------------------
beds = []
for bfi in range(len(bedFiles)):
	print(" · Reading '%s' ..." % (bedFiles[bfi],))
	bed = read_bed(bedFiles[bfi])
	bed['bfi'] = bfi
	beds.append(bed)
beds = pd.concat(beds)
beds.index = range(beds.shape[0])
#print(beds)

# Get chromosome sizes ---------------------------------------------------------
print(" · Retrieving chromosome sizes ...")
chr_list = list(set(beds['chr']))
chr_sizes = get_chr_size(beds, chr_list)
#print(chr_sizes)

# Prepare bins -----------------------------------------------------------------
if chrWide:
	print(" · Generating bins (chrWide) ...")
else:
	print(" · Generating bins (size: %d; step: %d) ..." % (binSize, binStep))
bins = {}
for chrlab in chr_sizes.keys():
	if chrWide:
		bins[chrlab] = 0
	else:
		bins[chrlab] = np.array(range(
			0, chr_sizes[chrlab] - binSize + 1, binStep))
#print(bins)

# Assign to bins ---------------------------------------------------------------
print (" · Assigning reads to bins ...")
binned = []
for bfi in range(len(bedFiles)):
	bed = beds.iloc[np.where(bfi == beds['bfi'])[0], :]
	binned.append({})
	if chrWide:
		for chrlab in bins.keys():
			binned[bfi][chrlab] = bed['score'][np.where(chrlab == bed['chr'])[0]]
	else:
		for chrlab in bins.keys():
			print assign_to_bin(bins[chrlab], binSize, bed)
print(binned)

# 2) Calculate centrality estimates --------------------------------------------

# 3) Output --------------------------------------------------------------------

# END ==========================================================================

################################################################################
