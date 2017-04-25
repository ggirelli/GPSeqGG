#!/usr/bin/python
# -*- coding: utf-8 -*-

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 1.0.1
# Description: compare original global centrality output and shuffled ones.
# 
# ------------------------------------------------------------------------------



# DEPENDENCIES =================================================================

import argparse
import numpy as np
import os
import pandas as pd

# INPUT ========================================================================

# Add script description
parser = argparse.ArgumentParser(
	description = 'Compare global centrality outputs.'
)

# Add params
parser.add_argument('oriGC', type = str, nargs = 1,
	help = 'Path to proper global centrality output.')
parser.add_argument('shuffledDir', type = str, nargs = 1,
	help = 'Path to folder with shuffled orders.')
parser.add_argument('outFile', type = str, nargs = 1,
	help = 'Path to output file name.')

# Add flags
parser.add_argument('-p', metavar = 'perc', type = int, nargs = 1,
	default = [10],
	help = 'Percentage of reshuffled reads. Default: 10')

# Parse arguments
args = parser.parse_args()

# Retrieve arguments
proper = args.oriGC[0]
inDir = args.shuffledDir[0]
outName = args.outFile[0]
perc = args.p[0]

# Default order
default = ['chr' + str(i) for i in range(1, 23)]
default.append('chrX')

# Proper order
proper = pd.read_csv(proper, '\t')
cols = proper.columns.values

# FUNCTIONS ====================================================================

def str2list(s):
	'''Convert string to list of characters.

	Args:
		s (strings)

	Returns:
		list: if the input was a string, its conversion to list.
		any: the original input if it was not a string.
	'''
	if type('') == type(s):
		return [c for c in s]
	return s

def unique(l):
	'''Remove duplicated elements from a list.

	Args:
		l (list): list to be de-duplicated.

	Returns:
		list: l without duplicated elements.
	'''
	return [e for e in set(l)]

def dKendall(l1, l2):
	'''Calculate the Kendall tau distance between two lists of
	elements, as the number of swaps needed to make them identical.
	(identical) 0 <= d <= 1 (totally different)

	Args:
		l1 (list, string): first list.
		l2 (list, string): second list.

	Returns:
		float: the Kendall tau distance between l1 and l2.
	'''

	# Convert strings to lists
	l1 = str2list(l1)
	l2 = str2list(l2)

	# Lists MUST be lists
	if any([not type(l) == type([]) for l in [l1, l2]]):
		return None


	# The lists MUST have the same length
	if len(l1) != len(l2):
		msg = 'The Kendall tau distance is not defined '
		msg += 'for sets of different lengths.'
		print(msg)
		return None

	# They must have the same number of elements
	if not set(l1) == set(l2):
		msg = 'The Kendall tau distance is not defined '
		msg += 'for sets of different elements.'
		print(msg)
		return None

	# Check list 2 ranking
	ranked = [l1.index(e) for e in l2]

	# Count swaps
	steps = 0
	for i in range(len(ranked)):
		for j in range(i + 1, len(ranked)):
			if ranked[i] > ranked[j]:
				steps += 1

	# Output normalization
	return steps / ( len(l1) * (len(l1) - 1 ) / 2.)

# RUN ==========================================================================

flist = [f for f in os.listdir(inDir) if os.path.isfile(os.path.join(inDir, f))]
flist = [f for f in flist if '%dperc.recap.txt' % (perc,) in f]

cmp_matrix = pd.DataFrame(index = np.arange(len(flist)), columns = cols)

for fi in range(len(flist)):
	fname = flist[fi]

	t = pd.read_csv(inDir + '/' + fname, '\t')
	for ki in range(len(cols)):
		k = cols[ki]
		cmp_matrix[k][fi] = dKendall(proper[k].tolist(), t[k].tolist())

cmp_matrix.to_csv(outName, sep = '\t', header = True, index = False)

# END --------------------------------------------------------------------------

################################################################################
