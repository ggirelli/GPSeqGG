#!/usr/bin/python
# -*- coding: utf-8 -*-

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 0.1.1
# Description: groups UMIs based on their location
# 
# ------------------------------------------------------------------------------



# DEPENDENCIES =================================================================

import argparse
import sys

# INPUT ========================================================================

# Add script description
parser = argparse.ArgumentParser(
	description = 'Group UMIs based on their location.'
)

# Add params
parser.add_argument('dirpath', metavar = 'dirpath', type = str, nargs = 1,
	help = 'Experiment condition directory, contains UMI counts.')
parser.add_argument('experiment', metavar = 'expID', type = str, nargs = 1,
	help = 'Experiment ID (e.g., TK26).')
parser.add_argument('condition', metavar = 'cond', type = str, nargs = 1,
	help = 'Condition folder name (e.g., 400U2h).')
parser.add_argument('umi_length', metavar = 'umiLen', type = int, nargs = 1,
	help = 'Length of the UMIs in nt (e.g., 8).')

# Add flags
parser.add_argument('--mask-file', metavar = 'MF', type = str, nargs = 1,
	default = [''],
	help = 'File with regions to be masked.')

# Parse arguments
args = parser.parse_args()

# Retrieve arguments
dirpath = args.dirpath[0]
experiment = args.experiment[0]
condition = args.condition[0]
umi_length = args.umi_length[0]
maskfile = args.mask_file[0]

# Set variables ----------------------------------------------------------------

# UMI.pos file columns
chr_col_id = 0
pos_col_id = 1
lseq_col_id = 2
lqual_col_id = 3

# RUN ==========================================================================

print('\nGrouping UMIs per location:')

# Read file --------------------------------------------------------------------
print(' · Retrieving filtered mapped reads ...')
fname = dirpath + condition + '.filtered.umi.pos.txt'
f = open(fname, 'r+')
reads = f.readlines()
f.close()


# Group UMIs per position ------------------------------------------------------
print(' · Grouping UMIs per position ...')

# Will contain grouped UMIs
d = {}

# For every read
for read in reads:
	read = read.strip().split('\t')

	# Identify position
	k = read[chr_col_id] + '~' + read[pos_col_id]

	# Store UMI (seq, qual) at position
	if d.has_key(k):
		d[k][0].append(read[lseq_col_id][0:umi_length])
		d[k][1].append(read[lqual_col_id][0:umi_length])
	else:
		d[k] = [[read[lseq_col_id][0:umi_length]],
			[read[lqual_col_id][0:umi_length]]]

# Output position-grouped UMIs -------------------------------------------------
print(' · Saving grouped UMIs ...')

# Will contain the file content
s = ''

# For every UMi group
for (pos, umi) in d.items():
	# Retrieve position
	c = pos.split('~')

	# Update chr
	if c[0] == 'X':
		c[0] = '23'
	if c[0] == 'Y':
		c[0] = '24'

	# Add sequence and qual to output
	s += '\t'.join(c) + '\t' + ' '.join([str(x) for x in umi[0]]) + '\t' 
	s += ' '.join([str(x) for x in umi[1]]) + '\n'

# Decide filename
if 0 != len(maskfile):
	fname = dirpath + 'UMIpos.all.txt'
else:
	fname = dirpath + 'UMIpos.txt'

# Write output to file
f = open(fname, 'w+')
f.write(s)
f.close()

# Masking UMIs -----------------------------------------------------------------
if 0 != len(maskfile):
	print(' · Masking locations ...')

	# Read masked locations
	print(' >>> Reading maskfile ...')
	f = open(maskfile, 'r+')
	mrf = f.readlines()
	if 0 != len(mrf):
		mrf.pop(0)
	f.close()

	# Restructure
	mrl = {}
	for read in mrf:
		read = read.strip().split('\t')

		# Save regions
		if mrl.has_key(read[0]):
			mrl[read[0]].append([int(read[2]), int(read[3])])
		else:
			mrl[read[0]] = [[int(read[2]), int(read[3])]]

	# Log
	print(' >>> Found ' + str(len(mrf)) + ' regions to mask.')
	print(' >>> Masking ...')

	# Set counter and output variable
	dm = {}
	mc = 0
	for pos in d.keys():
		pos = pos.split('~')

		# Mask based on position
		if mrl.has_key(pos[0]):
			masked = False
			for m in mrl[pos[0]]:
				if int(float(pos[1])) <= m[1] and int(float(pos[1])) >= m[0]:
					masked = True
					break

			if not masked:
				dm['~'.join(pos)] = d['~'.join(pos)]
			else:
				mc += 1
		else:
			dm['~'.join(pos)] = d['~'.join(pos)]

	print(' >>> Masked ' + str(mc) + ' locations.')

	# Output
	print(' · Saving masked UMIs ...')

	# Will contain the file content
	s = ''

	# For every not-masked UMI
	for (pos, umi) in dm.items():
		# Retrieve position
		c = pos.split('~')

		# Update chr
		if c[0] == 'X':
			c[0] = '23'
		if c[0] == 'Y':
			c[0] = '24'

		# Add sequence and qual to output
		s += '\t'.join(c) + '\t' + ' '.join([str(x) for x in umi[0]]) + '\t'
		s += ' '.join([str(x) for x in umi[1]]) + '\n'

	# Write output to file
	fname = dirpath + 'UMIpos.txt'
	f = open(fname, 'w+')
	f.write(s)
	f.close()

# END --------------------------------------------------------------------------

################################################################################
