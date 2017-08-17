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

# RUN ==========================================================================

# 1) Bin bed files -------------------------------------------------------------

# 2) Calculate centrality estimates --------------------------------------------

# 3) Output --------------------------------------------------------------------

# END ==========================================================================

################################################################################
