#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 0.1.0
# Description: calculates global centrality metrics
# 
# ------------------------------------------------------------------------------



# ENV VAR ======================================================================

export LC_ALL=C

# INPUT ========================================================================

# Help string
helps="
 usage: ./global_centrality.sh [-h] -c csList [BEDFILE]...

 Description:
  Calculate global centrality metrics.

 Mandatory arguments:
  BEDFILE	Bed file(s). Expected to be ordered per condition.
  -c csList	Cutsite list file.

 Optional arguments:
  -h		Show this help page.
"

# Parse options
while getopts hc: opt "${bedfiles[@]}"; do
	case $opt in
		h)
			echo -e "$helps\n"
			exit 0
		;;
		c)
			if [ -e $OPTARG ]; then
				csList=$OPTARG
			else
				msg="!!! Invalid -c option, file not found.\n    File: $OPTARG"
				echo -e "$helps\n$msg"
				exit 1
			fi
		;;
		\?)
		;;
	esac
done

# Check mandatory options
if [ -z "$csList" ]; then
	msg="!!! Missing mandatory -c option."
	echo -e "$helps\n$msg"
	exit
fi

# Read list of bedfile paths
for bf in $bedfiles; do
	echo 1
	# Check that specified bed files exist
	if [ ! -e $bf ]; then
		msg="!!! Specified file not found.\n    File: $bf"
		echo -e "$helps\n$msg"
		exit 1
	fi
done

# RUN ==========================================================================

# Cycle over chromosomes
for chr in $(echo $(seq 1 22) X); do
	echo $chr
done

# mk cumulative over conditions
# 	calc #cs per cchr
# 	cycle over bedfiles
# 		count #reads in the cond in the chr
# 		sum normalized read counts per condition per chr
# 
# build rankings by comparing consecutive conditions
# 	compare with ratio B/A
# Make table with position and rankings
# Sum rankings and sort

# End --------------------------------------------------------------------------

################################################################################
