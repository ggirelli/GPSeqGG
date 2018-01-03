#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 1.0.0
# Description: filter alignment output (SAM file)
# 
# ------------------------------------------------------------------------------



# ENV VAR ======================================================================

export LC_ALL=C

# INPUT ========================================================================

# Help string
helps="
 usage: ./sam_filter.sh [-h] -i samfile

 Description:
  Run FASTQC tool.

 Mandatory arguments:
  -i samfile	Input SAM file to be filterd.

 Optional arguments:
  -h		Show this help page.
  -p 		Keep only primary alignments.
  -1 		Keep only R1, for PE-seq.
  -c 		Remove chimeric reads, for PE-seq.
  -m 		Keep only mapped reads.
  -q qual 	Keep only alignments with MAPQ >= qual.
  -C chrlist 	Comma-separated list of chromosomes to be removed.
"

# Default values
keep_mapped=false
keep_r1=false
keep_primary=false
rm_chimeric=false
min_mapq=0

# Parse options
while getopts hp1cmi:q:C: opt; do
	case $opt in
		h)
			echo -e "$helps\n"
			exit 1
		;;
		p) # Keep primary alignments only
			keep_primary=true
		;;
		1) # Keep R1 only
			keep_r1=true
		;;
		c) # Remove chimeras
			rm_chimeric=true
		;;
		m) # Keep only mapped
			keep_mapped=true
		;;
		i) # Input SAM file
			if [ -f "$OPTARG" ]; then
				samfile="$OPTARG"
			else
				msg="!ERROR! Invalid value for option -f."
				msg="$msg\n        File not found: $OPTARG."
				echo -e "$msg"
				exit 0
			fi
		;;
		q) # MAPQ filter
			if [ $OPTARG -ge 0 ]; then
				min_mapq=$OPTARG
			else
				msg="!ERROR! Invalid value for option -q."
				msg="$msg\n        qual must be a positive integer."
				echo -e "$msg"
				exit 0
			fi
		;;
		C) # Comma-separated list of chromosomes
			#...
		;;
	esac
done

# Check mandatory options


# RUN ==========================================================================



# END --------------------------------------------------------------------------

################################################################################
