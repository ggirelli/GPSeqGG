#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 1.0.1
# Description: Shuffle a certain percentage of reads in the given bed files.
# 
# ------------------------------------------------------------------------------



# ENV VAR ======================================================================

export LC_ALL=C

# DEPENDENCIES =================================================================

function join_by { local IFS="$1"; shift; echo "$*"; }
IFS=', ' read -r -a array <<< "$string"

# INPUT ========================================================================

# Help string
helps="
 usage: ./beds_shuffle.sh [-h][-n nIter][-p perc][-o outDir] -s seed [BEDFILE]...

 Description:
  Shuffle a certain percentage of reads in the given bed files.

 Mandatory arguments:
  -s seed	Seed for random number generation.
  BEDFILE	Bed file(s). In any order.

 Optional arguments:
  -h	Show this help page.
  -n nIter	Number of iterations. Default: 100
  -p perc	Percentage of reads to shuffle. Default: 10
  -o outDir	Output directory. Default: ./shuffled/
"

# Default values
nIter=100
perc=10
outDir='./shuffled/'

# Parse options
while getopts hn:p:o:s: opt "${bedfiles[@]}"; do
	case $opt in
		h)
			echo -e "$helps\n"
			exit 0
		;;
		n)
			if [ 1 -le $OPTARG]; then
				nIter=$OPTARG
			fi
		;;
		p)
			if [ 1 -le $OPTARG ]; then
				perc=$OPTARG
			fi
		;;
		o)
			outDir=$OPTARG
		;;
		s)
			seed=$OPTARG
		;;
	esac
done

# Check mandatory options
if [ -z "$seed" ]; then
	msg="!!! Missing mandatory -s option."
	echo -e "$helps\n$msg"
	exit
fi

# Read bedfile paths
shift $(($OPTIND - 1))
bedfiles=()
for bf in $*; do
	if [ -e $bf ]; then
		bedfiles+=("$bf")
	else
		msg="!!! Invalid bedfile, file not found.\n    File: $bf"
		echo -e " $helps\n$msg"
		exit 1
	fi
done
if [ 0 -eq ${#bedfiles[@]} ]; then
	msg="!!! No bedfiles provided."
	echo -e " $helps\n$msg"
	exit 1
fi

# TEST =========================================================================

# RUN ==========================================================================

# Create output directory if it does not exist
if [ ! -d $outDir ]; then
	mkdir -p $outDir
fi

# Run single-reshuffling script
for bfi in $(seq 0 `bc <<< "${#bedfiles[@]}-1"`); do
	bf=${bedfiles[$i]}
	if [ 0 -eq $bfi ]; then
		./bed_shuffle.R $seed $bf -n $nIter -p $perc -o $outDir -k T
	else
		./bed_shuffle.R $seed $bf -n $nIter -p $perc -o $outDir
	fi
done

if [ -e $outDir"/.Random.seed.RData" ]; then
	rm $outDir"/.Random.seed.RData"
fi

# End --------------------------------------------------------------------------

################################################################################
