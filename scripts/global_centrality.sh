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

# DEPENDENCIES =================================================================

function join_by { local IFS="$1"; shift; echo "$*"; }

# INPUT ========================================================================

# Help string
helps="
 usage: ./global_centrality.sh [-h] -c csList -o outFile [BEDFILE]...

 Description:
  Calculate global centrality metrics.

 Mandatory arguments:
  BEDFILE	Bed file(s). Expected to be ordered per condition.
  -c csList	Cutsite list file.
  -o outFile	Output matrix file.

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
	esac
done

# Check mandatory options
if [ -z "$csList" ]; then
	msg="!!! Missing mandatory -c option."
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
		echo -e " $helps<n$msg"
		exit 1
	fi
done

# RUN ==========================================================================

# Default empty matrix
matrix=""

# Cycle over chromosomes
for chr in $(echo $(seq 1 22) X Y); do
	chr=" chr$chr"
	#echo "Working on $chr..."

	# Number of cutsite in the chromosome
	ncs=`cat $csList | grep $chr | wc -l`

	if [ 0 -eq $ncs ]; then
		continue
	fi

	# Array of normalized counts
	counts=(0)

	for bf in ${bedfiles[@]}; do
		#echo "  BED file: $bf"

		# Number of reads in the condition
		ncc=`cat $bf | sed 1d | cut -f 5 | paste -sd+ | bc`
		
		if [ 0 -eq $ncc ]; then
			continue
		fi

		# Number of reads in the condition in the chromosome
		n=`cat $bf | sed 1d | grep $chr | cut -f 5 | paste -sd+ | bc`
		
		if [ -z "$n" ]; then
			n=0
		fi
		
		# Normalize
		r=`bc -l <<< "$n / $ncc / $ncs"`

		# Save normalized counts
		counts+=(`bc <<< "${counts[${#counts[@]} - 1]}+$r"`)
	done

	# Remove starting value (0)
	unset counts[0]

	# Merging cumulative set
	merged=`join_by " " "${counts[@]}"`

	# Forming new matrix row
	newrow=`echo "$chr $merged" | sed "s/^ //" | tr -s " " | tr " " "\t"`

	# Adding row to matrix
	matrix="$matrix$newrow\n"
done

# Writing matrix
echo -e "$matrix" > matrix.tmp.dat




# build rankings by comparing consecutive conditions
# 	compare with ratio B/A
# Make table with position and rankings
# Sum rankings and sort

# End --------------------------------------------------------------------------

################################################################################
