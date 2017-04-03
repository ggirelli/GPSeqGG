#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 1.0.1
# Description: 	merges bedfiles into a matrix.
# 				The score column is merged based on the position or the name.
# 
# ------------------------------------------------------------------------------



# ENV VAR ======================================================================

export LC_ALL=C

# DEPENDENCIES =================================================================

function join_by { local IFS="$1"; shift; echo "$*"; }

# INPUT ========================================================================

# Help string
helps="
 usage: ./bed2matrix.sh [-hn] -o outFile [BEDFILEs]...

 Description:
  Merge bedfiles into a matrix. The score column is merged based on the positon
  given by the chr+start+end columns (default) or by the name column (-n option)

 Mandatory arguments:
  BEDFILEs	Bed file(s). Expected to be ordered per condition.
  -o outFile	Output matrix file.

 Optional arguments:
  -h		Show this help page.
  -n		Merge bedfiles based on name instead of location.
"

# Default options
byName=false

# Parse options
while getopts hno: opt "${bedfiles[@]}"; do
	case $opt in
		h)
			echo -e "$helps\n"
			exit 0
		;;
		n)
			byName=true
		;;
		o)
			outFile=$OPTARG
		;;
	esac
done

# Check mandatory options
if [ -z "$outFile" ]; then
	msg="!!! Missing mandatory -o option."
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

# TEST =========================================================================

# RUN ==========================================================================

if $byName; then
	# Merge by name ------------------------------------------------------------
	echo -e " > Merging by name..."	
else
	# Merge by location --------------------------------------------------------
	echo -e " > Merging by location..."

	awkprogram='
	(FNR==NR){
		a[]
	}
	'
fi

# End --------------------------------------------------------------------------

################################################################################
