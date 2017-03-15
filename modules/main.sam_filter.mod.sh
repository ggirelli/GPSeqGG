#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 0.1.1
# Description: module for SAM filtering.
# 
# ------------------------------------------------------------------------------



# START ========================================================================

function filter_sam() {
	echo -e 'SAM filtering\n=====================\n'

	# Get mapq threshold by input
	input_int 30 'MAPQ threshold' $mapqThr
	mapqthr=$v

	for condition in "${condv[@]}"; do
		echo -e "\nAnalyzing UMIs from condition '$condition'..."

		# From https://goo.gl/MNXp5o
		function join_by { local IFS="$1"; shift; echo "$*"; }

		# Select chromosomes to remove
		flags=()
		if $rmX; then
			flags+=('X')
		fi
		if $rmY; then
			flags+=('Y')
		fi
		flags=`join_by , "${flags[@]}"`
		if [ -n $flags ]; then
			flags="-r $flags"
		fi

		# Filter SAM
		time $scriptdir/sam_filter.R $cout/$condition/ $expID $condition \
			-mt $mapqThr -cs $cutsite -c $threads $flags & pid0=$!
		wait $pid0

	done
}

################################################################################
