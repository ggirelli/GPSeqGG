#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 1.1.0
# Description: module for plotting.
# 
# ------------------------------------------------------------------------------



# START ========================================================================

function plot_umi() {
	echo -e 'Plotting\n=====================\n'
	can_plot=true

	cslbool=0
	if [[ -n $csList ]]; then
		cslbool=1
	fi

	if [ -e "$chrLengths" -a -e "$csList" ]; then
		# Make multi-condition plots -------------------------------------------
		conds=`grep -P "^$expID" $indir/patterns.tsv | cut -f 2 \
			| tr '\n' ',' | sed -r 's/^(.*),$/\1/'`

		# Prepare flags for heterochromosomes removal
		flags=""
		if [ "$rmX" = true ]; then flags=$flags" --rmChrX"; fi
		if [ "$rmY" = true ]; then flags=$flags" --rmChrY"; fi
		if [[ -n $neg ]]; then flags=$flags" --neg $neg"; fi
		if [ -e "$maskFile" ]; then flags=$flags" -m "$maskFile; fi

		# Run plotting script
		$scriptdir/umi_plot.R $flags -i $binSize -t $binStep -c $threads \
			$out/ $expID $conds $csList $chrLengths & pid0=$!
		wait $pid0
	else
		echo -e "Skipped plot step ..."
	fi
}

################################################################################
