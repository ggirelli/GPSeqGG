#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 0.1.1
# Description: module for plotting.
# 
# ------------------------------------------------------------------------------



# START ========================================================================

function plot_umi() {
	echo -e 'Plotting\n=====================\n'

	# Get cutsite list file by input
	input_fname 'Cutsite list file' 'a list of cutsite positions' $csList
	cutsitelist=$v

	cslbool=0
	if [[ -n $csList ]]; then
		cslbool=1
	fi

	# Get binSize by input
	input_int 1e6 'Binsize' $binSize
	binSize=$v

	# Get binStep by input
	input_int 1e6 'Binstep' $binStep
	binStep=$v

	# Get maskfile by input
	input_fname 'Mask file' 'a list of regions to be masked' $maskFile
	maskFile=$v
	if [ -z "$maskFile" ]; then
		msg="A list of masked regions is required for the plot step."
		msg="$msg\nExit."
		echo -e "$msg"
		exit 1
	fi

	# Get chrlengths by input
	input_fname 'Chr length file' \
		'lengths of chromosomes in the specified genome version' $chrLengths
	chrLengths=$v
	if [ -z "$chrLengths" ]; then
		msg="A list of chromosome lengths is required for the plot step."
		msg="$msg\nExit."
		echo -e "$msg"
		exit 1
	fi

	# Make multi-condition plots -----------------------------------------------
	
	# Prepare flags for heterochromosomes removal
	flags=""
	if [ "$rmX" = true ]; then flags="$flags --rmChrX"; fi
	if [ "$rmY" = true ]; then flags="$flags --rmChrY"; fi
	if [[ -n $neg ]]; then flags="$flags --neg $neg"; fi

	# Run plotting script
	$scriptdir/umi_plot.R $flags -i $binSize -t $binStep -c $threads \
		$out/ $expID $conds $csList $chrLengths $maskFile & pid0=$!
	wait $pid0
}

################################################################################
