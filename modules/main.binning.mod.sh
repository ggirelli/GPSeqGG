#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 1.0.0
# Description: module for binning.
# 
# ------------------------------------------------------------------------------



# START ========================================================================

function bin_step() {
	echo -e 'Binning\n=====================\n'

	# Get binSize by input
	input_int 1e6 'Binsize' $binSize
	binSize=$v

	# Get binStep by input
	input_int 1e6 'Binstep' $binStep
	binStep=$v

	# Get cutsite list file by input
	input_fname 'Cutsite list file' 'a list of cutsite positions' $csList
	cutsitelist=$v

	# Get maskfile by input
	input_fname 'Mask file' 'a list of regions to be masked' $maskFile
	maskfile=$v

	# Get chrlengths by input
	input_fname 'Chr length file' \
		'lengths of chromosomes in the specified genome version' $chrLengths
	chrlengths=$v

	cslbool=0
	if [[ -n $csList ]]; then
		cslbool=1
	fi

	if [ -e "$csList" ]; then
		# Bin cutsites -------------------------------------------------------------
		echo -e "\nBinning cutsites ..."
		time $scriptdir/cs_bin.R -i $binSize -t $binStep -c $threads \
			$xout $csList $chrLengths & pid0=$!
		wait $pid0
	else
		echo -e "\nSkipped cutsite binning ..."
	fi

	if [ -e "$chrLengths" ]; then
		# Bin UMIs -----------------------------------------------------------------
		echo -e "\nBinning UMIs ..."
		if [ -z "$csList" ]; then
			$scriptdir/umi_bin.R -i $binSize -t $binStep -p $threads \
				$out/ $expID $conds $chrLengths & pid0=$!
		else
			$scriptdir/umi_bin.R -c -i $binSize -t $binStep -p $threads \
				$out/ $expID $conds $chrLengths & pid0=$!
		fi
		wait $pid0
	else
		echo -e "\nSkipped UMI binning ..."
	fi
}

################################################################################
