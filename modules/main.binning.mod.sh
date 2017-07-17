#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 1.1.0
# Description: module for binning.
# 
# ------------------------------------------------------------------------------



# START ========================================================================

function bin_step() {
	echo -e 'Binning\n=====================\n'

	cslbool=0
	if [[ -n $csList ]]; then
		cslbool=1
	fi

	if [ -e "$csList" ]; then
		# Bin cutsites -------------------------------------------------------------
		echo -e "Binning cutsites ..."
		time $scriptdir/cs_bin.R -i $binSize -t $binStep -c $threads \
			$xout $csList $chrLengths & pid0=$!
		wait $pid0
	else
		echo -e "Skipped cutsite binning ...\n"
	fi

	if [ -e "$chrLengths" ]; then
		conds=`grep -P "^$expID" $indir/patterns.tsv | cut -f 2 \
			| tr '\n' ',' | sed -r 's/^(.*),$/\1/'`

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
		echo -e "Skipped UMI binning ..."
	fi
}

################################################################################
