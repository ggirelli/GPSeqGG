#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 0.1.1
# Description: module for UMI preparation.
# 
# ------------------------------------------------------------------------------



# START ========================================================================

function prepare_umi() {
	echo -e 'UMI grouping & deduplicating\n=====================\n'

	# Get cutsite list file by input
	input_fname 'Cutsite list file' 'a list of cutsite positions' $csList
	cutsitelist=$v

	# Get maskfile by input
	input_fname 'Mask file' 'a list of regions to be masked' $maskFile
	maskfile=$v

	function prepare_umi_single_condition() {
		echo -e "\nPreparing UMIs from condition '$condition'..."
		cslbool=0
		if [[ -n $csList ]]; then
			cslbool=1
		fi

		# Group UMIs -----------------------------------------------------------
		if [[ -n $maskFile ]]; then
			$scriptdir/umi_group.py $cout/$condition/ $expID $condition \
				$umiLength --mask-file $maskFile & pid0=$!
			wait $pid0
		else
			$scriptdir/umi_group.py $cout/$condition/ $expID $condition \
				$umiLength & pid0=$!
			wait $pid0
		fi

		if [[ -n $csList ]]; then
			# Assign UMIs to cutsites ------------------------------------------
			$scriptdir/pos2cutsites.R $cout/$condition/ $expID $condition \
				$csList -i $csRange -c $threads & pid0=$!
			wait $pid0
		fi

		if [[ $umiLength -ne 0 ]]; then
			# Deduplicating UMIs -----------------------------------------------
			echo -e "\nDeduplicating UMIs ..."
			$scriptdir/umi_dedupl.R $cout/$condition/ $expID $condition \
				-p $platform -co $pthr -c $threads -cs $cslbool \
				-em $emax -ep $eperc & pid0=$!
			wait $pid0
		else
			if [[ -n $csList ]]; then
				cp $cout/$condition/UMIpos.atcs.txt \
					$cout/$condition/UMIpos.unique.atcs.txt
			else
				cp $cout/$condition/UMIpos.txt \
					$cout/$condition/UMIpos.unique.txt
			fi
		fi
	}
	for condition in "${condv[@]}"; do
		time prepare_umi_single_condition
	done
}

################################################################################
