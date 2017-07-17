#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 1.0.1
# Description: module for UMI preparation.
# 
# ------------------------------------------------------------------------------



# START ========================================================================

function prepare_umi() {
	echo -e 'UMI grouping & deduplicating\n=====================\n'

	# Update summary header
	header=`head -n 1 $outcontrol/summary_sam_filter`
	header="$header\tunmasked\tunmasked/umis"
	header="$header\tnon_orphan\tnon_orph/unmasked"
	header="$header\tpass_read_qc\tread_qc/non_orph"
	header="$header\tunique_umis\tunique/read_qc"
	echo -e "$header" > $outcontrol/summary_umi_prep

	function prepare_umi_single_condition() {
		echo -e "\nPreparing UMIs from condition '$condition'..."

		cslbool=0
		if [[ -n "$csList" ]]; then
			cslbool=1
		fi

		# Group UMIs -----------------------------------------------------------
		if [[ -n "$maskFile" ]]; then
			time $scriptdir/umi_group.py $cout/$condition/ $expID $condition \
				$umiLength --mask-file $maskFile & pid0=$!
			wait $pid0
		else
			time $scriptdir/umi_group.py $cout/$condition/ $expID $condition \
				$umiLength & pid0=$!
			wait $pid0
		fi

		if [[ -n "$csList" ]]; then
			# Assign UMIs to cutsites ------------------------------------------
			time $scriptdir/pos2cutsites.R $cout/$condition/ $expID $condition \
				$csList -i $csRange -c $threads & pid0=$!
			wait $pid0
		fi

		if [[ $umiLength -ne 0 ]]; then
			# Deduplicating UMIs -----------------------------------------------
			echo -e "\nDeduplicating UMIs ..."
			time $scriptdir/umi_dedupl.R $cout/$condition/ $expID $condition \
				-p $platform -co $pthr -c $threads -cs $cslbool \
				-em $emax -ep $eperc & pid0=$!
			wait $pid0
		else
			if [[ -n "$csList" ]]; then
				cp $cout/$condition/UMIpos.atcs.txt \
					$cout/$condition/UMIpos.unique.atcs.txt
			else
				cp $cout/$condition/UMIpos.txt \
					$cout/$condition/UMIpos.unique.txt
			fi
		fi

		# Generate bed file ----------------------------------------------------

		# Standard bed header
		trackName="track name=\"$expID.$condition.dedupUMIs\" "
		trackName="$trackName description=\"emax=$emax,eperc=$eperc"
		trackName="$trackName,csRange=$csRange,MAPQthr=$mapqThr,pthr=$pthr"

		if [ -n "$csList" ]; then
			# Update bed header
			trackName="$trackName,atcs=T\""

			# Generate bed
			$scriptdir/bed_make.sh -o $cout/$condition/ -c $cutsite \
				-f $csList -t "$trackName"

			# Save bed
			cp $cout/$condition/UMIpos.unique.atcs.bed \
				$out/$expID"_"$condition"_GG__cutsiteLoc-umiCount.bed"
		else
			# Update bed header
			trackName="$trackName,atcs=F\""

			# Generate bed
			$scriptdir/bed_make.sh -o $cout/$condition/ -t "$trackName"

			# Save bed
			cp $cout/$condition/UMIpos.unique.bed \
				$out/$expID"_"$condition"_GG__cutsiteLoc-umiCount.bed"
		fi

		# Update summary -------------------------------------------------------

		# Input reads
		ir=`cat $cout/"$condition"/"$condition".umi_prep_notes.txt | \
			grep 'input reads' | head -n 1 | cut -d ' ' -f 1`

		# Masked reads
		c=`cat $cout/"$condition"/"$condition".umi_prep_notes.txt | \
			grep 'reads masked' | head -n 1 | cut -d ' ' -f 1`
		if [ -z "$c" ]; then c=$ir; else c=`expr $ir - $c`; fi
		p=`printf "%.2f%%" "$(bc <<< "scale = 4; $c / $ir * 100")"`

		# Assigned to cutsite (non-orphan)
		c2=`cat $cout/"$condition"/"$condition".umi_prep_notes.txt | \
			grep 'non-orphan' | head -n 1 | cut -d ' ' -f 1`
		p2=`printf "%.2f%%" "$(bc <<< "scale = 4; $c2 / $c * 100")"`
		
		# Pass the read quality filter	
		c3=`cat $cout/"$condition"/"$condition".umi_prep_notes.txt | \
			grep 'pass the read quality filter' | head -n 1 | cut -d ' ' -f 1`
		p3=`printf "%.2f%%" "$(bc <<< "scale = 4; $c3 / $c2 * 100")"`

		# Unique UMIs
		c4=`cat $cout/"$condition"/"$condition".umi_prep_notes.txt | \
			grep 'UMIs left after deduplication' | head -n 1 | cut -d ' ' -f 1`
		p4=`printf "%.2f%%" "$(bc <<< "scale = 4; $c4 / $c3 * 100")"`

		# Add to summary
		new_fields="\t$c\t$p\t$c2\t$p2\t$c3\t$p3\t$c4\t$p4"
		grep "$condition" "$outcontrol/summary_sam_filter" | \
			awk -v nf="$new_fields" "{ print \$0 nf }" - \
			>> $outcontrol/summary_umi_prep


		# Clean ----------------------------------------------------------------
		
		if [ 2 -le $neatness ]; then
			echo -e "\n~ Cleaning..."
			rm -v $cout/$condition/$condition* \
				$cout/$condition/*.log \
				$cout/$condition/summary \
				$cout/$condition/UMIpos.*
		fi
	}
	for condition in ${condv[@]}; do
		time prepare_umi_single_condition
	done

	cp $outcontrol/summary_umi_prep $out/summary
}

################################################################################
