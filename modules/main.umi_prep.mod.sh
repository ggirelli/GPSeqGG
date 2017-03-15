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

	# Update summary header
	new_fields="\tunmasked\tunmasked/umis"
	new_fields="$new_fields\tassigned_to_cs\tto_cs/unmasked"
	new_fields="$new_fields\tpass_read_qc\tread_qc/to_cs"
	new_fields="$new_fields\tunique_umis\tunique/read_qc"
	head -n 1 $outcontrol/summary_sam_filter | \
		awk -v nf="$new_fields" "{ print \$0 nf }" - \
		> $outcontrol/summary_umi_prep

	function prepare_umi_single_condition() {
		echo -e "\nPreparing UMIs from condition '$condition'..."
		cslbool=0
		if [[ -n $csList ]]; then
			cslbool=1
		fi

		# Group UMIs -----------------------------------------------------------
		if [[ -n $maskFile ]]; then
			time $scriptdir/umi_group.py $cout/$condition/ $expID $condition \
				$umiLength --mask-file $maskFile & pid0=$!
			wait $pid0
		else
			time $scriptdir/umi_group.py $cout/$condition/ $expID $condition \
				$umiLength & pid0=$!
			wait $pid0
		fi

		if [[ -n $csList ]]; then
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
			if [[ -n $csList ]]; then
				cp $cout/$condition/UMIpos.atcs.txt \
					$cout/$condition/UMIpos.unique.atcs.txt
			else
				cp $cout/$condition/UMIpos.txt \
					$cout/$condition/UMIpos.unique.txt
			fi
		fi

		# Generate bed file ----------------------------------------------------

		cslength=${#cutsite}
		if [ -n $csList ]; then
			# Obtain cutsites from list and compare with uniqued UMIs
			awkprogram='
			FNR==NR
			{
				FS=OFS="\t";
				k=$1"~"$2;
				a[k]=$0"\t"NR;
				next
			}

			{
				FS=OFS="\t";
				k="chr"$1"~"$2;
			}
			(k in a){
				split(a[k], cs, "\t");
				split($0, umi, "\t");
				n=split(umi[3], umis, " ");
				{
				FS=OFS="";
				print cs[1],"\t",cs[2],"\t",cs[2]+cslen-1,"\tcs_",cs[3],"\t",n;
				}
			}'
			awk -v cslen=$cslength "$awkprogram" \
				<(cat "$csList") \
				<(cat "$cout/$condition/UMIpos.unique.atcs.txt") \
				> "$cout/$condition/UMIpos.unique.atcs.bed"

			# Copy to main directory
			cp "$cout/$condition/UMIpos.unique.atcs.bed" \
				> "$out"/"$expID"_"$condition"_GG__cutsiteLoc-umiCount.atcs.bed
		else
			# Without cutsite assignment
			awkprogram='
			{
				split($0, r, "\t");
				n=split(r[3], u, " ");
			}

			23==r[1]{ r[1]="X" }
			24==r[1]{ r[1]="Y" }

			0!=n{
				FS=OFS="";
				r[1]="chr"r[1];
				print r[1],"\t",r[2],"\t",r[2]+cslen-1,"\tloc_",NR,"\t",n
			}
			'
			awk -v cslen=$cslength "$awkprogram" \
				<(cat "$cout/$condition/UMIpos.unique.txt") \
				> "$cout/$condition/UMIpos.unique.bed"

			# Copy to main directory
			cp "$cout/$condition/UMIpos.unique.bed" \
				> "$out"/"$expID"_"$condition"_GG__cutsiteLoc-umiCount.bed
		fi

		# Update summary -------------------------------------------------------

		# Input reads
		ir=`cat $cout/"$condition"/"$condition".umi_prep_notes.txt | \
			grep 'input reads' | head -n 1 | cut -d ' ' -f 1`

		# Masked reads
		c=`cat $cout/"$condition"/"$condition".umi_prep_notes.txt | \
			grep 'reads masked' | head -n 1 | cut -d ' ' -f 1`
		c=`expr $ir - $c`
		p=`printf "%.2f%%" "$(bc <<< "scale = 4; $c / $ir * 100")"`

		# Assigned to cutsite
		c2=`cat $cout/"$condition"/"$condition".umi_prep_notes.txt | \
			grep 'assigned to a cutsite' | head -n 1 | cut -d ' ' -f 1`
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
	}
	for condition in "${condv[@]}"; do
		time prepare_umi_single_condition
	done

	cp $outcontrol/summary_umi_prep $out/summary
}

################################################################################
