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

	# Update summary header
	new_fields="\tsecondary_aln\tchimeras\tunmapped\tr2\tmapq < $mapqthr\trmChr"
	new_fields="$new_fields\tafterSAMfilter"
	head -n 1 $outcontrol/summary_pattern | \
		awk -v nf="$new_fields" "{ print \$0 nf }" - \
		> $outcontrol/summary_sam_filter

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

		# Update summary -------------------------------------------------------

		# Secondary alignments
		c=`cat $cout/"$condition"/"$condition".sam_filter_notes.txt | \
			grep 'secondary alignments' | head -n 1 | cut -d ' ' -f 1`

		# Chimeric reads
		c2=`cat $cout/"$condition"/"$condition".sam_filter_notes.txt | \
			grep 'chimeric reads' | head -n 1 | cut -d ' ' -f 1`
		
		# Unmapped reads	
		c3=`cat $cout/"$condition"/"$condition".sam_filter_notes.txt | \
			grep 'unmapped reads' | head -n 1 | cut -d ' ' -f 1`

		# R2 reads	
		c4=`cat $cout/"$condition"/"$condition".sam_filter_notes.txt | \
			grep 'R2 reads' | head -n 1 | cut -d ' ' -f 1`

		# MAPQ < thr reads	
		c5=`cat $cout/"$condition"/"$condition".sam_filter_notes.txt | \
			grep 'reads with MAPQ <' | head -n 1 | cut -d ' ' -f 1`

		# Chr removed reads	
		c6=`cat $cout/"$condition"/"$condition".sam_filter_notes.txt | \
			grep 'reads from removed chromosomes' | head -n 1 | cut -d ' ' -f 1`

		# Surviving reads	
		c7=`cat $cout/"$condition"/"$condition".sam_filter_notes.txt | \
			grep 'reads left after filtering' | head -n 1 | cut -d ' ' -f 1`

		# Add to summary
		new_fields="\t$c\t$c2\t$c3\t$c4\t$c5\t$c6\t$c7"
		grep "$condition" "$outcontrol/summary_pattern" | \
			awk -v nf="$new_fields" "{ print \$0 nf }" - \
			>> $outcontrol/summary_sam_filter
	done

}

################################################################################
