#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 1.0.0
# Description: module for pattern filtering.
# 
# ------------------------------------------------------------------------------



# START ========================================================================

function pattern_filtering() {
	echo -e 'Pattern filtering\n=====================\n'

	# Count total reads
	echo -e "Counting total reads..."
	total_read_count=`wc -l $in/r1oneline.fa | cut -d " " -f 1`
	header="condition\tversion\tpattern\ttrim_len\ttotal_read_count"
	header="$header\treads_with_prefix\tprefix/total"
	echo -e $header > $out/summary

	# Work on single conditions
	patfiles="$indir/patterns.tsv"
	for condition in ${condv[@]}; do
		echo -e "\n> Working on condition '$condition'..."

		# Select pattern
		needle="^$expID\t$condition\t"
		pattern=`grep -P $needle $patfiles | cut -f 3`
		tlen=`grep -P "$needle" $patfiles | cut -f 4`

		# Save condition-specific patfile
		patfile="$cout/$condition/pat_file"
		mkdir -p $cout/$condition
		grep -P "^$expID\t$condition\t" $patfiles \
			> "$cout/$condition/patterns.tsv"
		echo "$pattern" > $patfile
		echo -e "Pattern: $pattern\n"

		# Identify condition-specific reads
		time $scriptdir/pattern_filter.sh -t $threads -i $in \
			-o "$cout/$condition" -p $patfile & pid0=$!
		wait $pid0

		# Print condition-specific read count in the summary
		count=`cat $cout/"$condition"/filtered.r1.fa | paste - - | wc -l`

		convstr="scale=4;perc=$count/$total_read_count;scale=2;(perc*100)/1"
		perc=`echo $convstr | bc`

		echo -e $header > $cout/"$condition"/summary
		header="$condition\t$version\t$pattern\t$tlen"
		header="$header\t$total_read_count\t$count\t$perc%"
		echo -e $header >> $cout/"$condition"/summary
		echo -e $header >> $out/summary
	done

	cp $out/summary $outcontrol/summary_pattern
}

################################################################################
