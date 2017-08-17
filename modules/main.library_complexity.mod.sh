#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 1.1.0
# Description: module for library complexity estimation.
# 
# ------------------------------------------------------------------------------



# START ========================================================================

function library_complexity() {

	echo -e 'Library complexity estimation\n=====================\n'

	# Check if preseq is executable
	preseq="$libdir/preseq/preseq"
	if [ -z "$(command -v $preseq)" ]; then
		echo -e "preseq must be compiled to estimate library complexity."
		return
	fi

	# Check if datamash is available
	if [ -z "$(command -v datamash)" ]; then
		echo -e "datamash must be installed to estimate library complexity."
		return
	fi

	# Remove merged files if they exist already
	kind=("c_curve" "lc_extrap" "bound_pop")
	for k in ${kind[@]}; do
		if [ -e $xout"/lc."$k".txt" ]; then rm $xout"/lc."$k".txt"; fi
	done

	# Run per condition
	patfiles="$indir/patterns.tsv"
	for condition in ${condv[@]}; do
		echo -e " · Working on condition '$condition' ..."

		# Identify proper input file
		prefix="UMIpos.unique"
		if [[ -n "$csList" ]]; then
			prefix=$prefix".atcs"
		fi

		max_count=$(cut -f 4 $cout/$condition/$prefix".txt" | tr ' ' '\n' \
			| datamash max 1)
		if [ 4 -gt $max_count ]; then
			msg="Skipping condition '$condition': "
			msg=$msg"not enough duplicates to estimate library complexity."
			echo -e "$msg"
			continue
		fi

		# Create UMI count file
		cut -f 4 $cout/$condition/$prefix".txt" | tr ' ' '\n' \
			> $cout/$condition/$prefix".umi_counts.tmp"

		# Run preseq
		$preseq c_curve -s 100 -o $cout/$condition/"lc."$prefix".c_curve.tmp" \
			-V $cout/$condition/$prefix".umi_counts.tmp"
		$preseq lc_extrap -o $cout/$condition/"lc."$prefix".lc_extrap.tmp" \
			-V $cout/$condition/$prefix".umi_counts.tmp"
		$preseq bound_pop -o $cout/$condition/"lc."$prefix".bound_pop.tmp" \
			-V $cout/$condition/$prefix".umi_counts.tmp"

		# Add condition column
		for k in ${kind[@]}; do
			header=$(head -n1 $cout/$condition/"lc."$prefix"."$k".tmp")
			header=$header"\tCONDITION"
			echo -e "$header" > $cout/$condition/"lc."$prefix"."$k".txt"
			cat $cout/$condition/"lc."$prefix"."$k".tmp" | sed 1d \
				| awk -v c="$condition" '{ print $0"\t"c }' \
				>> $cout/$condition/"lc."$prefix"."$k".txt"
			rm $cout/$condition/"lc."$prefix"."$k".tmp"
		done

		# Merging
		for k in ${kind[@]}; do
			infile=$cout/$condition/"lc."$prefix"."$k".txt"
			outfile=$xout"/lc."$k".txt"
			if [ ! -e $outfile ]; then
				head -n1 $infile > $outfile
			fi
			cat $infile | sed 1d >> $outfile
			rm $infile
		done
		
		# Remove tmp input
		rm "$cout/$condition/$prefix.umi_counts.tmp"
	done

	# Generate plots
	echo -e " · Generating plots."
	$scriptdir/library_complexity_plot.R $out $expID
}

################################################################################
