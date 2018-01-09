#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 1.1.0
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
	header=`head -n 1 $outcontrol/summary_pattern`
	header="$header\tsecondary_aln\tchimeras\tr2\tunmapped\tmapq<$mapqthr"
	header="$header\trmChr\tumis\tumis/prefix"
	echo -e "$header" > "$outcontrol/summary_sam_filter"

	for condition in ${condv[@]}; do
		echo -e "\nAnalyzing UMIs from condition '$condition'..."

		# Count condition total reads
		echo -e "Counting condition reads..."
		cond_count=`grep "^$condition" $outcontrol/summary_align | cut -f 5 \
			| head -n 1`

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
		if [ -n "$flags" ]; then
			flags="-C $flags"
		fi
		if [ -n "$cutsite" ]; then
			flags="$flags -l ${#cutsite}"
		fi

		# Additional SAM filters
		time $scriptdir/sam_filter.sh -p1cm -t $threads -q $mapqthr $flags \
			-i "$cout/$condition/$condition.linkers.sam" \
			2> "$cout/$condition/$condition.sam_filter_notes.txt" & pid0=$!
		wait $pid0

		# Extract linker information from SAM file
		awkprg='
		BEGIN { OFS = FS = "\t"; }
		{
			# Extract linker information
			match($0, "LS:Z:([^[:blank:]]*)", lseq)
			match($0, "LQ:Z:([^[:blank:]]*)", lqual)

			# Output
			print "chr" $3 OFS $4 OFS $5 OFS lseq[1] OFS lqual[1] OFS $1;
		}'
		cat "$cout/$condition/$condition.linkers.filtered.sam" | \
			awk "$awkprg" > "$cout/$condition/$condition.filtered.umi.pos.txt"

		# Update summary -------------------------------------------------------

		# Secondary alignments
		c=`cat $cout/"$condition"/"$condition".sam_filter_notes.txt | \
			grep 'secondary alignments' | head -n 1 | cut -d ' ' -f 1`

		# Chimeric reads
		c2=`cat $cout/"$condition"/"$condition".sam_filter_notes.txt | \
			grep 'chimeric reads' | head -n 1 | cut -d ' ' -f 1`

		# R2 reads	
		c3=`cat $cout/"$condition"/"$condition".sam_filter_notes.txt | \
			grep 'R2 reads' | head -n 1 | cut -d ' ' -f 1`

		# Unmapped reads	
		c4=`cat $cout/"$condition"/"$condition".sam_filter_notes.txt | \
			grep 'unmapped reads' | head -n 1 | cut -d ' ' -f 1`

		# MAPQ < thr reads	
		c5=`cat $cout/"$condition"/"$condition".sam_filter_notes.txt | \
			grep 'reads with MAPQ <' | head -n 1 | cut -d ' ' -f 1`

		# Chr removed reads	
		c6=`cat $cout/"$condition"/"$condition".sam_filter_notes.txt | \
			grep 'reads from removed chromosomes' | head -n 1 | cut -d ' ' -f 1`

		# Surviving reads	
		c7=`cat $cout/"$condition"/"$condition".sam_filter_notes.txt | \
			grep 'reads left after filtering' | head -n 1 | cut -d ' ' -f 1`
		p7=`printf "%.2f%%" "$(bc <<< "scale = 4; $c7 / $cond_count * 100")"`

		# Add to summary
		new_fields="\t$c\t$c2\t$c3\t$c4\t$c5\t$c6\t$c7\t$p7"
		grep "^$condition" "$outcontrol/summary_pattern" | \
			awk -v nf="$new_fields" "{ print \$0 nf }" - \
			>> "$outcontrol/summary_sam_filter"

		# Clean ----------------------------------------------------------------
		
		if [ 0 -lt $neatness ]; then
			echo -e "\n~ Cleaning..."
			rm -v "$cout/$condition/$condition.sam"
		fi

	done

	cp "$outcontrol/summary_sam_filter" "$out/summary"
}

################################################################################
