#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 1.0.0
# Description: module for read alignment.
# 
# ------------------------------------------------------------------------------



# START ========================================================================

function alignment() {
	echo -e 'Alignment\n====================='

	# Update summary header
	header=`head -n 1 $outcontrol/summary_pattern`
	header="$header\tmapped\tmapped/prefix\tproperly_paired"
	header="$header\tproperly_paired/prefix\tinter_chr"
	echo -e "$header" > $outcontrol/summary_align

	patfiles="$indir/patterns.tsv"
	for condition in ${condv[@]}; do
		echo -e "\nAligning reads from condition '$condition'..."

		# Run trimmer ----------------------------------------------------------
		time $scriptdir/reads_trim.sh -o $cout \
			-e $expID -c "$condition" -p $patfiles

		# Run aligner ----------------------------------------------------------
		if [ -n "$bwaIndex" -a "bwa" == "$aligner" ]; then
			bwaOpt="-i $bwaIndex"
		fi
		if [ $numb_of_files -eq 2 ]; then
			# Paired-end
			time $scriptdir/reads_align.sh -t $threads -o "$cout/$condition" \
				-c "$condition" -p -r $refGenome -a $aligner $bwaOpt
		else
			# Single-end
			time $scriptdir/reads_align.sh -t $threads -o "$cout/$condition" \
				-c "$condition" -r $refGenome -a $aligner $bwaOpt
		fi

		# Update summary -------------------------------------------------------
		echo -e " · Retrieving flagstats ..."
		samtools flagstat $cout/"$condition"/"$condition".sorted.bam \
			> $cout/"$condition"/"$condition".bam_notes.txt

		# Mapped reads
		count=`cat $cout/"$condition"/"$condition".bam_notes.txt | \
			grep 'mapped' | head -n 1 | cut -d ' ' -f 1`
		perc=`cat $cout/"$condition"/"$condition".bam_notes.txt | \
			grep 'mapped' | head -n 1 | cut -d '(' -f 2 | cut -d ':' -f 1`

		# Properly paired reads
		count2=`cat $cout/"$condition"/"$condition".bam_notes.txt | \
			grep 'properly paired' | head -n 1 | cut -d ' ' -f 1`
		perc2=`cat $cout/"$condition"/"$condition".bam_notes.txt | \
			grep 'properly paired' | head -n 1 | cut -d '(' -f 2 | \
			cut -d ':' -f 1`
		
		# Chimeric reads	
		count3=`cat $cout/"$condition"/"$condition".bam_notes.txt | \
			grep 'with mate mapped to a different chr' | head -n 1 | \
			cut -d ' ' -f 1`

		# Add to summary
		new_fields="\t$count\t$perc\t$count2\t$perc2\t$count3"
		grep "$condition" "$outcontrol/summary_pattern" | \
			awk -v nf="$new_fields" "{ print \$0 nf }" - \
			>> $outcontrol/summary_align

		# Add back the UMIs to the SAM file ------------------------------------
		echo -e " · Adding linkers to SAM file ..."

		awkprogram='@include "join";
		FNR==NR{a[$1]=$0;next} ($1 in a) { OFS="\t";

		split(a[$1], t1, OFS);
		split($0, t2, OFS);

		t2[10]=t1[2] "\t" t2[10];
		t2[11]=t1[4] "\t" t2[11];

		print join(t2, 1, 11, OFS) }'
		awk "$awkprogram" \
			<(cat "$cout/$condition/filtered.r1.linkers.oneline.fq") \
			<(cat "$cout/$condition/$condition.sam" | \
				grep -v "^\@" | tr -s ' ') \
			> "$cout/$condition/$condition.linkers.sam"

		# Clean ----------------------------------------------------------------
		
		if [ 0 -lt $neatness ]; then
			echo -e "\n~ Cleaning..."
			rm -v $cout/$condition/filtered*
		fi

	done

	cp $outcontrol/summary_align $out/summary
}

################################################################################
