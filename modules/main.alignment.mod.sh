#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 0.1.1
# Description: module for read alignment.
# 
# ------------------------------------------------------------------------------



# START ========================================================================

function alignment() {
	echo -e 'Alignment\n=====================\n'

	new_fields="\tmapped\tmapped/prefix\tproperly_paired"
	new_fields="$new_fields\tproperly_paired/prefix\tinter_chr"
	head -n 1 $outcontrol/summary_pattern | \
		awk -v nf="$new_fields" "{ print \$0 nf }" - \
		> $outcontrol/summary_align

	patfiles="$indir/pat_files"
	for condition in "${condv[@]}"; do
		# echo -e "Aligning reads from condition '$condition'..."

		# Run trimmer ----------------------------------------------------------
		$scriptdir/reads_trim.sh -o $cout -c "$condition" -p $patfiles

		# Run aligner ----------------------------------------------------------
		if [ $numb_of_files -eq 2 ]; then
			# Paired-end
			$scriptdir/reads_align.sh -t $threads -o "$cout/$condition" \
				-c "$condition" -p -r $refGenome -a $aligner
		else
			# Single-end
			$scriptdir/reads_align.sh -t $threads -o "$cout/$condition" \
				-c "$condition" -r $refGenome -a $aligner
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

		split(a[$1], t, OFS);

		t[10]=$2 "\t" t[10];
		t[11]=$4 "\t" t[11];

		print join(t, 1, length(t)-1, OFS) }'
		awk "$awkprogram" \
			<(cat "$cout/$condition/$condition.sam" | \
				grep -v "^\@" | tr -s ' ') \
			<(cat "$cout/$condition/filtered.r1.linkers.oneline.fq") \
			> "$cout/$condition/$condition.linkers.sam"
	done

	cp $outcontrol/summary_align $out/summary
}

################################################################################
