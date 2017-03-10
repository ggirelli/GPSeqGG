#!/usr/bin/env bash

# THIS SCRIPT CAN BE CALLED AS
# ./find_seq.sh expName conditions
# The seqs file MUST end with an empty line
################################################################################
clear

if [ "$#" -ne 2 ]; then
	echo -e "Correct usage:\n./find_seq.sh expName conditions"
	exit 1
fi

# DEFINING VARIABLES
experiment=$1		# e.i. TK20
conds=$2
numbproc=8

################################################################################
# PREPARE DIRECTORY STRUCTURE

datadir=$DATA/BiCro-Data/Sequencing
indir=$datadir/$experiment
in=$datadir/$experiment/indata

resdir=$DATA/BiCro-Analysis/Sequencing
outdir=$resdir/$experiment
out=$outdir/outdata
outcontrol=$outdir/outdata.control

aux=$outdir/auxdata
auxcontrol=$outdir/auxdata.control

refgen=$DATA/BiCro-Resources/genomes/$genome/list.fa

scriptdir=$REPO/GPSeq-Sequencing

################################################################################
# LOAD DATA FILES

find $indir -maxdepth 1 -type f -iname "*$experiment*" | sort > filelist
numb_of_files=`cat filelist | wc -l`
r1=`cat filelist | head -n1`
if [ $numb_of_files -eq 2 ]; then
    r2=`cat filelist | tail -n1`
fi
rm filelist

################################################################################

# Initialize summaries
seqs="$indir/seqs"
echo "" > $out/summary_seqfinder
while read s; do
	echo -e "#seq$counter\t$s" >> $out/summary_seqfinder
done < $seqs
echo -e "\ncondition\tseq\tcount\tperc" >> $out/summary_seqfinder
cp $out/summary_seqfinder $out/summary_seqfinder_conditions

# Count total reads
total_read1_count=`wc -l $in/r1oneline.fq | cut -d " " -f 1`
total_read2_count=`wc -l $in/r2oneline.fq | cut -d " " -f 1`

# Search seq in unfiltered files
counter=1
while read s; do
	seqfile=$outcontrol/seq_tmp
	echo "$s" > $seqfile

	echo -e "Searching seq$counter..."

	# Search in R1
	echo -e "\tWorking on R1..."
	cat "$in/r1.fa" | parallel --tmpdir $HOME/tmp --block 100M -k --pipe -L 2 "scan_for_matches $seqfile - " > $outcontrol/seq"$counter".r1.fa & pid1=$!
	wait $pid1
	count=`wc -l $outcontrol/"seq$counter".r1.fa | cut -d " " -f 1`
	perc=`echo "scale=4; perc = $count / $total_read1_count; scale=2; (perc * 100)/1" | bc`
	echo -e "R1\tseq$counter\t$count\t$perc" >> $out/summary_seqfinder

	# Search in R2
	if [ $numb_of_files -eq 2 ]; then
		echo -e "\tWorking on R2..."
		cat $in/r2.fa | parallel --tmpdir $HOME/tmp --block 100M -k --pipe -L 2 "scan_for_matches $seqfile - " > $outcontrol/"seq$counter".r2.fa & pid1=$!
		wait $pid1
		count=`wc -l $outcontrol/"seq$counter".r2.fa | cut -d " " -f 1`
		perc=`echo "scale=4; perc = $count / $total_read2_count; scale=2; (perc * 100)/1" | bc`
		echo -e "R2\tseq$counter\t$count\t$perc" >> $out/summary_seqfinder
	fi

	counter=`echo "$counter + 1" | bc`
	rm $seqfile
done < $seqs

# Work on single conditions
IFS=',' read -r -a condv <<< "$conds"
for condition in "${condv[@]}"; do
	echo -e "\nWorking on condition '$condition'..."
	mkdir -p $outcontrol/"$condition"

	# Search sequences
	counter=1
	while read s; do
		echo -e "\tSearching seq$counter..."

		seqfile=$outcontrol/"$condition"/seq_tmp
		echo "$s" > $seqfile

		cat $out/"$condition"/filtered.r1.fa | parallel --tmpdir $HOME/tmp --block 100M -k --pipe -L 2 "scan_for_matches $seqfile - " > $outcontrol/"$condition"/seq"$counter"_filtered.r1.fa & pid1=$!
		wait $pid1
		count=`wc -l $outcontrol/"$condition"/seq"$counter"_filtered.r1.fa | cut -d " " -f 1`
		perc=`echo "scale=4; perc = $count / $total_read1_count; scale=2; (perc * 100)/1" | bc`
		echo -e "$condition\tseq$counter\t$count\t$perc" >> $out/summary_seqfinder_conditions

		counter=`echo "$counter + 1" | bc`
		rm $seqfile
	done < $seqs
done