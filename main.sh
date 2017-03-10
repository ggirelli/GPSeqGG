#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 0.1.1
# Description: analyze GPSeq sequencing data
# 
# Help page: ./main.sh -h
# 
# ------------------------------------------------------------------------------




# ENV VAR ======================================================================

export LC_ALL=C

# DEPENDENCIES =================================================================

# Script folder
scriptdir="`dirname ${BASH_SOURCE}`/"

# Load functions
source $scriptdir/main.functions.sh

# CHECK OPTIONS ================================================================

# Run HELP ---------------------------------------------------------------------

opt=$1
if [[ ${opt:0:1} == '-' ]]; then
	if [ $opt == '-h' ]; then
		cat 'docs/main.help'
		exit 1
	fi
fi

# Check params -----------------------------------------------------------------

if [ "$#" -lt 1 ]; then
	echo -e "Correct usage:\n./main [-opt] settingsFile"
	exit 1
fi

# Read settings file -----------------------------------------------------------

opt=$1
if [[ ${opt:0:1} == '-' ]]; then
	if [ "$#" -lt 2 ]; then
		echo -e "Correct usage:\n./main [-opt] settingsFile"
		exit 1
	fi
	settingsFile=$2
else
	settingsFile=$1
fi
. $settingsFile

# Check settings ---------------------------------------------------------------

# Print settings
echo -e '# SETTINGS ========================================================'
cat $settingsFile

# Check settings
check_settings # from main.functions

# PREPARE DIRECTORY STRUCTURE ==================================================

# Input folders
datadir=$DATA/BiCro-Data/Sequencing
indir=$datadir/$experiment && mkdir -p $indir
in=$datadir/$experiment/indata && mkdir -p $in

# Output folders
resdir=$DATA/BiCro-Analysis/Sequencing
outdir=$resdir/$experiment
out=$outdir/ && mkdir -p $out

# Additional outputs
outcontrol=$outdir/tmp && mkdir -p $outcontrol
logpath="$out/log"

# Genome reference
refgen=$DATA/BiCro-Resources/genomes/$genome/list.fa

# START LOG ====================================================================

{
echo -e '\n'
echo -e '# START ===========================================================\n'

# LOAD DATA FILES ==============================================================

find $indir -maxdepth 1 -type f -iname "*$experiment*R[12]*" | sort > filelist
numb_of_files=`cat filelist | wc -l`
r1=`cat filelist | head -n1`
echo "R1 is " $r1
if [ $numb_of_files -eq 2 ]; then
    r2=`cat filelist | tail -n1`
    echo "R2 is " $r2
fi
rm filelist

# START PIPELINE ===============================================================

# QUALITY CONTROL --------------------------------------------------------------
function quality_control() {
	echo -e '# Quality control ==============================================\n'
	# Produce quality control summarie(s)
	time $scriptdir/quality_control.sh -t $numbproc -o $out -1 $r1 -2 $r2
}
execute_step $opt 'quality control' quality_control

# FILE GENERATION --------------------------------------------------------------
function file_generation() {
	echo -e '# File generation ==============================================\n'
	# Generate necessary files
	time $scriptdir/files_prepare.sh -t $numbproc -o $in -1 $r1 -2 $r2
}
execute_step $opt 'file generation' file_generation

# PATTERN FILTERING ------------------------------------------------------------
function pattern_filtering() {
	echo -e '# Pattern filtering ============================================\n'

	# Count total reads
	echo -e "\nCounting total reads..."
	total_read_count=`wc -l $in/r1oneline.fa | cut -d " " -f 1`
	header="condition\tpattern\ttotal_read_count\treads_with_prefix"
	header="$header\tprefix/total"
	echo -e $header > $out/summary

	# Work on single conditions
	patfiles="$indir/pat_files"
	IFS=',' read -r -a condv <<< "$conds"
	for condition in "${condv[@]}"; do
		echo -e "\nWorking on condition '$condition'..."

		# Select pattern
		pattern=`grep -P "$condition\t" $patfiles | cut -f 2`

		# Save condition-specific patfile
		patfile="$out/$condition/pat_file"
		mkdir -p $out/$condition
		echo "$pattern" > $patfile
		echo -e "\nPattern: $pattern"

		# Identify condition-specific reads
		time $scriptdir/pattern_filter.sh -t $numbproc -i $in \
			-o "$out"/$condition -p $patfile & pid0=$!
		wait $pid0

		# Print condition-specific read count in the summary
		count=`cat $out/"$condition"/filtered.r1.fa | paste - - | wc -l`

		convstr="scale=4;perc=$count/$total_read_count;scale=2;(perc*100)/1"
		echo $constr
		perc=`echo $convstr | bc`

		echo -e $header > $out/"$condition"/summary
		header="$condition\t$pattern\t$total_read_count\t$count\t$perc%"
		echo -e $header >> $out/"$condition"/summary
		echo -e $header >> $out/summary
	done

}
execute_step $opt 'pattern_filtering' pattern_filtering

# ALIGNMENT --------------------------------------------------------------------
function alignment() {
	echo -e '# Alignment ====================================================\n'

	header=`head -n 1 $out/summary`
	echo -e "$header\tmapped/prefix\tproperly_paired\tinter_chr" \
		> $outcontrol/summarytmp

	IFS=',' read -r -a condv <<< "$conds"
	patfiles="$indir/pat_files"
	for condition in "${condv[@]}"; do
		echo -e "\nAligning reads from condition '$condition'..."

		# Run trimmer ----------------------------------------------------------
		$scriptdir/reads_trim.sh \
			-t $numbproc -o $out $c "$condition" -p $patfiles

		# Run aligner ----------------------------------------------------------
		if [ $numb_of_files -eq 2 ]; then
			# Paired-end
			$scriptdir/reads_align.sh -t $numbproc -o $out -c "$condition" \
				-p -r $genome -a $aligner
		else
			# Single-end
			$scriptdir/reads_align.sh -t $numbproc -o $out -c "$condition" \
				-r $genome -a $aligner
		fi

		# Update summary -------------------------------------------------------
		echo -e " · Retrieving flagstats ..."
		samtools flagstat $out/"$condition"/"$condition".sorted.bam \
			> $out/"$condition"/"$condition".bam_notes.txt

		perc=`cat $out/"$condition"/"$condition".bam_notes.txt | \
			grep 'mapped' | head -n 1 | cut -d '(' -f 2 | cut -d ':' -f 1`
		perc2=`cat $out/"$condition"/"$condition".bam_notes.txt | \
			grep 'properly paired' | head -n 1 | cut -d '(' -f 2 | \
			cut -d ':' -f 1`
		chim=`cat $out/"$condition"/"$condition".bam_notes.txt | \
			grep 'with mate mapped to a different chr' | head -n 1 | \
			cut -d ' ' -f 1`
		row=`grep "$condition" "$out/summary"`
		echo -e "$row\t$perc\t$perc2\t$chim" >> $outcontrol/summarytmp

		# Add back the UMIs to the SAM file ------------------------------------
		echo -e " · Adding linkers to SAM file ..."
		grep -v "^\@" $out/"$condition"/"$condition".sam | tr -s ' ' | \
			tr '\t' ' ' | sort --parallel=$numbproc \
			--temporary-directory=$HOME/tmp -k1,1 | \
			join - $out/"$condition"/filtered.r1.linkers.oneline.fq \
			-j 1 -o 0,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.2,1.10,2.3,1.11 \
			> $out/"$condition"/"$condition".linkers.sam
	done

	cp $out/summary $outcontrol/summary_pre_align
	mv $outcontrol/summarytmp $out/summary

}
execute_step $opt 'alignment' alignment

# Filter SAM -------------------------------------------------------------------
function filter_sam() {
	echo -e '# SAM filtering ================================================\n'

	# Get mapq threshold by input
	input_int 30 'MAPQ threshold' $mapqthr
	mapqthr=$v

	IFS=',' read -r -a condv <<< "$conds"
	for condition in "${condv[@]}"; do
		echo -e "\nAnalyzing UMIs from condition '$condition'..."

		# Filter SAM
		time $scriptdir/sam_filter.R $out/$condition/ $experiment $condition \
			-mt $mapqthr -cs $cutsite -c $numbproc & pid0=$!
		wait $pid0

	done

}
execute_step $opt 'SAM filtering' filter_sam


# PREPARE UMI ------------------------------------------------------------------
function prepare_umi() {
	echo -e '# UMI grouping&deduplicating ===================================\n'

	# Get cutsite list file by input
	input_fname 'Cutsite list file' 'a list of cutsite positions' $cutsitelist
	cutsitelist=$v

	# Get maskfile by input
	input_fname 'Mask file' 'a list of regions to be masked' $maskfile
	maskfile=$v

	function prepare_umi_single_condition() {
		echo -e "\nPreparing UMIs from condition '$condition'..."
		cslbool=0
		if [[ -n $cutsitelist ]]; then
			cslbool=1
		fi

		# Group UMIs -----------------------------------------------------------
		if [[ -n $maskfile ]]; then
			$scriptdir/umi_group.py $out/$condition/ $experiment $condition \
				$umi_length --mask-file $maskfile & pid0=$!
			wait $pid0
		else
			$scriptdir/umi_group.py $out/$condition/ $experiment $condition \
				$umi_length & pid0=$!
			wait $pid0
		fi

		if [[ -n $cutsitelist ]]; then
			# Assign UMIs to cutsites ------------------------------------------
			$scriptdir/pos2cutsites.R $out/$condition/ $experiment $condition \
				$cutsitelist -i $csrange -c $numbproc & pid0=$!
			wait $pid0
		fi

		if [[ $umi_length -ne 0 ]]; then
			# Deduplicating UMIs -----------------------------------------------
			echo -e "\nDeduplicating UMIs ..."
			$scriptdir/umi_dedupl.R $out/$condition/ $experiment $condition \
				-p $platform -co $pthr -c $numbproc -cs $cslbool \
				-em $emax -ep $eperc & pid0=$!
			wait $pid0
		else
			if [[ -n $cutsitelist ]]; then
				cp $out/$condition/UMIpos.atcs.txt \
					$out/$condition/UMIpos.unique.atcs.txt
			else
				cp $out/$condition/UMIpos.txt $out/$condition/UMIpos.unique.txt
			fi
		fi
	}
	IFS=',' read -r -a condv <<< "$conds"
	for condition in "${condv[@]}"; do
		time prepare_umi_single_condition
	done

}
execute_step $opt 'UMI preparation' prepare_umi

# BIN UMI ----------------------------------------------------------------------
function bin_step() {
	echo -e '# Binning ======================================================\n'

	# Get binsize by input
	input_int 1e6 'Binsize' $binsize
	binsize=$v

	# Get binstep by input
	input_int 1e6 'Binstep' $binstep
	binstep=$v

	# Get cutsite list file by input
	input_fname 'Cutsite list file' 'a list of cutsite positions' $cutsitelist
	cutsitelist=$v

	# Get maskfile by input
	input_fname 'Mask file' 'a list of regions to be masked' $maskfile
	maskfile=$v

	# Get chrlengths by input
	input_fname 'Chr length file' \
		'lengths of chromosomes in the specified genome version' $chrlengths
	chrlengths=$v

	cslbool=0
	if [[ -n $cutsitelist ]]; then
		cslbool=1
	fi

	# Bin cutsites -------------------------------------------------------------
	echo -e "\nBinning cutsites ..."
	time $scriptdir/cs_bin.R -i $binsize -t $binstep -c $numbproc \
		$out $cutsitelist $chrlengths & pid0=$!
	wait $pid0

	# Bin UMIs -----------------------------------------------------------------
	echo -e "\nBinning UMIs ..."
	$scriptdir/umi_bin.R -i $binsize -t $binstep -c $numbproc \
		$out/ $experiment $conds $cutsitelist $chrlengths & pid0=$!
	wait $pid0

}
execute_step $opt 'binning' bin_step

# ANALYZE UMI ------------------------------------------------------------------
function analyze_umi() {
	echo -e '# Plotting =====================================================\n'

	# Get cutsite list file by input
	input_fname 'Cutsite list file' 'a list of cutsite positions' $cutsitelist
	cutsitelist=$v

	cslbool=0
	if [[ -n $cutsitelist ]]; then
		cslbool=1
	fi

	# Get binsize by input
	input_int 1e6 'Binsize' $binsize
	binsize=$v

	# Get binstep by input
	input_int 1e6 'Binstep' $binstep
	binstep=$v

	# Get maskfile by input
	input_fname 'Mask file' 'a list of regions to be masked' $maskfile
	maskfile=$v

	# Get chrlengths by input
	input_fname 'Chr length file' \
		'lengths of chromosomes in the specified genome version' $chrlengths
	chrlengths=$v

	# Make multi-condition plots -----------------------------------------------
	
	# Prepare flags for heterochromosomes removal
	flags=""
	if [ "$rmX" = true ]; then flags="$flags --rmChrX"; fi
	if [ "$rmY" = true ]; then flags="$flags --rmChrY"; fi
	if [[ -n $neg ]]; then flags="$flags --neg $neg"; fi

	# Run plotting script
	$scriptdir/umi_plot.R $flags -i $binsize -t $binstep -c $numbproc \
		$out/ $experiment $conds $cutsitelist $chrlengths $maskfile & pid0=$!
	wait $pid0

}
time execute_step $opt 'UMI analysis' analyze_umi

# END --------------------------------------------------------------------------

echo -e "\n\n~~ fin ~~"

# END LOG ======================================================================

} &> >(tee $logpath)

################################################################################
