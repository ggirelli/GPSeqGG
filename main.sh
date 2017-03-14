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
scriptdir="`dirname ${BASH_SOURCE}`/scripts/"
moddir="`dirname ${BASH_SOURCE}`/modules/"

# Load functions
source $moddir/main.functions.sh

# CHECK OPTIONS ================================================================

# Help string
helps="
usage: ./main.sh [-h][-w][-t threads] -i inDir -o outDir -e expID -c conditions
 [-n][-a aligner][-g refGenome][-d bwaIndex][-x][-y][-f cutsite][-q mapqThr]
 [-p platform][-u umiLength][-r csRange][-j emax][-k eperc][-z binSize]
 [-b binStep][-l csList][-m maskFile][-s chrLengths]

 Description:
  Run a step-by-step interactive GPSeq sequencing data analysis.

 Required files:
  Requires R1 (and R2 if paired-end sequencing) and a pattern files in the
  input directory (inDir). The patterns file should have a condition per row
  with condition name, pattern (scan_for_mateches format) and pattern length
  separated by tabulations.

 Mandatory arguments:
  -i indir	Input directory.
  -o outdir	Output directory. Created if not found.
  -e expID	Experiment ID.
  -c conditions	Comma-separated conditions.

 Optional arguments:
  -h	Show this help page.
  -w	Perform every step of the pipeline without asking.
  -n	Negative present. Expected label: neg;
  -x	Remove X chromosome after alignment.
  -y	Remove Y chromosome after alignment.
  -t threads	Number of threads for parallelization.
  -a aligner	Aligner. Either 'bwa' (default) or 'bowtie2'.
  -g refGenome	Path to reference genome file. Default: 'hg19'.
  -d bwaIndex	Path to BWA index file. Required if BWA is the selected aligner.
  -f cutsite	Cutsite sequence. Default: 'AAGCTT' (HindIII).
  -q mapqThr	Mapping quality threshold. Default: 30.
  -p platform	Sequencing platform. Default: 'L'.
  -u umilength	UMI sequence length. Default: 8.
  -r csRange	Range around cutsite for UMI assignment. Default: 40.
  -j emax	Maximum error probability for read quality filtering. Default: 1e-3.
  -k eperc	Maximum % of bases with emax error probability. Default: 20.
  -z binSize	Bin size. Default: 1e6.
  -b binStep	Bin step. Default: 1e5.
  -l csList	File with cutsite list. Columns: chr|pos. No header.
  -m maskFile	File with masked regions. Columns: id|chr|start|end. No header.
  -s chrLengths	File with chromosome lengths. chr|len. No header.
"

# Default values
dontask=0
threads=1
rmX=false
rmY=false
aligner='bwa'
refGenome='hg19'
cutsite='AAGCTT'
mapqThr=30
platform='L'
umiLength=8
csRange=40
binSize=1e6
binStep=1e5
pthr=0
emax=1e-3
eperc=20

# Parse options
while getopts hwt:i:o:e:c:ng:a:d:xyf:q:p:u:r:z:b:j:k:l:m:s: opt; do
	case $opt in
		h)
			# Help
			echo -e "$helps"
			exit 0
		;;
		w)
			# Automatic
			dontask=1
		;;
		t)
			# Threads
			if [ 0 -ge "$OPTARG" ]; then
				echo -e "Enforcing a minimum of 1 thread.\n"
			else
				threads=$OPTARG
			fi
		;;
		i)
			# Input directory
			if [ -d "$OPTARG" ]; then
				indir=$OPTARG
			else
				msg="Invalid -i option, folder not found.\nFolder: $OPTARG"
				echo -e "$helps\n!!! $msg"
				exit 1
			fi
		;;
		o)
			# Output directory
			outdir=$OPTARG
			if [ ! -d "$OPTARG" ]; then
				msg="Output folder not found, creating it."
				mkdir -p $outdir
			fi
		;;
		e)
			# Experiment ID
			expID=$OPTARG
		;;
		c)
			# Comma-separated condition string
			conds=$OPTARG
			# Condition array
			IFS=',' read -r -a condv <<< "$conds"
		;;
		n)
			# Negative condition label
			neg='neg'
		;;
		g)
			# File with cutsite list
			if [ -e $OPTARG ]; then
				refGenome=$OPTARG
			else
				msg="Invalid -g option, file not found.\nFile: $OPTARG"
				echo -e "$helps\n!!! $msg"
				exi 1
			fi
		;;
		a)
			# Aligner
			if [ 'bwa' == "$OPTARG" -o 'bowtie2' == "$OPTARG" ]; then
				aligner=$OPTARG
			else
				msg="Invalid -a option. Available values: 'bwa', 'bowtie2'."
				echo -e "$helps\n!!! $msg"
				exit 1
			fi
		;;
		d)
			# File with cutsite list
			if [ -e $OPTARG ]; then
				bwaIndex=$OPTARG
			else
				msg="Invalid -d option, file not found.\nFile: $OPTARG"
				echo -e "$helps\n!!! $msg"
				exi 1
			fi
		;;
		x)
			# Remove chrX after alignment
			rmX=true
		;;
		y)
			# Remove chrY after alignment
			rmY=true
		;;
		f)
			# Cutsite sequence
			cutsite=$OPTARG
		;;
		q)
			# Mapping quality threshold
			mapqThr=$OPTARG
		;;
		p)
			# Sequencing platform
			platform=$OPTARG
		;;
		u)
			# UMI length in nt
			umiLength=$OPTARG
		;;
		r)
			# Range around cutsite
			csRange=$OPTARG
		;;
		z)
			# Bin size in nt
			binSize=$OPTARG
		;;
		b)
			# Bin step nin nt
			binStep=$OPTARG
		;;
		j)
			emax=$OPTARG
		;;
		k)
			eperc=$OPTARG
		;;
		l)
			# File with cutsite list
			if [ -e $OPTARG ]; then
				csList=$OPTARG
			else
				msg="Invalid -e option, file not found.\nFile: $OPTARG"
				echo -e "$helps\n!!! $msg"
				exi 1
			fi
		;;
		m)
			# File with masked regions
			if [ -e $OPTARG ]; then
				maskFile=$OPTARG
			else
				msg="Invalid -m option, file not found.\nFile: $OPTARG"
				echo -e "$helps\n!!! $msg"
				exi 1
			fi
		;;
		s)
			# File with chromosome lengths
			if [ -e $OPTARG ]; then
				chrLengths=$OPTARG
			else
				msg="Invalid -s option, file not found.\nFile: $OPTARG"
				echo -e "$helps\n!!! $msg"
				exi 1
			fi
		;;
	esac
done

# Check mandatory options
if [ -z "$indir" ]; then
	echo -e "$helps\n!!! Missing mandatory -i option.\n"
	exit 1
fi
if [ -z "$outdir" ]; then
	echo -e "$helps\n!!! Missing mandatory -o option.\n"
	exit 1
fi
if [ -z "$expID" ]; then
	echo -e "$helps\n!!! Missing mandatory -e option.\n"
	exit 1
fi
if [ -z "$conds" ]; then
	echo -e "$helps\n!!! Missing mandatory -c option.\n"
	exit 1
fi

# Additional checks
if [ "bwa" == "$aligner" -a -z "$bwaIndex" ]; then
	echo -e "$helps\n!!! Missing mandatory -d option.\n"
	exit 1
fi

# Check settings ---------------------------------------------------------------

# Settings header
settings="
# SETTINGS
"
if [ 1 -eq $dontask ]; then
	settings="$settings
 Perform EVERY step of the pipeline.
"
fi

# Settings mandatory options
settings="$settings
 Input directory:\t$indir
 Output directory:\t$outdir
 
 Experiment ID: $expID
 Conditions: $conds"
if [ -n "$neg" ]; then
	settings="$settings\n Negative: $neg"
else
	settings="$settings\n!!! No negative condition."
fi

# Other settings
settings="$settings

 Threads: $threads

 Cutsite sequence: $cutsite
 Reference genome: $refGenome
 Aligner: $aligner"
if [ -n $bwaIndex ]; then
	settings="$settings\n BWA index: $bwaIndex"
fi

settings="$settings
 MAPQ threshold: $mapqThr
 Platforhm: $platform

 UMI length: $umiLength nt
 Cutsite range: $csRange nt
 Max error probability: $emax
 Max emax bases percentage: $eperc%

 Bin size: $binSize nt
 Bin step: $binStep nt
"

# Additional files
if [ -n "$csList" ]; then
	settings="$settings\n Cutsites:\n $csList\n"
fi
if [ -n "$maskFile" ]; then
	settings="$settings\n Masked regions:\n $maskFile\n"
fi
if [ -n "$chrLengths" ]; then
	settings="$settings\n Chromosome lengths:\n $chrLengths\n"
fi

# Print settings
clear
echo -e "$settings"

# Ask the user to double-check everything
check_settings

# PREPARE DIRECTORY STRUCTURE ==================================================

# Input folders
in=$indir/indata && mkdir -p $in

# Output folders
out=$outdir/ && mkdir -p $out
cout=$outdir/conditions && mkdir -p $cout
pout=$outdir/plots && mkdir -p $pout
xout=$outdir/aux && mkdir -p $xout

# Additional outputs
outcontrol=$outdir/tmp && mkdir -p $outcontrol
logpath="$out/$expID.log"

# START LOG ====================================================================
clear
{

# Start
echo -e "$settings
START\n=====================\n"

# LOAD DATA FILES ==============================================================

find $indir -maxdepth 1 -type f -iname "*$expID*R[12]*" | sort > filelist
numb_of_files=`cat filelist | wc -l`

if [ 0 -eq $numb_of_files ]; then
	echo -e "!!! ERROR. No R1/R2 files found for $expID.\n"
	exit 1
fi

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
	echo -e 'Quality control\n=====================\n'
	# Produce quality control summarie(s)
	if [ -z "$r2" ]; then
		time $scriptdir/quality_control.sh -t $threads -o $xout -1 $r1
	else
		time $scriptdir/quality_control.sh -t $threads -o $xout -1 $r1 -2 $r2
	fi
}
execute_step $dontask 'quality control' quality_control

# FILE GENERATION --------------------------------------------------------------
function file_generation() {
	echo -e 'File generation\n=====================\n'
	# Generate necessary files
	if [ -z "$r2" ]; then
		time $scriptdir/files_prepare.sh -t $threads -o $in -1 $r1
	else
		time $scriptdir/files_prepare.sh -t $threads -o $in -1 $r1 -2 $r2
	fi
}
execute_step $dontask 'file generation' file_generation

# PATTERN FILTERING ------------------------------------------------------------
function pattern_filtering() {
	echo -e 'Pattern filtering\n=====================\n'

	# Count total reads
	echo -e "Counting total reads..."
	total_read_count=`wc -l $in/r1oneline.fa | cut -d " " -f 1`
	header="condition\tpattern\ttotal_read_count\treads_with_prefix"
	header="$header\tprefix/total"
	echo -e $header > $out/summary

	# Work on single conditions
	patfiles="$indir/pat_files"
	for condition in "${condv[@]}"; do
		echo -e "\n> Working on condition '$condition'..."

		# Select pattern
		pattern=`grep -P "$condition\t" $patfiles | cut -f 2`

		# Save condition-specific patfile
		patfile="$cout/$condition/pat_file"
		mkdir -p $cout/$condition
		echo "$pattern" > $patfile
		echo -e "\nPattern: $pattern"

		# Identify condition-specific reads
		time $scriptdir/pattern_filter.sh -t $threads -i $in \
			-o "$cout/$condition" -p $patfile & pid0=$!
		wait $pid0

		# Print condition-specific read count in the summary
		count=`cat $cout/"$condition"/filtered.r1.fa | paste - - | wc -l`

		convstr="scale=4;perc=$count/$total_read_count;scale=2;(perc*100)/1"
		perc=`echo $convstr | bc`

		echo -e $header > $cout/"$condition"/summary
		header="$condition\t$pattern\t$total_read_count\t$count\t$perc%"
		echo -e $header >> $cout/"$condition"/summary
		echo -e $header >> $out/summary
	done

	cp $out/summary $outcontrol/summary_pattern
}
execute_step $dontask 'pattern_filtering' pattern_filtering

# ALIGNMENT --------------------------------------------------------------------
function alignment() {
	echo -e 'Alignment\n=====================\n'

	new_fields="\tmapped\tmapped/prefix\tproperly_paired"
	new_fields="$new_fields\tproperly_paired/prefix\tinter_chr"
	head -n 1 $outcontrol/summary_pattern | \
		awk -v nf="$new_fields" "{ print \$0 nf }" - \
		> $outcontrol/summary_align

	patfiles="$indir/pat_files"
	for condition in "${condv[@]}"; do
		echo -e "Aligning reads from condition '$condition'..."

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
		grep -v "^\@" $cout/"$condition"/"$condition".sam | tr -s ' ' | \
			tr '\t' ' ' | sort --parallel=$threads \
			--temporary-directory=$HOME/tmp -k1,1 | \
			join - $cout/"$condition"/filtered.r1.linkers.oneline.fq \
			-j 1 -o 0,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.2,1.10,2.3,1.11 \
			> $cout/"$condition"/"$condition".linkers.sam
	done

	cp $outcontrol/summary_align $out/summary

}
execute_step $dontask 'alignment' alignment

# Filter SAM -------------------------------------------------------------------
function filter_sam() {
	echo -e 'SAM filtering\n=====================\n'

	# Get mapq threshold by input
	input_int 30 'MAPQ threshold' $mapqThr
	mapqthr=$v

	for condition in "${condv[@]}"; do
		echo -e "\nAnalyzing UMIs from condition '$condition'..."

		# Filter SAM
		time $scriptdir/sam_filter.R $cout/$condition/ $expID $condition \
			-mt $mapqThr -cs $cutsite -c $threads & pid0=$!
		wait $pid0

	done

}
execute_step $dontask 'SAM filtering' filter_sam


# PREPARE UMI ------------------------------------------------------------------
function prepare_umi() {
	echo -e 'UMI grouping & deduplicating\n=====================\n'

	# Get cutsite list file by input
	input_fname 'Cutsite list file' 'a list of cutsite positions' $csList
	cutsitelist=$v

	# Get maskfile by input
	input_fname 'Mask file' 'a list of regions to be masked' $maskFile
	maskfile=$v

	function prepare_umi_single_condition() {
		echo -e "\nPreparing UMIs from condition '$condition'..."
		cslbool=0
		if [[ -n $csList ]]; then
			cslbool=1
		fi

		# Group UMIs -----------------------------------------------------------
		if [[ -n $maskFile ]]; then
			$scriptdir/umi_group.py $cout/$condition/ $expID $condition \
				$umiLength --mask-file $maskFile & pid0=$!
			wait $pid0
		else
			$scriptdir/umi_group.py $cout/$condition/ $expID $condition \
				$umiLength & pid0=$!
			wait $pid0
		fi

		if [[ -n $csList ]]; then
			# Assign UMIs to cutsites ------------------------------------------
			$scriptdir/pos2cutsites.R $cout/$condition/ $expID $condition \
				$csList -i $csRange -c $threads & pid0=$!
			wait $pid0
		fi

		if [[ $umiLength -ne 0 ]]; then
			# Deduplicating UMIs -----------------------------------------------
			echo -e "\nDeduplicating UMIs ..."
			$scriptdir/umi_dedupl.R $cout/$condition/ $expID $condition \
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
	}
	for condition in "${condv[@]}"; do
		time prepare_umi_single_condition
	done

}
execute_step $dontask 'UMI preparation' prepare_umi

# BIN UMI ----------------------------------------------------------------------
function bin_step() {
	echo -e 'Binning\n=====================\n'

	# Get binSize by input
	input_int 1e6 'Binsize' $binSize
	binSize=$v

	# Get binStep by input
	input_int 1e6 'Binstep' $binStep
	binStep=$v

	# Get cutsite list file by input
	input_fname 'Cutsite list file' 'a list of cutsite positions' $csList
	cutsitelist=$v

	# Get maskfile by input
	input_fname 'Mask file' 'a list of regions to be masked' $maskFile
	maskfile=$v

	# Get chrlengths by input
	input_fname 'Chr length file' \
		'lengths of chromosomes in the specified genome version' $chrLengths
	chrlengths=$v

	cslbool=0
	if [[ -n $csList ]]; then
		cslbool=1
	fi

	# Bin cutsites -------------------------------------------------------------
	echo -e "\nBinning cutsites ..."
	time $scriptdir/cs_bin.R -i $binSize -t $binStep -c $threads \
		$xout $csList $chrLengths & pid0=$!
	wait $pid0

	# Bin UMIs -----------------------------------------------------------------
	echo -e "\nBinning UMIs ..."
	if [ -z "$csList" ]; then
		$scriptdir/umi_bin.R -i $binSize -t $binStep -p $threads \
			$out/ $expID $conds $chrLengths & pid0=$!
	else
		$scriptdir/umi_bin.R -c -i $binSize -t $binStep -p $threads \
			$out/ $expID $conds $chrLengths & pid0=$!
	fi
	wait $pid0

}
execute_step $dontask 'binning' bin_step

# ANALYZE UMI ------------------------------------------------------------------
function analyze_umi() {
	echo -e 'Plotting\n=====================\n'

	# Get cutsite list file by input
	input_fname 'Cutsite list file' 'a list of cutsite positions' $csList
	cutsitelist=$v

	cslbool=0
	if [[ -n $csList ]]; then
		cslbool=1
	fi

	# Get binSize by input
	input_int 1e6 'Binsize' $binSize
	binSize=$v

	# Get binStep by input
	input_int 1e6 'Binstep' $binStep
	binStep=$v

	# Get maskfile by input
	input_fname 'Mask file' 'a list of regions to be masked' $maskFile
	maskFile=$v
	if [ -z "$maskFile" ]; then
		msg="A list of masked regions is required for the plot step."
		msg="$msg\nExit."
		echo -e "$msg"
		exit 1
	fi

	# Get chrlengths by input
	input_fname 'Chr length file' \
		'lengths of chromosomes in the specified genome version' $chrLengths
	chrLengths=$v
	if [ -z "$chrLengths" ]; then
		msg="A list of chromosome lengths is required for the plot step."
		msg="$msg\nExit."
		echo -e "$msg"
		exit 1
	fi

	# Make multi-condition plots -----------------------------------------------
	
	# Prepare flags for heterochromosomes removal
	flags=""
	if [ "$rmX" = true ]; then flags="$flags --rmChrX"; fi
	if [ "$rmY" = true ]; then flags="$flags --rmChrY"; fi
	if [[ -n $neg ]]; then flags="$flags --neg $neg"; fi

	# Run plotting script
	$scriptdir/umi_plot.R $flags -i $binSize -t $binStep -c $threads \
		$out/ $expID $conds $csList $chrLengths $maskFile & pid0=$!
	wait $pid0

}
time execute_step $dontask 'UMI analysis' analyze_umi

# END --------------------------------------------------------------------------

echo -e "\n\n~~ fin ~~"

# END LOG ======================================================================

} &> >(tee $logpath)

################################################################################
