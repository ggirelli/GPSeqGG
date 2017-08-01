#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 1.1.0
# Description: module for main script input option definition.
# 
# ------------------------------------------------------------------------------



# START ========================================================================

# Help string
helps="
usage: ./main.sh [-h][-w][-t threads] -i inDir -o outDir
 [-a aligner][-g refGenome][-d bwaIndex][-x][-y][-q mapqThr]
 [-p platform][-u umiLength][-r csRange][-j emax][-k eperc][-z binSize]
 [-b binStep][-l csList][-m maskFile][-s chrLengths][-n neatness]

 Description:
  Run a step-by-step interactive GPSeq sequencing data analysis.

 Required files:
  Requires R1 (and R2 if paired-end sequencing) and a pattern files in the
  input directory (inDir). The patterns file should have a condition per row
  with condition name, pattern (scan_for_matches format), cutsite sequence and
  non-genomic portion length, separated by tabulations.

 Mandatory arguments:
  -i indir	Input directory.
  -o outdir	Output directory. Created if not found.

 Optional arguments:
  -h	Show this help page.
  -w	Perform every step of the pipeline without asking.
  -x	Remove X chromosome after alignment.
  -y	Remove Y chromosome after alignment.
  -t threads	Number of threads for parallelization.
  -g refGenome	Path to reference genome file. Default: 'hg19'.
  -a aligner	Aligner. Either 'bwa' (default) or 'bowtie2'.
  -d bwaIndex	Path to BWA index file. Required if BWA is the selected aligner.
  -q mapqThr	Mapping quality threshold. Default: 30.
  -p platform	Sequencing platform. Default: 'L'.
  -u umilength	UMI sequence length. Default: 8.
  -r csRange	Range around cutsite for UMI assignment. Default: 40.
  -j emax	Maximum error probability for read quality filtering. Default: 1e-3.
  -k eperc	Maximum % of bases with emax error probability. Default: 20.
  -z binSize	Bin size. Default: 1e6.
  -b binStep	Bin step. Default: 1e5.
  -c cutsite	Cutsite sequence, needed for read re-position. Default: AAGCTT.
  -l csList	File with cutsite list. Columns: chr|pos. No header.
  -m maskFile	File with masked regions. Columns: id|chr|start|end. No header.
  -s chrLengths	File with chromosome lengths. chr|len. No header.
  -n neatness	Neatness level: 0 (heavy), 1 (light), 2 (lightest). Default: 1.
"

# Default values
dontask=false
threads=1
rmX=false
rmY=false
aligner='bwa'
refGenome='hg19'
mapqThr=30
platform='L'
umiLength=8
csRange=40
binSize=1e6
binStep=1e5
pthr=0
emax=1e-3
eperc=20
cutsite="AAGCTT"
neatness=1

neatness_labels=()
neatness_labels+=('heavy')
neatness_labels+=('light')
neatness_labels+=('lightest')


# Parse options
while getopts hwxyi:o:t:g:a:d:q:p:u:r:j:k:z:b:c:l:m:s:n: opt; do
	case $opt in
		# Help
		h)
			echo -e "$helps"
			exit 0
		;;
		# Automatic
		w) dontask=true;;
		# Remove chrX after alignment
		x) rmX=true;;
		# Remove chrY after alignment
		y) rmY=true;;
		# Input directory
		i)
			if [ -d "$OPTARG" ]; then
				indir=$OPTARG
			else
				msg="Invalid -i option, folder not found.\nFolder: $OPTARG"
				echo -e "$helps\n!!! ERROR! $msg"
				exit 1
			fi
		;;
		# Output directory
		o)
			outdir=$OPTARG
			if [ ! -d "$OPTARG" ]; then
				msg="Output folder not found, creating it."
				mkdir -p $outdir
			fi
		;;
		# Threads
		t)
			if [ 0 -ge "$OPTARG" ]; then
				echo -e "Enforcing a minimum of 1 thread.\n"
			else
				threads=$OPTARG
			fi
		;;
		# Reference genome
		g) refGenome=$OPTARG;;
		# Aligner
		a)
			if [ 'bwa' == "$OPTARG" -o 'bowtie2' == "$OPTARG" ]; then
				aligner=$OPTARG
			else
				msg="Invalid -a option. Available values: 'bwa', 'bowtie2'."
				echo -e "$helps\n!!! ERROR! $msg"
				exit 1
			fi
		;;
		# File with cutsite list
		d)
			if [ -e $OPTARG ]; then
				bwaIndex=$OPTARG
			else
				msg="Invalid -d option, file not found.\nFile: $OPTARG"
				echo -e "$helps\n!!! ERROR! $msg"
				exi 1
			fi
		;;
		# Mapping quality threshold
		q) mapqThr=$OPTARG;;
		# Sequencing platform
		p) platform=$OPTARG;;
		# UMI length in nt
		u) umiLength=$OPTARG;;
		# Range around cutsite
		r) csRange=$OPTARG;;
		# Maximum error probability for read quality filtering
		j) emax=$OPTARG;;
		# Maximum % of bases with emax error probability
		k) eperc=$OPTARG;;
		# Bin size in nt
		z) binSize=$OPTARG;;
		# Bin step nin nt
		b) binStep=$OPTARG;;
		# Cutsite sequence
		c) cutsite=$OPTARG;;
		# File with cutsite list
		l)
			if [ -e $OPTARG ]; then
				csList=$OPTARG
			else
				msg="Invalid -l option, file not found.\nFile: $OPTARG"
				echo -e "$helps\n!!! ERROR! $msg"
				exit 1
			fi
		;;
		# File with masked regions
		m)
			if [ -e $OPTARG ]; then
				maskFile=$OPTARG
			else
				msg="Invalid -m option, file not found.\nFile: $OPTARG"
				echo -e "$helps\n!!! ERROR! $msg"
				exit 1
			fi
		;;
		# File with chromosome lengths
		s)
			if [ -e $OPTARG ]; then
				chrLengths=$OPTARG
			else
				msg="Invalid -s option, file not found.\nFile: $OPTARG"
				echo -e "$helps\n!!! ERROR! $msg"
				exit 1
			fi
		;;
		# Neatness level
		n)
			if [ 0 -gt "$OPTARG" ]; then
				neatness=0
			else
				if [ 2 -lt "$OPTARG" ]; then
					neatness=2
				else
					neatness=$OPTARG
				fi
			fi
		;;
		# Stop if unknown option is provided
		?)
			echo -e "$helps\n!!! ERROR! Illegal option -- ${@:$OPTIND-1:1}"
			exit 1 
		;;
	esac
done

# Stop if positional option was found
if [[ -n "${@:$OPTIND:1}" ]]; then
	echo -e "$helps\n!!! ERROR! Illegal option -- ${@:$OPTIND:1}"
	exit 1
fi

# Check mandatory options
if [ -z "$indir" ]; then
	echo -e "$helps\n!!! ERROR! Missing mandatory -i option.\n"
	exit 1
fi
if [ -z "$outdir" ]; then
	echo -e "$helps\n!!! ERROR! Missing mandatory -o option.\n"
	exit 1
fi

# Additional checks
if [ "bwa" == "$aligner" -a -z "$bwaIndex" ]; then
	echo -e "$helps\n!!! ERROR! Missing mandatory -d option.\n"
	exit 1
fi

# Check patterns.tsv --------------------------------------------------------------

if [ -e "$indir/patterns.tsv" ]; then

	# Count columns
	for i in $(seq 1 `cut -f 1 $indir/patterns.tsv | wc -l`); do
		# Retrieve pattern row
		row=`cat $indir/patterns.tsv | head -n $i | tail -n 1`

		# Count fields
		nc=`echo "$row" | awk '{ n=split($0, t, "\t"); print n; }' | uniq`

		# Check field number
		if [ 4 -ne $nc ]; then
			msg="!!! Wrong column number in patterns.tsv, row $i.\n"
			msg="$msg    Expected 4 columns, found $nc."
			echo -e "$helps\n$msg"
			exit 1
		fi
	done

	expv=$(cut -f 1 $indir/patterns.tsv | sort | uniq)
else
	fomt="experiment_id\tcondition_label\tlinker_pattern"
	fomt=$fomt"\ttrim_length\tcutsite_seq"
	echo -e $fomt > $indir/patterns.tsv

	msg="!!! ERROR! Missing patterns.tsv.\n"
	msg=$msg"    Created patterns.tsv with standard format:\n"
	msg=$msg"    $fomt\n"
	echo -e "$helps\n$msg"
	exit 1
fi

# Check settings ---------------------------------------------------------------

# Settings header
settings="
# SETTINGS
"
if $dontask; then
	settings="$settings
 Perform EVERY step of the pipeline.
"
fi

# Settings mandatory options
settings="$settings
 Input directory:\t$indir
 Output directory:\t$outdir"

settings=$settings"\n\n Pattern instructions:\n"

patterns=" experiment_id\tcondition_label\tlinker_pattern\ttrim_length\n"
patterns=$(echo -e $patterns$(cat -t $indir/patterns.tsv | tr '\n' '\t' | \
	sed 's/\t/\\n/g' | sed 's/\^I/\\t/g') | column -c 5 -s $'\t' -t)
settings=$settings"$patterns"

# Other settings
settings="$settings

 Threads: $threads

 Reference genome: $refGenome
 Aligner: $aligner"
if [ -n $bwaIndex ]; then
	settings="$settings\n BWA index: $bwaIndex"
fi

settings="$settings
 MAPQ threshold: $mapqThr
 Platforhm: $platform

 UMI length: $umiLength nt
 Cutsite seq: $cutsite
 Cutsite range: $csRange nt
 Max error probability: $emax
 Max emax bases percentage: $eperc%

 Bin size: $binSize nt
 Bin step: $binStep nt

 Neatness level: $neatness (${neatness_labels[$neatness]})
"
# Chromosome removal
if $rmX; then
	settings="$settings!!! Removing chrX after alignment.\n"
fi
if $rmY; then
	settings="$settings!!! Removing chrY after alignment.\n"
fi

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

################################################################################
