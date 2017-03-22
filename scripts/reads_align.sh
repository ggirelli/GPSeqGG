#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 0.1.0
# Description: align reads
# 
# ------------------------------------------------------------------------------




# ENV VAR ======================================================================

export LC_ALL=C

# INPUT ========================================================================

# Help string
helps="
 usage: ./reads_align.sh [-h][-t threads] -o outDir -c cond
 	[-p][-r ref][-a aln][-i index]

 Description:
  Align reads using either bowtie2 or bwa mem.

 Mandatory arguments:
  -o outDir	Output directory. Created if not found.
  -c cond	Condition to analyze.

 Optional arguments:
  -h		Show this help page.
  -t threads	Number of threads for parallelization.
  -p		Option for paired-end sequencing.
  -r ref	Reference ref. Default: hg19.
  -a aln	Aligner. Either 'bwa' or 'bowtie2'. Default: 'bwa'.
  -i index	Path to BWA index.
  		Default: '\$DATA'/BiCro-Resources/genomes/'\$ref'bwt/'\$ref'.fa
"

# Default values
threads=1
paired=0
ref='hg19'
aligner='bwa'

# Parse options
while getopts ht:o:c:pr:a:i: opt; do
	case $opt in
		h)
			echo -e "$helps\n"
			exit 1
		;;
		t)
			if [ 0 -ge "$OPTARG" ]; then
				echo -e "Enforcing a minimum of 1 thread.\n"
			else
				threads=$OPTARG
			fi
		;;
		o)
			out_dir=$OPTARG
			if [ ! -d "$OPTARG" ]; then
				msg="Output folder not found, creating it."
				mkdir -p $out_dir
			fi
		;;
		c)
			condition=$OPTARG
		;;
		p)
			paired=1
		;;
		r)
			ref=$OPTARG
		;;
		a)
			if [ 'bwa' == "$OPTARG" -o 'bowtie2' == "$OPTARG" ]; then
				aligner=$OPTARG
			else
				msg="Invalid -a option. Available values: 'bwa', 'bowtie2'."
				echo -e "$helps\n$msg"
				exit 1
			fi
		;;
		i)
			if [ -e $OPTARG ]; then
				bwaIndex=$OPTARG
			else
				msg="Invalid -i option, file not found.\nFile: $OPTARG"
				echo -e "$helps\n$msg"
				exit 1
			fi
		;;
	esac
done

# Check mandatory options
if [ -z "$out_dir" ]; then
	echo -e "$helps\n!!! Missing mandatory -o option.\n"
	exit 1
fi
if [ -z "$condition" ]; then
	echo -e "$helps\n!!! Missing mandatory -c option.\n"
	exit 1
fi

# Additional checks
if [ -z "$bwaIndex" ]; then
	bwaIndex="$DATA"/BiCro-Resources/genomes/"$ref"bwt/"$ref".fa
fi

# RUN ==========================================================================

# Paired-end alignment ---------------------------------------------------------
if [[ 1 -eq $paired ]]; then
	echo " · Performing paired-end alignment ..."
	# BWA alignment
	if [ -n "$aligner" -a "bwa" == "$aligner" ]; then
		bwa mem -t $threads $bwaIndex $out_dir/filtered.r1.noLinker.fq \
			$out_dir/filtered.r2.fq > $out_dir/"$condition".sam \
			2> $out_dir/bwa.log

	# Bowtie2 alignment
	elif [ "bowtie2" == "$aligner" ]; then
		bowtie2 -q -p $threads -1 $out_dir/filtered.r1.noLinker.fq \
			-2 $out_dir/filtered.r2.fq -S $out_dir/"$condition".sam -x $ref \
			2> $out_dir/bowtie2.log
	else
		echo -e "ERROR: unrecognized aligner."
		exit 1
	fi

# Single-end alignment ---------------------------------------------------------
else
	echo " · Performing single-end alingment ..."

	# BWA alignment
	if [ -n "$aligner" -a "bwa" == "$aligner" ]; then
		bwa mem -t $threads $bwaIndex $out_dir/filtered.r1.noLinker.fq \
			> $out_dir/"$condition".sam 2> $out_dir/bwa.log

	# Bowtie2 alignment
	elif [[ 'bowtie2' == "$aligner" ]]; then
		
		bowtie2 -q -p $threads -U $out_dir/filtered.r1.noLinker.fq \
			-S $out_dir/"$condition".sam -x $ref 2> $out_dir/bowtie2.log
	else
		echo -e "ERROR: unrecognized aligner."
		exit 1
	fi
fi

# Save log ---------------------------------------------------------------------
if [ -n "$aligner" -o "bwa" == "$aligner" ]; then
	cat $out_dir/bwa.log
elif [[ 'bowtie2' == "$aligner" ]]; then
	cat $out_dir/bowtie2.log
else
	echo -e "ERROR: unrecognized aligner."
	exit 1
fi

# Generate BAM -----------------------------------------------------------------
echo -e " · Generating and sorting BAM file ..."
samtools sort -@ $threads -o $out_dir/"$condition".sorted.bam \
	$out_dir/"$condition".sam

echo -e " · Indexing BAM file ..."
samtools index $out_dir/"$condition".sorted.bam $out_dir/"$condition".sorted.bai

# END --------------------------------------------------------------------------

################################################################################
