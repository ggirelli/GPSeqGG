#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 0.1.0
# Description: trim linker from reads
# 
# ------------------------------------------------------------------------------



# ENV VAR ======================================================================

export LC_ALL=C

# INPUT ========================================================================

# Help string
helps="
 usage: ./reads_trim.sh [-h] -o outDir -c cond -p patFile

 Description:
  Trim linker from reads.

 Mandatory arguments:
  -o outDir	Output directory. Created if not found.
  -c cond	Condition to analyze.
  -p patFile	Pattern file.

 Optional arguments:
  -h		Show this help page.
"

# Parse options
while getopts ho:c:p: opt; do
	case $opt in
		h)
			echo -e "$helps\n"
			exit 1
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
			if [ -e "$OPTARG" ]; then
				patFile=$OPTARG
			else
				msg="Invalid -p option, vile not found.\nFile: $OPTARG"
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
if [ -z "$patFile" ]; then
	echo -e "$helps\n!!! Missing mandatory -p option.\n"
	exit 1
fi

# Additional checks
if [ ! -d "$out_dir/$condition/" ]; then
	msg="Invalid -c option, folder not found.\nFolder: $out_dir/$condition/"
	echo -e "$helps\n$msg"
	exit 1
fi

# RUN ==========================================================================

# Select non-genomic region length
length=`grep "$condition" $patFile | cut -f 3`
echo -e " · Trimming the first $length bases (pattern file)."

# TRIM -----------------------------------------------------------------

# Remove first $length bases/cigar-chars from R1 FASTA file
echo -e " >>> Trimming R1 FASTA file..."
sed -r "s/^[ACTGN]{$length}//" $out_dir/"$condition"/filtered.r1.fa \
	> $out_dir/"$condition"/filtered.r1.noLinker.fa

# Generate temporary oneline fastq file
cat $out_dir/"$condition"/filtered.r1.fq | paste - - - - \
	> $out_dir/"$condition"/filtered.r1.oneline.fq

# Remove first $length bases/cigar-chars from R1 FASTQ file
echo -e " >>> Trimming R1 FASTQ file..."

cat $out_dir/"$condition"/filtered.r1.oneline.fq | cut -f 1 \
	> $out_dir/"$condition"/tmpbase

cat $out_dir/"$condition"/filtered.r1.oneline.fq | cut -f 2 | \
	sed -r "s/^(.{$length})(.*)\$/\2\t+/" - | \
	paste $out_dir/"$condition"/tmpbase - > $out_dir/"$condition"/tmp

cat $out_dir/"$condition"/filtered.r1.oneline.fq | cut -f 4 | \
	sed -r "s/^(.{$length})(.*)\$/\2/" - | \
	paste $out_dir/"$condition"/tmp - | tr "\t" "\n" \
	> $out_dir/"$condition"/filtered.r1.noLinker.fq

# Retain only the first $length bases/cigar-chars from R1 FASTQ file
echo -e " >>> Saving linkers..."

cat $out_dir/"$condition"/filtered.r1.oneline.fq | cut -f 2 | \
	sed -r "s/^(.{$length})(.*)\$/\1\t+/" - | \
	paste $out_dir/"$condition"/tmpbase - > $out_dir/"$condition"/tmp

cat $out_dir/"$condition"/filtered.r1.oneline.fq | cut -f 4 | \
	sed -r "s/^(.{$length})(.*)\$/\1/" - | \
	paste $out_dir/"$condition"/tmp - \
	> $out_dir/"$condition"/filtered.r1.linkers.oneline.fq

rm $out_dir/"$condition"/tmp*

# Remove temporary oneline fastq file
rm $out_dir/"$condition"/filtered.r1.oneline.fq

echo -e " · Trimmed."

# END --------------------------------------------------------------------------

################################################################################
