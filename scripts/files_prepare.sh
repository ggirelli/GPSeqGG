#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 1.0.0
# Description: generates files necessary for GPSeq sequencing data analysis
# 
# ------------------------------------------------------------------------------



# ENV VAR ======================================================================

export LC_ALL=C

# INPUT ========================================================================

# Help string
helps="
 usage: ./files_prepare.sh [-h][-t threads] -o outDir -1 r1file [-2 r2file]

 Description:
  Prepare files for analysis.

 Mandatory arguments:
  -o outDir	Output folder, created if not found.
  -1 r1file	R1 fastq file.

 Optional arguments:
  -h		Show this help page.
  -t threads	Number of threads for parallelization
  -2 r2file	R2 fastq file. 
"

# Default values
threads=1

# Parse options
while getopts ht:o:1:2: opt; do
	case $opt in
		h)
			echo -e "$helps\n"
			exit 0
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
		1)
			if [ -e "$OPTARG" ]; then
				r1=$OPTARG
			else
				msg="Invalid -1 option, file not found.\nFile: $OPTARG"
				echo -e "$helps\n$msg"
				exit 1
			fi
		;;
		2)
			if [ -e "$OPTARG" ]; then
				r2=$OPTARG
			else
				msg="Invalid -2 option, file not found.\nFile: $OPTARG"
				echo -e "$helps\n$msg"
				exit 1
			fi
		;;
	esac
done

# Check mandatory options
if [ -z "$r1" ]; then
	echo -e "$helps\n!!! Missing mandatory -1 option.\n"
	exit 1
fi
if [ -z "$out_dir" ]; then
	echo -e "$helps\n!!! Missing mandatory -o option.\n"
	exit 1
fi

# RUN ==========================================================================

# Decompress -------------------------------------------------------------------
echo -e "Decompress tar.gz..."

# Decompress R1
gunzip -c $r1 | cut -d ' ' -f 1 > $out_dir/r1.fq & pid0=$!

# Decompress R2
if [ -n "$r2" ]; then
	gunzip -c $r2 | cut -d ' ' -f 1 > $out_dir/r2.fq & pid4=$!
fi

# Remove QUALs -----------------------------------------------------------------
echo -e "Remove quality values..."

# R1
wait $pid0
cat $out_dir/r1.fq | paste - - - - | cut -f 1,2 | sed 's/^@/>/' | tr "\t" "\n" | \
	cut -d ' ' -f 1 > $out_dir/r1.fa & pid1=$!

# R2
if [ -n "$r2" ]; then
	wait $pid4
	cat $out_dir/r2.fq | paste - - - - | cut -f 1,2 | sed 's/^@/>/' | \
		tr "\t" "\n" | cut -d ' ' -f 1 > $out_dir/r2.fa & pid5=$!
fi

# One-line FA ------------------------------------------------------------------
echo -e "Generate oneline fasta file..."

# R1
wait $pid1
cat $out_dir/r1.fa | paste - - > $out_dir/r1oneline.fa & pid2=$!

# R2
if [ -n "$r2" ]; then
	wait $pid5
	cat $out_dir/r2.fa | paste - - > $out_dir/r2oneline.fa & pid6=$!
fi

# One-line FQ ------------------------------------------------------------------
echo -e "Generate online fastq file..."
mkdir -p $HOME/tmp

# R1
wait $pid2
cat $out_dir/r1.fq | tr ' ' '_' | paste - - - - | \
	sort --parallel=$threads --temporary-directory=$HOME/tmp -k1,1 \
	> $out_dir/r1oneline.fq & pid3=$!
wait $pid3

# R2
if [ -n "$r2" ]; then
	wait $pid6
	cat $out_dir/r2.fq | tr ' ' '_' | paste - - - - | \
	sort --parallel=$threads --temporary-directory=$HOME/tmp -k1,1 \
	> $out_dir/r2oneline.fq & pid7=$!
	wait $pid7
fi

# End --------------------------------------------------------------------------

################################################################################
