#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 0.1.0
# Description: run fastqc tool
# 
# Usage:
# ./quality_control numbOfiles ncores outFolder r1file [r2file]
# 
# ------------------------------------------------------------------------------



# ENV VAR ======================================================================

export LC_ALL=C

# INPUT ========================================================================

# Help string
helps="
 usage: ./quality_control.sh [-h][-t threads] -o outDir -1 r1file [-2 r2file]

 Description:
  Run FASTQC tool.

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
	exit 0
fi
if [ -z "$out_dir" ]; then
	echo -e "$helps\n!!! Missing mandatory -o option.\n"
	exit 0
fi

# RUN ==========================================================================

if [ -z $r2 ]; then

	# Single-end FASTQC
	fastqc -t $threads -o $out_dir --nogroup $r1

else

	# Pair-end FASTQC
	parallel fastqc -o $out_dir -t $threads --nogroup {} ::: $r1 $r2

fi

# END --------------------------------------------------------------------------

################################################################################
