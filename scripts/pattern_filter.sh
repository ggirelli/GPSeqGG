#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 0.1.0
# Description: selects reads based on the provided pattern
# Requires: scan_for_matches
# 
# ------------------------------------------------------------------------------



# ENV VAR ======================================================================

export LC_ALL=C

# INPUT ========================================================================

# Help string
helps="
 usage: ./pattern_filter.sh [-h][-t threads] -i inDir -o outDir -p patFile

 Description:
  Select reads based on the provided pattern by runnin scan_for_matches.

 Mandatory arguments:
  -i inDir	Input directory.
  -o outDir	Output directory. Created if not found.
  -p patFile	Pattern file.

 Optional arguments:
  -h		Show this help page.
  -t threads	Number of threads for parallelization.
"

# Default values
threads=1

# Parse options
while getopts ht:i:o:p: opt; do
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
		i)
			if [ -d "$OPTARG" ]; then
				in_dir=$OPTARG
			else
				msg="Invalid -i option, folder not found.\nFolder: $OPTARG"
				echo -e "$helps\n$msg"
				exit 1
			fi
		;;
		o)
			out_dir=$OPTARG
			if [ ! -d "$OPTARG" ]; then
				msg="Output folder not found, creating it."
				mkdir -p $out_dir
			fi
		;;
		p)
			if [ -e "$OPTARG" ]; then
				pat_file=$OPTARG
			else
				msg="Invalid -p option, file not found.\nFile: $OPTARG"
				echo -e "$helps\n$msg"
				exit 1
			fi
		;;
	esac
done

# Check mandatory options
if [ -z "$in_dir" ]; then
	echo -e "$helps\n!!! Missing mandatory -i option.\n"
	exit 1
fi
if [ -z "$out_dir" ]; then
	echo -e "$helps\n!!! Missing mandatory -o option.\n"
	exit 1
fi
if [ -z "$pat_file" ]; then
	echo -e "$helps\n!!! Missing mandatory -p option.\n"
	exit 1
fi

# RUN ==========================================================================

echo "Filtering R1 reads based on patterns ..."

# Identify reads that match the provided pattern
cat $in_dir/r1.fa | parallel --tmpdir $HOME/tmp --block 100M -k \
	--pipe -L 2 "scan_for_matches $pat_file - " | tr -d ' '  \
	> $out_dir/filtered.r1.fa & pid3=$!
wait $pid3

echo "Generating filtered R1 fastq file ..."
# Use FA to generate FQ
cat $out_dir/filtered.r1.fa | grep '>' | tr '>' '@' | \
	cut -d ':' -f -7 | sort --parallel=$threads \
	--temporary-directory=$HOME/tmp -k1,1 | \
	join -t $'\t' - $in_dir/r1oneline.fq | tr '\t' '\n' \
	> $out_dir/filtered.r1.fq & pid4=$!
wait $pid4

# Filter R2 if present
nr2=`find $in_dir -maxdepth 1 -type f -iname "*r2*" | wc -l`
if [[ 0 -lt "$nr2" ]]; then

	echo "Filtering R2 reads based on patterns ..."
	cat $out_dir/filtered.r1.fa | grep '>' | cut -d ':' -f -7 | \
		sed 's/$/\t/' | grep -F -f - $in_dir/r2oneline.fa | \
		sort --parallel=$threads \
		--temporary-directory=$HOME/tmp -k1,1 | tr '\t' '\n' \
		> $out_dir/filtered.r2.fa & pid5=$!

	wait $pid5
	echo "Generating filtered R2 fastq file ..."
	cat $out_dir/filtered.r2.fa | grep '>' | tr '>' '@' | cut -d ':' -f -7 | \
		sort --parallel=$threads --temporary-directory=$HOME/tmp -k1,1 | \
		join -t $'\t' - $in_dir/r2oneline.fq | \
		tr '\t' '\n' > $out_dir/filtered.r2.fq & pid6=$!
	wait $pid6
fi

# END --------------------------------------------------------------------------

################################################################################
