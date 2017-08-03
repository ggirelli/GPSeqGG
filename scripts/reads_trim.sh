#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 1.1.0
# Description: trim linker from reads
# 
# ------------------------------------------------------------------------------



# ENV VAR ======================================================================

export LC_ALL=C

# INPUT ========================================================================

# Help string
helps="
 usage: ./reads_trim.sh [-h] -o outDir -e expID -c cond -p patFile

 Description:
  Trim linker from reads.

 Mandatory arguments:
  -o outDir	Folder containing conditions.
  -e expID	Experiment ID (same as in patFile).
  -c cond	Condition to analyze.
  -p patFile	Pattern file.

 Optional arguments:
  -h		Show this help page.
"

# Parse options
while getopts ho:e:c:p: opt; do
	case $opt in
		h)
			echo -e "$helps\n"
			exit 1
		;;
		o)
			out_dir=$OPTARG
			if [ ! -d "$OPTARG" ]; then
				msg="\n!!! Invalid -o value.\n"
				echo -e $helps$msg"Output folder not found!"
				exit 1
			fi
		;;
		e)
			expID=$OPTARG
		;;
		c)
			condition=$OPTARG
		;;
		p)
			if [ -e "$OPTARG" ]; then
				patFile=$OPTARG
			else
				msg="Invalid -p option, file not found.\nFile: $OPTARG"
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
if [ -z "$expID" ]; then
	echo -e "$helps\n!!! Missing mandatory -e option.\n"
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

# Shortcuts
d="$out_dir/$condition/"

# RUN ==========================================================================

# Select non-genomic region length
length=`grep -P "^$expID\t$condition\t" $patFile | cut -f 4`
echo -e " · Trimming the first $length bases (pattern file)."

# TRIM -----------------------------------------------------------------

# Generate temporary oneline fastq file
echo -e " >>> Generate temporary oneline FASTQ file..."
cat $d/filtered.r1.fq | paste - - - - \
	> $d/filtered.r1.oneline.fq

# Remove first $length bases/cigar-chars from R1 FASTQ file
echo -e " >>> Trimming R1 FASTQ/FASTA files & saving linkers..."

# Setup AWK trimming program
awkcommand='{
print substr($1,2) "\t" substr($2,0,len) "\t" $3 "\t" substr($4,0,len) > o1;
print $1 "\t" substr($2,len+1) "\t" $3 "\t" substr($4,len+1) > o2;
print $1 "\t" substr($2,len+1) > o3;
}'

# Run AWK
awk -v len=$length \
	-v o1="$d/filtered.r1.linkers.oneline.fq" \
	-v o2="$d/filtered.r1.noLinker.oneline.fq" \
	-v o3="$d/filtered.r1.noLinker.oneline.fa" \
	"$awkcommand" $d/filtered.r1.oneline.fq

# Split lines
cat $d/filtered.r1.noLinker.oneline.fq | tr "\t" "\n" \
	> $d/filtered.r1.noLinker.fq
cat $d/filtered.r1.noLinker.oneline.fa | tr "\t" "\n" \
	> $d/filtered.r1.noLinker.fa

# Remove temporary oneline files
rm $d/filtered.r1.noLinker.oneline.f*
rm $d/filtered.r1.oneline.fq

echo -e " · Trimmed."

# END --------------------------------------------------------------------------

################################################################################
