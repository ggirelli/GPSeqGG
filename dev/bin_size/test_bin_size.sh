#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 1.0.0
# Date: 20170821
# Project: GPSeq - centrality estimation
# Description: estimate genomic region nuclear centrality.
# 
# ------------------------------------------------------------------------------



# ENV VAR ======================================================================

export LC_ALL=C

# PARAMS =======================================================================


# Help string
helps="
usage: ./estimate_centrality.sh [-h][-n binSize_min][-x binSize_max]
                                [-e binSize_step][-p binStep]
                                -o outdir -i bedfile

 Description:
  Identify smallest bin-size that provides a restriction probability
  distribution comparable with similar bin-sizes. Regions in the bed file are
  assigned to a bin if their middle point is inside it. Left-boundary is
  inclusive. Require bedtools for intersection

 Mandatory arguments:
  -o outdir        Output folder.
  -i bedfile       GPSeq condition bedfile.

 Optional arguments:
  -h	Show this help page.
  -n binSize_min   Default to 100 kbp.
  -x binSize_max   Default to 10 Mbp.
  -e binSize_step  Step between bin sizes. Default to 100e3 bp (100 kbp).
  -p binStep       Percentage of bin size to use for the bin step, as moving
                   window. Set to 1 for non-overlapping bins. Default to 0.1
                   (10% bin size).
"

# Default values
binSize_min=100000
binSize_max=10000000
binSize_step=100000
binStep=0.1

# Parse options
while getopts hi:o:n:x:e:p: opt; do
	case $opt in
		h)
			# Help page
			echo -e "$helps"
			exit 0
		;;
		i)
			if [ -e $OPTARG ]; then
				bedfile=$OPTARG
			else
				msg="!!!ERROR! The provided bed file does not exist.\n"
				msg="$msg          File: $OPTARG"
				echo -e "$help\n$msg"
				exit 1
			fi
		;;
		o)
			if [ ! -d $OPTARG ]; then
				mkdir -p $OPTARG
			fi
			outdir=$OPTARG
		;;
		n)
			if [ 0 -lt $OPTARG ]; then
				binSize_min=$OPTARG
			else
				msg="!!!ERROR! Invalid -n option value."
				msg="$msg binSize_min must be > 0."
				echo -e "$help\n$msg"
				exit 1
			fi
		;;
		x)
			if [ 0 -lt $OPTARG ]; then
				binSize_max=$OPTARG
			else
				msg="!!!ERROR! Invalid -x option value."
				msg="$msg binSize_max must be > 0."
				echo -e "$help\n$msg"
				exit 1
			fi
		;;
		e)
			if [ 0 -lt $OPTARG ]; then
				binSize_step=$OPTARG
			else
				msg="!!!ERROR! Invalid -e option value."
				msg="$msg binSize_step must be > 0."
				echo -e "$help\n$msg"
				exit 1
			fi
		;;
		p)
			if [ 0 -lt $OPTARG ]; then
				binStep=$OPTARG
			else
				msg="!!!ERROR! Invalid -p option value."
				msg="$msg binStep must be > 0."
				echo -e "$help\n$msg"
				exit 1
			fi
		;;
		?)
			msg="!!! Unrecognized option."
			echo -e "$help\n$msg"
			exit 1
		;;
	esac
done

# Check mandatory options
if [ -z "$outdir" ]; then
	echo -e "$helps\n!!! ERROR! Missing mandatory -o option.\n"
	exit 1
fi
if [ -z "$bedfile" ]; then
	echo -e "$helps\n!!! ERROR! Missing mandatory -i option.\n"
	exit 1
fi
if [ ! -x "$(command -v bedtools)" -o -z "$(command -v bedtools)" ]; then
	echo -e "$helps\n!!! ERROR! Missing bedtools.\n"
	exit 1
fi

# Additional checks
if [ $binSize_max -le $binSize_min ]; then
	msg="!!!ERROR! -x must be greater than -n.\n"
	echo -e "$help\n$msg"
	exit 1
fi

# Print settings

settings="
  Bin size min : $binSize_min
  Bin size max : $binSize_max
 Bin size step : $binSize_step
      Bin step : $binStep

    Output dir : $outdir
      Bed file : $bedfile
"

echo -e "$settings\n"

# RUN ==========================================================================

# 0) Identify chromosome sizes -------------------------------------------------
echo -e " Retrieving chromosome sizes ..."
chrSize=$(cat $bedfile | grep -v 'track' | datamash -sg1 -t$'\t' max 3)

# Sort chromosomes
awk_add_chr_id='
BEGIN {
	OFS = FS = "\t";
	convert["X"] = 23;
	convert["Y"] = 24;
}
{
	chrid = substr($1, 4);
	if ( chrid in convert ) {
		chrid = convert[chrid];
	}
	print chrid OFS $0;
}'
echo -e "$chrSize" | awk "$awk_add_chr_id" | sort -k1,1n | cut -f2,3 \
	> "$outdir/chr_size.tsv"


# 1) Generate bin bed file -----------------------------------------------------
echo -e " Generating bins ..."

for binSize in $(seq $binSize_min $binSize_step $binSize_max); do
	current_step=$(bc <<< "$binSize * $binStep")

	prefix="bins.size$binSize.step$current_step"
	awk_mk_bins='
	BEGIN {
		OFS = FS = "\t";
	}
	{
		for ( i = 0; i < $2; i += step ) {
			print $1 OFS i OFS i+size;
		}
	}'
	cat "$outdir/chr_size.tsv" | \
		awk -v size=$binSize -v step=$current_step "$awk_mk_bins" \
		> "$outdir/$prefix.bed"
done

# 2) Intersect with bedtools ---------------------------------------------------
echo -e " Intersecting ..."

for binSize in $(seq $binSize_min $binSize_step $binSize_max); do
	current_step=$(bc <<< "$binSize * $binStep")
	prefix="bins.size$binSize.step$current_step"
	fname=$(echo -e "$bedfile" | tr "/" "\t" | awk '{ print $NF }')
	
	echo -e " > Intersecting $prefix ..."
	bedtools intersect -a "$outdir/$prefix.bed" \
		 -b "$bedfile" -wa -wb \
		> "$outdir/intersected.$prefix.$fname.tsv"
done

# 3) Sum up and normalize ------------------------------------------------------

# 4) Compare distributions -----------------------------------------------------

# 5) Plot ----------------------------------------------------------------------

# END ==========================================================================

################################################################################
