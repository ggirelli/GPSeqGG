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
usage: ./estimate_centrality.sh [-h][-s binSize][-p binStep]
                                -o outdir -c csBed [BEDFILE]...

 Description:
  Calculate global centrality metrics. Requires bedtools for bin assignment.

 Mandatory arguments:
  -o outdir     Output folder.
  -c csBed      Cutsite bedfile.
  BEDFILE       At least two (2) GPSeq condition bedfiles, space-separated and
                in increasing order of restriction conditions intensity.
                Expected to be ordered per condition.

 Optional arguments:
  -h	Show this help page.
  -s binSize    Bin size in bp. Default to chromosome-wide bins.
  -p binStep    Bin step in bp. Default to bin sizeinStep.
"

# Default values
binSize=0
binStep=0
chrWide=true

# Parse options
while getopts hs:p:o:c: opt; do
	case $opt in
		h)
			# Help page
			echo -e "$helps"
			exit 0
		;;
		s)
			# Bin size
			if [ $OPTARG -le 0 ]; then
				msg="!!! ERROR! Invalid -s option. Bin size must be > 0."
				echo -e "$help\n$msg"
				exit 1
			else
				binSize=$OPTARG
				chrWide=false
			fi
		;;
		p)
			# Bin step
			if [ $OPTARG -le 0 ]; then
				msg="!!! ERROR! Invalid -p option. Bin step must be > 0."
				echo -e "$help\n$msg"
				exit 1
			else
				binStep=$OPTARG
			fi
		;;
		o)
			# Output directory
			if [ ! -d $OPTARG ]; then
				mkdir -p $OPTARG
			fi
			outdir=$OPTARG
		;;
		c)
			# Cutsite bedfile
			if [ ! -e $OPTARG ]; then
				msg="!!! ERROR! Invalid -c option. File not found: $OPTARG"
				echo -e "$help\n$msg"
				exit 1
			else
				csBed=$OPTARG
			fi
		;;
		?)
			msg="!!! ERROR! Unrecognized option."
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
if [ -z "$csBed" ]; then
	echo -e "$helps\n!!! ERROR! Missing mandatory -c option.\n"
	exit 1
fi
if [ ! -x "$(command -v bedtools)" -o -z "$(command -v bedtools)" ]; then
	echo -e "$helps\n!!! ERROR! Missing bedtools.\n"
	exit 1
fi

# Read bedfile paths
shift $(($OPTIND - 1))
bedfiles=()
for bf in $*; do
	if [ -e $bf -a -n $bf ]; then
		bedfiles+=("$bf")
	else
		msg="!!! Invalid bedfile, file not found.\n    File: $bf"
		echo -e " $helps\n$msg"
		exit 1
	fi
done
if [ 0 -eq ${#bedfiles[@]} ]; then
	msg="!!! No bedfile was specified!\n"
	echo -e " $helps\n$msg"
	exit 1
fi

# Additional checks
if [ ! $binStep -eq 0 -a $binSize -eq 0 ]; then
	echo -e "WARNING: missing bin size, ignoring -p option.\n"
fi
if [ ! $binSize -eq 0 -a $binStep -eq 0 ]; then
	binStep=$binSize
fi

# Print settings

settings=""
if $chrWide; then
	settings="$settings
 Using chr-wide bins."
else
	settings="$settings
   Bin size : $binSize
   Bin step : $binStep"
fi
settings="$settings
 
 Output dir : $outdir
   Cutsites : $csBed
  Bed files :
   $(echo ${bedfiles[@]} | sed 's/ /\n   /g')"

echo -e "$settings\n"

# RUN ==========================================================================

# 0) Identify chromosome sizes -------------------------------------------------
echo -e " Retrieving chromosome sizes ..."
chrSize=$(cat ${bedfiles[@]} | grep -v 'track' | datamash -sg1 -t$'\t' max 3)

# Sort chromosomes
echo -e "$chrSize" | awk -f "add_chr_id.awk" | sort -k1,1n | cut -f2,3 \
	> "$outdir/chr_size.tsv"


# 1) Generate bin bed file -----------------------------------------------------
echo -e " Generating bins ..."

prefix="bins.size$binSize.step$binStep"
if $chrWide; then
	cat "$outdir/chr_size.tsv" | awk "{ print $1 '\t' 0 '\t' $2 }" \
		> "$outdir/bins.size$binSize.step$binStep.bed"
else
	cat "$outdir/chr_size.tsv" | \
		awk -v size=$binSize -v step=$binStep -f "mk_bins.awk" \
		> "$outdir/$prefix.bed"
fi

# 2) Intersect with bedtools ---------------------------------------------------
echo -e " Intersecting ..."

# Assign bed reads to bins
for bfi in $(seq 0 $(bc <<< "${#bedfiles[@]} - 1")); do
	fname=$(echo -e "${bedfiles[$bfi]}" | tr "/" "\t" | awk '{ print $NF }')
	echo -e " > Intersecting $fname ..."
	bedtools intersect -a "$outdir/$prefix.bed" \
		 -b "${bedfiles[$bfi]}" -wa -wb | cut -f 1-3,8\
		> "$outdir/intersected.$prefix.$fname.tsv"
done

# Assign cutsites to bins
echo -e " > Intersecting $csBed ..."
bedtools intersect -a "$outdir/$prefix.bed" -b "$csBed" -c \
	> "$outdir/intersected.$prefix.cutsites.tsv"

# 3) Calculate bin statistics --------------------------------------------------
echo -e " Calculating bin statistics ..."

# Stats of cutsites
echo -e " > Calculating for $fname ..."
cat "$outdir/intersected.$prefix.cutsites.tsv" | \
	datamash -sg1,2,3 sum 4 | awk -f "add_chr_id.awk" | sort -k1,1n -k3,3n | \
	cut -f2- > "$outdir/bin_stats.$prefix.cutsites.tsv"

# Stats of beds
for bfi in $(seq 0 $(bc <<< "${#bedfiles[@]} - 1")); do
	fname=$(echo -e "${bedfiles[$bfi]}" | tr "/" "\t" | awk '{ print $NF }')
	binned="$outdir/intersected.$prefix.$fname.tsv"

	# Calculate statistics
	echo -e " > Calculating for $fname ..."
	bin_stats=$(cat "$binned" | datamash -sg1,2,3 sum 4 mean 4 svar 4 | \
		awk -f "add_chr_id.awk" | sort -k1,1n -k3,3n | cut -f2-)

	# Add number of cutsites
	awk -f "add_cs_sum.awk" <(cat "$outdir/bin_stats.$prefix.cutsites.tsv") \
		<(echo -e "$bin_stats") | cut -f 1-6,10 \
		> "$outdir/bin_stats.$prefix.$fname.tsv"
	rm "$binned"
done


# END ==========================================================================

################################################################################
