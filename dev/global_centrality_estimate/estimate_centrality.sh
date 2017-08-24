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
  Calculate global centrality metrics. Requires bedtools for bin assignment,
  datamash for calculations, and gawk for text manipulation.

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
if [ ! -x "$(command -v datamash)" -o -z "$(command -v datamash)" ]; then
	echo -e "$helps\n!!! ERROR! Missing datamash.\n"
	exit 1
fi
if [ ! -x "$(command -v gawk)" -o -z "$(command -v gawk)" ]; then
	echo -e "$helps\n!!! ERROR! Missing gawk.\n"
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

# FUNCTIONS ====================================================================

# Add human chromosome numeric ID as first column of a bed-like table
add_chr_id='
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

# Merge two bed files based on first three columns 
merge_beds='
	BEGIN {
		OFS = FS = "\t";
		sep = "~";
	}

	( FNR == NR ) {
		k = $1 sep $2 sep $3;
		a[k] = $0;
		next;
	}

	{
		k = $1 sep $2 sep $3;
		if ( k in a ) {
			print $0 OFS a[k];
		}
	}'

# Produce a bed of bins
mk_bins='BEGIN { OFS = FS = "\t"; }
	{ for ( i = 0; i < $2; i += step ) { print $1 OFS i OFS i+size; } }'

add_cnr_bfi='BEGIN { OFS = FS = "\t"; }
	{
		$4 = cnr OFS $4;
		print $0 OFS bfi;
	}'

# Estimate centrality, requires two input variables:
#  calc: 'ratio' or 'diff'
#  type: '2p' (two point), 'f' (fixed) or 'g' (global)
estimate_centrality='
	BEGIN {
		OFS = FS = "\t";

	}

	function estimate(calc, a, b) {
		switch (calc) {
			case "ratio":
				return a / b;
				break;
			case "diff":
				return a - b;
				break;
		}
	}

	{
		if ( cumrat == 1) {
			# Sum probabilities
			for ( i = 2; i <= NF; i++ ) {
				$i = $i + $(i-1);
			}
		}

		if ( ratcum == 1 ) {
			# Build table
			for ( i = 1; i <= NF; i++ ) {
				nf=split($i, ff, ",");
				for ( j = 1; j <= nf; j++ ) {
					a[i, j] = ff[j];
				}
			}
			# Calculate ratio of cumulatives
			for ( i = 1; i <= NF; i++ ) {
				a[i, 2] += a[i - 1, 2];
				a[i, 1] += a[i - 1, 1];
				$i = a[i, 2] / (a[i, 1] * a[i, 3]);
			}
		}

		switch (type) {
			case "2p":
				print estimate(calc, $NF, $1);
				break;
			case "f":
				output = 0;
				for ( i = 2; i <= NF; i++ ) {
					output = output + estimate(calc, $i, $1);
				}
				print output;
				break;
			case "g":
				output = 0;
				for ( i = 2; i <= NF; i++ ) {
					output = output + estimate(calc, $i, $(i-1));
				}
				print output;
				break;
		}
	}'


# RUN ==========================================================================

# 0) Identify chromosome sizes -------------------------------------------------
echo -e " Retrieving chromosome sizes ..."
chrSize=$(cat ${bedfiles[@]} | grep -v 'track' | datamash -sg1 -t$'\t' max 3)

# Sort chromosomes
echo -e "$chrSize" | gawk "$add_chr_id" | sort -k1,1n | cut -f2,3 \
	> "$outdir/chr_size.tsv"


# 1) Generate bin bed file -----------------------------------------------------
echo -e " Generating bins ..."

# Set output prefix
if $chrWide; then
	prefix="bins.chrWide"
else
	prefix="bins.size$binSize.step$binStep"
fi

# Generate bins
if $chrWide; then
	cat "$outdir/chr_size.tsv" | gawk '{ print $1 "\t" 0 "\t" $2 }' \
		> "$outdir/$prefix.bed" & pid=$!
else
	cat "$outdir/chr_size.tsv" | \
		gawk -v size=$binSize -v step=$binStep "$mk_bins" \
		> "$outdir/$prefix.bed" & pid=$!
fi
wait $pid; rm "$outdir/chr_size.tsv"

# 2) Intersect with bedtools ---------------------------------------------------
echo -e " Intersecting ..."

# Assign bed reads to bins
for bfi in $(seq 0 $(bc <<< "${#bedfiles[@]} - 1")); do
	fname=$(echo -e "${bedfiles[$bfi]}" | tr "/" "\t" | gawk '{ print $NF }')
	echo -e " > Intersecting $fname ..."
	bedtools intersect -a "$outdir/$prefix.bed" \
		 -b "${bedfiles[$bfi]}" -wa -wb | cut -f 1-3,8\
		> "$outdir/intersected.$prefix.$fname.tsv"
done

# Assign cutsites to bins
echo -e " > Intersecting $csBed ..."
bedtools intersect -a "$outdir/$prefix.bed" -b "$csBed" -c \
	> "$outdir/intersected.$prefix.cutsites.tsv" & pid=$!
wait $pid; rm "$outdir/$prefix.bed"
# 3) Calculate bin statistics --------------------------------------------------
echo -e " Calculating bin statistics ..."

# Stats of cutsites
echo -e " > Calculating for $csBed ..."
cat "$outdir/intersected.$prefix.cutsites.tsv" | \
	datamash -sg1,2,3 sum 4 | gawk "$add_chr_id" | sort -k1,1n -k3,3n | \
	cut -f2- > "$outdir/bin_stats.$prefix.cutsites.tsv" & pid=$!
wait $pid; rm "$outdir/intersected.$prefix.cutsites.tsv"

# Stats of beds
for bfi in $(seq 0 $(bc <<< "${#bedfiles[@]} - 1")); do
	fname=$(echo -e "${bedfiles[$bfi]}" | tr "/" "\t" | gawk '{ print $NF }')
	binned="$outdir/intersected.$prefix.$fname.tsv"

	# Calculate statistics
	echo -e " > Calculating for $fname ..."
	bin_stats=$(cat "$binned" | datamash -sg1,2,3 sum 4 mean 4 svar 4 | \
		gawk "$add_chr_id" | sort -k1,1n -k3,3n | cut -f2-)

	# Add number of cutsites
	gawk "$merge_beds" <(cat "$outdir/bin_stats.$prefix.cutsites.tsv") \
		<(echo -e "$bin_stats") | cut -f 1-6,10 \
		> "$outdir/bin_stats.$prefix.$fname.tsv" & pid=$!
	wait $pid; rm "$binned"
done
rm "$outdir/bin_stats.$prefix.cutsites.tsv"

# 4) Assemble into bin data table ----------------------------------------------
echo -e " Combining information ..."

# Normalize read count by cutsite and condition
norm=""
for bfi in $(seq 0 $(bc <<< "${#bedfiles[@]} - 1")); do
	fname=$(echo -e "${bedfiles[$bfi]}" | tr "/" "\t" | gawk '{ print $NF }')
	stats="$outdir/bin_stats.$prefix.$fname.tsv"

	# Normalize
	#echo -e " > Normalizing for $fname ..."
	cond_n_reads=$(cat "${bedfiles[$bfi]}" | grep -v "track" | datamash sum 5)
	norm="$norm"$(cat "$stats" | gawk -v cnr=$cond_n_reads -v bfi=$bfi \
		"$add_cnr_bfi")"\n"
	rm "$stats"
done
# Columns:
# 1   2     3   4        5       6      7         8   9
# chr|start|end|condRead|readSum|readMu|readSigma|nCS|condID
echo -e "$norm" | gawk "$add_chr_id" | sort -k1,1n -k3,3n -k10,10n | \
	cut -f2- | sed 1d > "$outdir/normalized.$prefix.tsv"

# 5) Estimate centrality -------------------------------------------------------
echo -e " Estimating centrality ..."

# Prepare paste string
spaste=""; for i in $(seq 1 ${#bedfiles[@]}); do spaste="$spaste -"; done

#---------------------#
# Probability metrics #
#---------------------#

# Probability metric
prob_mat=$(cut -f4,5,8 "$outdir/normalized.$prefix.tsv" | \
	gawk '{ print $2 / ($1 * $3) }' | paste $spaste)
probability_two_points=$(echo -e "$prob_mat" | \
	gawk -v calc="ratio" -v type="2p" "$estimate_centrality")
probability_fixed=$(echo -e "$prob_mat" | \
	gawk -v calc="ratio" -v type="f" "$estimate_centrality")
probability_global=$(echo -e "$prob_mat" | \
	gawk -v calc="ratio" -v type="g" "$estimate_centrality")

# Cumulative ratio metric
cumrat=$(cut -f4,5,8 "$outdir/normalized.$prefix.tsv" | \
	gawk '{ print $2 / ($1 * $3) }' | paste $spaste)
cumrat_two_points=$(echo -e "$cumrat" | \
	gawk -v calc="ratio" -v cumrat=1 -v type="2p" "$estimate_centrality")
cumrat_fixed=$(echo -e "$cumrat" | \
	gawk -v calc="ratio" -v cumrat=1 -v type="f" "$estimate_centrality")
cumrat_global=$(echo -e "$cumrat" | \
	gawk -v calc="ratio" -v cumrat=1 -v type="g" "$estimate_centrality")

# Ratio cumulative metric
ratcum=$(cut -f4,5,8 "$outdir/normalized.$prefix.tsv" | \
	tr '\t' ',' | paste $spaste)
ratcum_two_points=$(echo -e "$ratcum" | \
	gawk -v calc="ratio" -v ratcum=1 -v type="2p" "$estimate_centrality")
ratcum_fixed=$(echo -e "$ratcum" | \
	gawk -v calc="ratio" -v ratcum=1 -v type="f" "$estimate_centrality")
ratcum_global=$(echo -e "$ratcum" | \
	gawk -v calc="ratio" -v ratcum=1 -v type="g" "$estimate_centrality")

#---------------------#
# Variability metrics #
#---------------------#

# Variance metric

# Fano factor metric

# Coefficient of variation metric

#--------#
# Output #
#--------#

# Prepare output table
metrics=$(cut -f1-3 "$outdir/normalized.$prefix.tsv" | uniq | paste -d$'\t' - \
	<(echo -e "$probability_two_points") \
	<(echo -e "$probability_fixed") \
	<(echo -e "$probability_global") \
	<(echo -e "$cumrat_two_points") \
	<(echo -e "$cumrat_fixed") \
	<(echo -e "$cumrat_global") \
	<(echo -e "$ratcum_two_points") \
	<(echo -e "$ratcum_fixed") \
	<(echo -e "$ratcum_global") \
	)

# Remove bin positions if chromosome wide
if $chrWide; then
	metrics=$(echo -e "$metrics" | cut -f1,4-)
fi

# Write output
echo -e "$metrics" > "$outdir/estimates.$prefix.tsv"

# END ==========================================================================

################################################################################
