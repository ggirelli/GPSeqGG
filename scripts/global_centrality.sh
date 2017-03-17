#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 0.1.0
# Description: calculates global centrality metrics
# 
# ------------------------------------------------------------------------------



# ENV VAR ======================================================================

export LC_ALL=C

# DEPENDENCIES =================================================================

function join_by { local IFS="$1"; shift; echo "$*"; }

# INPUT ========================================================================

# Help string
helps="
 usage: ./global_centrality.sh [-h] -c csList -o outFile [BEDFILE]...

 Description:
  Calculate global centrality metrics.

 Mandatory arguments:
  BEDFILE	Bed file(s). Expected to be ordered per condition.
  -c csList	Cutsite list file.
  -o outFile	Output matrix file.

 Optional arguments:
  -h		Show this help page.
"

# Parse options
while getopts hc:o: opt "${bedfiles[@]}"; do
	case $opt in
		h)
			echo -e "$helps\n"
			exit 0
		;;
		c)
			if [ -e $OPTARG ]; then
				csList=$OPTARG
			else
				msg="!!! Invalid -c option, file not found.\n    File: $OPTARG"
				echo -e "$helps\n$msg"
				exit 1
			fi
		;;
		o)
			outFile=$OPTARG
		;;
	esac
done

# Check mandatory options
if [ -z "$csList" ]; then
	msg="!!! Missing mandatory -c option."
	echo -e "$helps\n$msg"
	exit
fi
if [ -z "$outFile" ]; then
	msg="!!! Missing mandatory -o option."
	echo -e "$helps\n$msg"
	exit
fi

# Read bedfile paths
shift $(($OPTIND - 1))
bedfiles=()
for bf in $*; do
	if [ -e $bf ]; then
		bedfiles+=("$bf")
	else
		msg="!!! Invalid bedfile, file not found.\n    File: $bf"
		echo -e " $helps<n$msg"
		exit 1
	fi
done

# RUN ==========================================================================

# Generate cumulative probability distribution matrix --------------------------
# One chromosome per row, one condition per column.

# Default empty matrix
matrix=""

# Cycle over chromosomes
for chr in $(echo $(seq 1 22) X); do
	chr="chr$chr"

	# Number of cutsite in the chromosome
	ncs=`cat $csList | grep $chr | wc -l`

	if [ 0 -eq $ncs ]; then
		msg="!!! No cutsites found in $chr.\n    Please re-run without $chr."
		echo -e "$msg"
		exit 1
	fi

	# Array of normalized counts
	counts=(0)
	for i in $(seq 0 `expr ${#bedfiles[@]} - 1`); do
		bf=${bedfiles[$i]}

		# Number of reads in the condition
		ncc=`cat $bf | sed 1d | cut -f 5 | paste -sd+ | bc`
		
		if [ 0 -eq $ncc ]; then
			msg="!!! No reads found in condition #$i.\n"
			msg="$msg    Please re-run without the following bedFile: $bf"
			echo -e "$msg"
			exit 1
		fi

		# Number of reads in the condition in the chromosome
		n=`cat $bf | sed 1d | grep $chr | cut -f 5 | paste -sd+ | bc`
		
		if [ -z "$n" ]; then
			msg="!!! No reads found in condition #$i on $chr.\n"
			msg="$msg    Please re-run without the following bedFile: $bf\n"
			msg="$msg    Or without $chr."
			echo -e "$msg"
			exit 1
		fi
		
		# Normalize
		r=`bc -l <<< "$n / $ncc / $ncs"`
		
		# Sum previous ratio
		r=`bc <<< "${counts[${#counts[@]} - 1]}+$r"`

		# Add leading 0
		r=`echo "$r" | sed -r "s/^\.(.*)$/0.\1/"`

		# To check the calculations uncomment the following line
		# echo -e "$n / $ncc / $ncs + ${counts[${#counts[@]} - 1]} = $r"
		
		# Save normalized counts
		counts+=($r)
	done

	# Remove starting value (0)
	unset counts[0]

	# Merging cumulative set
	merged=`join_by " " "${counts[@]}"`

	# Forming new matrix row
	newrow=`echo "$chr $merged" | sed "s/^ //" | tr -s " " | tr " " "\t"`

	# Adding row to matrix
	matrix="$matrix$newrow\n"
done

# To save the matrix to file uncomment the following line
# echo -e "$matrix" > matrix.tmp.tsv

# Calculate ratio (B/A) between consecutive conditions -------------------------

awkprogram='@include "join";
NF {
	c=1;
    for (B = 3; B <= NF; B++) {
        A=B-1;
        a[c]=$B/$A;
        c++;
    }

	OFS=FS="\t";
	OFMT="%.4f"
    print $1 OFS join(a, 1, c, OFS);
}'
normatrix=`echo -e "$matrix" | awk "$awkprogram"`

# To save the normalized matrix to file uncomment the following line
# echo -e "$normatrix" > normatrix.tmp.tsv

# Sort & rank ------------------------------------------------------------------

ranked=""
for i in $(seq 2 `expr ${#bedfiles[@]}`); do
	ranking=`echo -e "$normatrix" | cut -f 1,$i | \
		awk '{ print NR OFS $1 OFS $2 }' | sort -k3,3n | \
		awk '{ print $1 OFS $2 OFS NR }' | sort -k1,1n | cut -d " " -f 2,3`

	if [ -z "$ranked" ]; then
		ranked=$ranking
	else
		ranked=`join --nocheck-order <(echo -e "$ranked") <(echo -e "$ranking")`
	fi
done

# To save the conditional rankings to file uncomment the following line
# echo -e "$ranked" > rankings.tmp.tsv

# Sum rankings & sort ----------------------------------------------------------

awkprogram='
{
	tot=0;
	for ( i = 2; i <= NF; i++ )
		tot+=$i;
	print $1 OFS tot;
}'
echo -e "$ranked" | tr " " "\t" | awk "$awkprogram" | sort -k2,2n > $outFile

# End --------------------------------------------------------------------------

################################################################################
