#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 1.0.0
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

# TEST =========================================================================

if [ 0 -eq 1 ]; then
	bedfiles=()
	bedfiles+=("#\nchr1\t0\t0\t0\t1\nchr2\t0\t0\t0\t2\nchr3\t0\t0\t0\t1")
	bedfiles+=("#\nchr1\t0\t0\t0\t1\nchr2\t0\t0\t0\t4\nchr3\t0\t0\t0\t3")
	bedfiles+=("#\nchr1\t0\t0\t0\t1\nchr2\t0\t0\t0\t8\nchr3\t0\t0\t0\t9")

	outFile="output"

	for i in $(seq 0 `bc <<< "${#bedfiles[@]} - 1"`); do
		echo -e ${bedfiles[$i]} > c$i.tmp
		bedfiles[$i]="c$i.tmp"
		echo -e "\n"${bedfiles[$i]}
		cat ${bedfiles[$i]}
	done
fi

# RUN ==========================================================================

# Generate cumulative probability distribution matrix --------------------------
# One chromosome per row, one condition per column.
echo -e " · Preparing matrices..."

# Default empty matrix
matrix_crs=""
matrix_rcs=""

# Cycle over chromosomes
for chr in $(echo $(seq 1 22) X); do
	chr="chr$chr"
	echo -e " >>> Working on $chr..."

	# Number of cutsite in the chromosome
	ncs=`cat $csList | awk -v chr=$chr '$1 == chr' | wc -l`

	if [ 0 -eq $ncs ]; then
		msg="!!! No cutsites found in $chr.\n    Please re-run without $chr."
		echo -e "$msg"
		exit 1
	fi

	# Array of normalized counts
	# crs: Cumulative of the RatioS
	# rcs: Ratio of the CumulativeS
	crs=(0)
	rcs=(0)
	counts=(0)
	totals=(0)
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
		n=`cat $bf | sed 1d | awk -v chr=$chr '$1 == chr' | \
			cut -f 5 | paste -sd+ | bc`
		
		if [ -z "$n" ]; then
			msg="!!! No reads found in condition #$i on $chr.\n"
			msg="$msg    Please re-run without the following bedFile: $bf\n"
			msg="$msg    Or without $chr."
			echo -e "$msg"
			exit 1
		fi

		# Ratio of the sums ----------------------------------------------------

		# Prepare rations
		prev_n=${counts[${#counts[@]}-1]}
		prev_ncc=${totals[${#totals[@]}-1]}
		r1=`bc -l <<< "($n + $prev_n) / ($ncc + $prev_ncc) / $ncs"`

		# Add leading 0
		r1=`echo "$r1" | sed -r "s/^\.(.*)$/0.\1/"`

		# To check the calculations uncomment the following line
		# echo -e "($n + $prev_n) / ($ncc + $prev_ncc) / $ncs = $r1"

		# Store in the arrays
		counts+=(`bc <<< "$n + $prev_n"`)
		totals+=(`bc <<< "$ncc + $prev_ncc"`)
		rcs+=($r1)

		# Sum of the ratios ----------------------------------------------------
		
		# Normalize
		r2=`bc -l <<< "$n / $ncc / $ncs"`
		
		# Sum previous ratio
		r2=`bc <<< "${crs[${#crs[@]} - 1]}+$r2"`

		# Add leading 0
		r2=`echo "$r2" | sed -r "s/^\.(.*)$/0.\1/"`

		# To check the calculations uncomment the following line
		# echo -e "$n / $ncc / $ncs = $r2"
		
		# Save normalized counts
		crs+=($r2)
	done

	# Remove starting value (0)
	unset crs[0]
	unset rcs[0]

	# Merging cumulative set
	merged_crs=`join_by " " "${crs[@]}"`
	merged_rcs=`join_by " " "${rcs[@]}"`

	# Forming new matrix row
	nrow_crs=`echo "$chr $merged_crs" | sed "s/^ //" | tr -s " " | tr " " "\t"`
	nrow_rcs=`echo "$chr $merged_rcs" | sed "s/^ //" | tr -s " " | tr " " "\t"`

	# Adding row to matrix
	matrix_crs="$matrix_crs$nrow_crs\n"
	matrix_rcs="$matrix_rcs$nrow_rcs\n"
done

# To save the matrix to file uncomment the following line
#echo -e "$matrix_crs" > matrix.crs.tmp.tsv
#echo -e "$matrix_rcs" > matrix.rcs.tmp.tsv

# Calculate ratio (B/A) between consecutive conditions -------------------------
echo -e " · Normalizing matrices..."

function normatrix() {
	matrix=$1

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
	echo -e "$matrix" | awk "$awkprogram"
}
normatrix_crs=`normatrix "$matrix_crs"`
normatrix_rcs=`normatrix "$matrix_rcs"`

# To save the normalized matrix to file uncomment the following line
#echo -e "$normatrix_crs" > normatrix.crs.tmp.tsv
#echo -e "$normatrix_rcs" > normatrix.rcs.tmp.tsv

# Sort & rank ------------------------------------------------------------------
echo -e " · Sorting and ranking..."

function rank_normatrix() {
	normatrix=$1
	nbeds=$2

	ranked=""
	for i in $(seq 2 $nbeds); do
		ranking=`echo -e "$normatrix" | cut -f 1,$i | \
			awk '{ print NR OFS $1 OFS $2 }' | sort -k3,3n | \
			awk '{ print $1 OFS $2 OFS NR }' | sort -k1,1n | cut -d " " -f 2,3`

		if [ -z "$ranked" ]; then
			ranked=$ranking
		else
			ranked=`join --nocheck-order <(echo -e "$ranked") \
				<(echo -e "$ranking")`
		fi
	done

	echo -e "$ranked"
}
ranked_crs=`rank_normatrix "$normatrix_crs" ${#bedfiles[@]}`
ranked_rcs=`rank_normatrix "$normatrix_rcs" ${#bedfiles[@]}`


# To save the conditional rankings to file uncomment the following line
#echo -e "$ranked_crs" > rankings.crs.tmp.tsv
#echo -e "$ranked_rcs" > rankings.rcs.tmp.tsv

# Sum rankings & sort ----------------------------------------------------------
echo -e " · Summing and sorting..."

function sumsort_ranks() {
	ranked=$1

	awkprogram='
	{
		tot=0;
		for ( i = 2; i <= NF; i++ )
			tot+=$i;
		print $1 OFS tot;
	}'
	echo -e "$ranked" | tr " " "\t" | awk "$awkprogram" | sort -k2,2n
}
sumsort_ranks "$ranked_crs"  > "$outFile".crs.txt
sumsort_ranks "$ranked_rcs"  > "$outFile".rcs.txt

# End --------------------------------------------------------------------------

################################################################################
