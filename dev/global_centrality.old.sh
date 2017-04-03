#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 1.0.1
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
 usage: ./global_centrality.sh [-hrdf] -c csList -o outFile [BEDFILE]...

 Description:
  Calculate global centrality metrics.

 Mandatory arguments:
  BEDFILE	Bed file(s). Expected to be ordered per condition.
  -c csList	Cutsite list file.
  -o outFile	Output matrix file.

 Optional arguments:
  -h		Show this help page.
  -r		Perform intermediate ranking.
  -d		Debug mode: write out intermediate results.
  -f		Calculate global centrality relative to first condition (fixed).
"

# Default options
interRank=false
debug=false
fixed=false

# Parse options
while getopts hdrfc:o: opt "${bedfiles[@]}"; do
	case $opt in
		h)
			echo -e "$helps\n"
			exit 0
		;;
		r)
			interRank=true
		;;
		d)
			debug=true
		;;
		f)
			fixed=true
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
	# ffs: Fano FactorS
	# cos: Coefficients Of Variation
	crs=(0)
	rcs=(0)
	ffs=()
	cos=()
	counts=(0)
	totals=(0)
	for i in $(seq 0 `expr ${#bedfiles[@]} - 1`); do
		bf=${bedfiles[$i]}

		# Select only current chromosome rows
		bfchr=`cat $bf | sed 1d | awk -v chr=$chr '$1 == chr'`
		
		# Number of reads in the condition
		ncc=`cat $bf | sed 1d | cut -f 5 | paste -sd+ | bc`
		
		if [ 0 -eq $ncc ]; then
			msg="!!! No reads found in condition #$i.\n"
			msg="$msg    Please re-run without the following bedFile: $bf"
			echo -e "$msg"
			exit 1
		fi

		# Number of reads in the condition in the chromosome
		n=`echo "$bfchr" | cut -f 5 | paste -sd+ | bc`
		
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

		# Variability based ---------------------------------------------------
		
		# Calculate mean
		mu=`echo "$bfchr" | datamash mean 5`
		
		# Calculate variance
		sd=`echo "$bfchr" | datamash sstdev 5`
		sigma=`bc -l <<< "$sd^2"`

		# Fano factor
		ffs+=(`bc -l <<< "$sigma^2 / $mu"`)

		# Coefficient of variation
		cov+=(`bc -l <<< "$sigma / $mu"`)		
	done

	# Remove starting value (0)
	unset crs[0]
	unset rcs[0]

	# Merging cumulative set
	merged_crs=`join_by " " "${crs[@]}"`
	merged_rcs=`join_by " " "${rcs[@]}"`
	merged_ffs=`join_by " " "${ffs[@]}"`
	merged_cov=`join_by " " "${cov[@]}"`

	# Forming new matrix row
	nrow_crs=`echo "$chr $merged_crs" | sed "s/^ //" | tr -s " " | tr " " "\t"`
	nrow_rcs=`echo "$chr $merged_rcs" | sed "s/^ //" | tr -s " " | tr " " "\t"`
	nrow_ffs=`echo "$chr $merged_ffs" | sed "s/^ //" | tr -s " " | tr " " "\t"`
	nrow_cov=`echo "$chr $merged_cov" | sed "s/^ //" | tr -s " " | tr " " "\t"`

	# Adding row to matrix
	matrix_crs="$matrix_crs$nrow_crs\n"
	matrix_rcs="$matrix_rcs$nrow_rcs\n"
	matrix_ffs="$matrix_ffs$nrow_ffs\n"
	matrix_cov="$matrix_cov$nrow_cov\n"
done

if $debug; then
	echo -e "$matrix_crs" > $outFile".matrix.crs.tmp.tsv"
	echo -e "$matrix_rcs" > $outFile".matrix.rcs.tmp.tsv"
	echo -e "$matrix_ffs" > $outFile".matrix.ffs.tmp.tsv"
	echo -e "$matrix_cov" > $outFile".matrix.cov.tmp.tsv"
fi

# Calculate ratio (B/A) between consecutive conditions -------------------------
echo -e " · Normalizing matrices..."

function normatrix_prev() {
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
function normatrix_first() {
	matrix=$1
	awkprogram='@include "join";
	NF {
		c=1;
	    for (B = 3; B <= NF; B++) {
	        a[c]=$B/$2;
	        c++;
	    }

		OFS=FS="\t";
		OFMT="%.4f"
	    print $1 OFS join(a, 1, c, OFS);
	}'
	echo -e "$matrix" | awk "$awkprogram"
}
function normatrix_2points() {
	matrix=$1
	awkprogram='@include "join";
	NF {
		OFS=FS="\t";
		OFMT="%.4f"
	    print $1 OFS $NF/$2;
	}'
	echo -e "$matrix" | awk "$awkprogram"
}
function difmatrix_prev() {
	matrix=$1
	awkprogram='@include "join";
	NF {
		c=1;
	    for (B = 3; B <= NF; B++) {
	        A=B-1;
	        a[c]=$B - $A;
	        c++;
	    }

		OFS=FS="\t";
		OFMT="%.4f"
	    print $1 OFS join(a, 1, c, OFS);
	}'
	echo -e "$matrix" | awk "$awkprogram"
}
function difmatrix_first() {
	matrix=$1
	awkprogram='@include "join";
	NF {
		c=1;
	    for (B = 3; B <= NF; B++) {
	        a[c]=$B - $2;
	        c++;
	    }

		OFS=FS="\t";
		OFMT="%.4f"
	    print $1 OFS join(a, 1, c, OFS);
	}'
	echo -e "$matrix" | awk "$awkprogram"
}
function difmatrix_2points() {
	matrix=$1
	awkprogram='@include "join";
	NF {
		OFS=FS="\t";
		OFMT="%.4f"
	    print $1 OFS $NF-$2;;
	}'
	echo -e "$matrix" | awk "$awkprogram"
}
if $fixed; then
	normatrix_crs=`normatrix_first "$matrix_crs"`
	normatrix_rcs=`normatrix_first "$matrix_rcs"`
	difmatrix_ffs=`difmatrix_first "$matrix_ffs"`
	difmatrix_cov=`difmatrix_first "$matrix_cov"`
else
	normatrix_crs=`normatrix_prev "$matrix_crs"`
	normatrix_rcs=`normatrix_prev "$matrix_rcs"`
	difmatrix_ffs=`difmatrix_prev "$matrix_ffs"`
	difmatrix_cov=`difmatrix_prev "$matrix_cov"`
fi
normatrix_crs_2p=`normatrix_2points "$matrix_crs"`
normatrix_rcs_2p=`normatrix_2points "$matrix_rcs"`
difmatrix_ffs_2p=`difmatrix_2points "$matrix_ffs"`
difmatrix_cov_2p=`difmatrix_2points "$matrix_cov"`

if $debug; then
	echo -e "$normatrix_crs" > $outFile".normatrix.crs.tmp.tsv"
	echo -e "$normatrix_rcs" > $outFile".normatrix.rcs.tmp.tsv"
	echo -e "$normatrix_crs_2p" > $outFile".normatrix.crs.2p.tmp.tsv"
	echo -e "$normatrix_rcs_2p" > $outFile".normatrix.rcs.2p.tmp.tsv"
	echo -e "$difmatrix_ffs" > $outFile".difmatrix.ffs.tmp.tsv"
	echo -e "$difmatrix_cov" > $outFile".difmatrix.cov.tmp.tsv"
	echo -e "$difmatrix_ffs_2p" > $outFile".difmatrix.ffs.2p.tmp.tsv"
	echo -e "$difmatrix_cov_2p" > $outFile".difmatrix.cov.2p.tmp.tsv"
fi

# Sort & rank ------------------------------------------------------------------
if $interRank; then
	echo -e " · Sorting and ranking..."

	function rank_normatrix() {
		normatrix=$1
		nbeds=$2

		ranked=""
		for i in $(seq 2 $nbeds); do
			ranking=`echo -e "$normatrix" | cut -f 1,$i | \
				awk '{ print NR OFS $1 OFS $2 }' | sort -k3,3n | \
				awk '{ print $1 OFS $2 OFS NR }' | sort -k1,1n | \
				cut -d " " -f 2,3`

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
	ranked_ffs=`rank_normatrix "$difmatrix_ffs" ${#bedfiles[@]}`
	ranked_cov=`rank_normatrix "$difmatrix_cov" ${#bedfiles[@]}`

	if $debug; then
		echo -e "$ranked_crs" > $outFile".rankings.crs.tmp.tsv"
		echo -e "$ranked_rcs" > $outFile".rankings.rcs.tmp.tsv"
		echo -e "$ranked_ffs" > $outFile".rankings.ffs.tmp.tsv"
		echo -e "$ranked_cov" > $outFile".rankings.cov.tmp.tsv"
	fi
else
	ranked_crs="$normatrix_crs"
	ranked_rcs="$normatrix_rcs"
	ranked_ffs="$difmatrix_ffs"
	ranked_cov="$difmatrix_cov"
fi

# Sum rankings & sort ----------------------------------------------------------
echo -e " · Summing and sorting..."

function sumsort_ranks() {
	ranked=$1

	awkprogram='
	{
		OFS=FS="\t";
		tot=0;
		for ( i = 2; i <= NF; i++ )
			tot+=$i;
		print $1 OFS tot;
	}'
	echo -e "$ranked" | tr " " "\t" | awk "$awkprogram" | sort -k2,2n
}
final_ranked_crs=`sumsort_ranks "$ranked_crs"`
final_ranked_rcs=`sumsort_ranks "$ranked_rcs"`
final_normatrix_crs_2p=`sumsort_ranks "$normatrix_crs_2p"`
final_normatrix_rcs_2p=`sumsort_ranks "$normatrix_rcs_2p"`
final_ranked_ffs=`sumsort_ranks "$ranked_ffs"`
final_ranked_cov=`sumsort_ranks "$ranked_cov"`
final_difmatrix_ffs_2p=`sumsort_ranks "$difmatrix_ffs_2p"`
final_difmatrix_cov_2p=`sumsort_ranks "$difmatrix_cov_2p"`

if $debug; then
	echo -e "$final_ranked_crs" > $outFile".crs.tmp.txt"
	echo -e "$final_ranked_rcs" > $outFile".rcs.tmp.txt"
	echo -e "$final_normatrix_crs_2p" > $outFile".crs.2p.tmp.txt"
	echo -e "$final_normatrix_rcs_2p" > $outFile".rcs.2p.tmp.txt"
	echo -e "$final_ranked_ffs" > $outFile".ffs.tmp.txt"
	echo -e "$final_ranked_cov" > $outFile".cov.tmp.txt"
	echo -e "$final_difmatrix_ffs_2p" > $outFile".ffs.2p.tmp.txt"
	echo -e "$final_difmatrix_cov_2p" > $outFile".cov.2p.tmp.txt"
fi

# Recap ------------------------------------------------------------------------
echo -e " · Recapping..."

# Print header
echo -e "CRS\tCRS2P\tRCS\tRCS2P\tFF\tFF2P\tCV\tCV2P" > $outFile".recap.txt"

# Merge and print rankings
paste \
	<(echo -e "$final_ranked_crs" | cut -f 1) \
	<(echo -e "$final_normatrix_crs_2p" | cut -f 1) \
	<(echo -e "$final_ranked_rcs" | cut -f 1) \
	<(echo -e "$final_normatrix_rcs_2p" | cut -f 1) \
	<(echo -e "$final_ranked_ffs" | cut -f 1) \
	<(echo -e "$final_difmatrix_ffs_2p" | cut -f 1) \
	<(echo -e "$final_ranked_cov" | cut -f 1) \
	<(echo -e "$final_difmatrix_cov_2p" | cut -f 1) \
	>> $outFile".recap.txt"

# End --------------------------------------------------------------------------

echo -e "\n        ~ FIN ~        "
echo -e "   └[∵┌]└[ ∵ ]┘[┐∵]┘   "

################################################################################
