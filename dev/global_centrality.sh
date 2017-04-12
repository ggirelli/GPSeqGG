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
 usage: ./global_centrality.sh	[-hidf][-r regfile]
				-c csList -o outFile [BEDFILE]...

 Description:
  Calculate global centrality metrics.

 Mandatory arguments:
  BEDFILE	Bed file(s). Expected to be ordered per condition.
  -c csList	Cutsite list file.
  -o outFile	Output matrix file.

 Optional arguments:
  -h		Show this help page.
  -i		Perform intermediate ranking.
  -d		Debug mode: write out intermediate results.
  -r regfile	Path to bedfile, containing regions to be assigned to.
"

# Default options
interRank=false
debug=false
fixed=false

# Parse options
while getopts hdifr:c:o: opt "${bedfiles[@]}"; do
	case $opt in
		h)
			echo -e "$helps\n"
			exit 0
		;;
		i)
			interRank=true
		;;
		d)
			debug=true
		;;
		r)
			if [ -e $OPTARG ]; then
				regFile=$OPTARG
			else
				msg="!!! Invalid -r option, file not found.\n    File: $OPTARG"
				echo -e "$helps\n$msg"
				exit 1
			fi
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
		echo -e " $helps\n$msg"
		exit 1
	fi
done

# RUN ==========================================================================

# Assign to regions ------------------------------------------------------------
# If a list of regions is specified, assign them to the reads
# Then, treat them as "chromosome" and run the analysis on them

if [ -n "$regFile" ]; then
	echo -e " · Assigning to regions of interest..."

	rois2chr='{
		OFS=FS="\t"
		nrois = split($6, rois, " ");
		for (i = 1; i <= nrois; i++) {
			print rois[i] OFS $2 OFS $3 OFS $4 OFS $5 OFS
		}
	}'

	# Make temporary bedfiles
	for bfi in $(seq 0 `bc <<< "${#bedfiles[@]} - 1"`); do
		bf=${bedfiles[$bfi]}

		# # Assign bedfile rows to ROIs
		# echo -e " >> Working on $bf ..."

		# # Use ROIs as chromosomes
		# echo -e " >>> Using regions as chromosomes..."
		# ./bed_addROIs.py $regFile $bf | awk "$rois2chr" \
		# 	> $bf".regions.tmp"

		# Use tmp file with fake chromosomes
		bedfiles[$bfi]=$bf".regions.tmp"
	done

	# Make temporary cutsite list
	echo -e " >> Working on $csList ..."
	# ./bed_addROIs.py $regFile $csList | awk "$rois2chr" \
	# 	> $csList".regions.tmp"
	csList=$csList".regions.tmp"
fi

# Identify regions -------------------------------------------------------------
echo -e " · Identifying regions..."

chr_list=""
for bf in ${bedfiles[@]}; do
	chr_list=$chr_list" "`cat $bf | sed 1d | cut -f 1 | sort | uniq`
done
chr_list=`echo $chr_list | tr ' ' '\n' | sort | uniq | tr '\n' ' '`
IFS=' ' read -r -a chr_list <<< "$chr_list"

# Generate cumulative probability distribution matrix --------------------------
# One chromosome per row, one condition per column.
echo -e " · Preparing matrices..."

# Default empty matrix
matrix_crs=""
matrix_rcs=""

# Cycle over chromosomes
for chr in ${chr_list[@]}; do
	echo -e " >>> Working on $chr..."

	# Number of cutsite in the chromosome
	ncs=`cat $csList | awk -v chr=$chr '$1 == chr' | wc -l`

	if [ 0 -eq $ncs ]; then
		msg="!!! No cutsites found in $chr.\n    Please re-run without $chr."
		echo -e "$msg"
		exit 1
	fi

	# Array of normalized counts
	# prs: ProbabilitieS
	# crs: Cumulative of the RatioS
	# rcs: Ratio of the CumulativeS
	# std: STandard Deviation
	# ffs: Fano FactorS
	# cov: Coefficients Of Variation
	prs=()
	crs=(0)
	rcs=(0)
	std=()
	ffs=()
	cov=()
	counts=(0)
	totals=(0)
	for i in $(seq 0 `expr ${#bedfiles[@]} - 1`); do
		bf=${bedfiles[$i]}

		# Select only current chromosome rows
		bfchr=`cat $bf | sed 1d | awk -v chr=$chr '$1 == chr'`
		
		# Number of reads in the condition
		tot=`cat $bf | sed 1d | cut -f 5 | paste -sd+ | bc`
		
		if [ 0 -eq $tot ]; then
			msg="!!! No reads found in condition #$i.\n"
			msg="$msg    Please re-run without the following bedFile: $bf"
			echo -e "$msg"
			exit 1
		fi

		# Number of reads in the condition in the chromosome
		n=`echo "$bfchr" | cut -f 5 | paste -sd+ | bc`
		if [ -z "$n" ]; then
			n=0
		fi

		# if [ -z "$n" ]; then
		# 	msg="!!! No reads found in condition #$i on $chr.\n"
		# 	msg="$msg    Please re-run without the following bedFile: $bf\n"
		# 	msg="$msg    Or without $chr."
		# 	echo -e "$msg"
		# 	exit 1
		# fi
		
		# Probability ----------------------------------------------------------

		r=`bc -l <<< "$n / $tot / $ncs"`
		prs+=($r)

		# Ratio of the sums ----------------------------------------------------

		# Prepare ratios
		prev_n=${counts[${#counts[@]}-1]}
		prev_t=${totals[${#totals[@]}-1]}
		r=`bc -l <<< "($n + $prev_n) / ($tot + $prev_t) / $ncs"`

		# Add leading 0
		r=`echo "$r" | sed -r "s/^\.(.*)$/0.\1/"`

		# To check the calculations uncomment the following line
		# echo -e "($n + $prev_n) / ($tot + $prev_t) / $ncs = $r"

		# Store in the arrays
		counts+=(`bc <<< "$n + $prev_n"`)
		totals+=(`bc <<< "$tot + $prev_t"`)
		rcs+=($r)

		# Sum of the ratios ----------------------------------------------------
		
		# Normalize
		r=`bc -l <<< "$n / $tot / $ncs"`
		
		# Sum previous ratio
		r=`bc <<< "${crs[${#crs[@]} - 1]}+$r"`

		# Add leading 0
		r=`echo "$r" | sed -r "s/^\.(.*)$/0.\1/"`

		# To check the calculations uncomment the following line
		# echo -e "$n / $tot / $ncs = $r"
		
		# Save normalized counts
		crs+=($r)

		# Variability based ----------------------------------------------------
		
		if [ 0 == $n ]; then
			# Calculate variance
			std+=("nan")

			# Fano factor
			ffs+=("nan")

			# Coefficient of variation
			cov+=("nan")
		else
			# Calculate mean
			mu=`echo "$bfchr" | datamash mean 5`
			
			# Calculate variance
			sd=`echo "$bfchr" | datamash sstdev 5`
			std+=($sd)

			if [ 0 == $sd ]; then
				echo -e "$bfchr" | cut -f 5
			fi

			sigma=`bc -l <<< "$sd^2"`

			# Fano factor
			ffs+=(`bc -l <<< "$sigma^2 / $mu"`)

			# Coefficient of variation
			cov+=(`bc -l <<< "$sigma / $mu"`)
		fi	
	done

	# Remove starting value (0)
	unset crs[0]
	unset rcs[0]

	# Merging cumulative set
	merged_prs=`join_by " " "${prs[@]}"`
	merged_crs=`join_by " " "${crs[@]}"`
	merged_rcs=`join_by " " "${rcs[@]}"`
	merged_std=`join_by " " "${std[@]}"`
	merged_ffs=`join_by " " "${ffs[@]}"`
	merged_cov=`join_by " " "${cov[@]}"`

	# Forming new matrix row
	nrow_prs=`echo "$chr $merged_prs" | sed "s/^ //" | tr -s " " | tr " " "\t"`
	nrow_crs=`echo "$chr $merged_crs" | sed "s/^ //" | tr -s " " | tr " " "\t"`
	nrow_rcs=`echo "$chr $merged_rcs" | sed "s/^ //" | tr -s " " | tr " " "\t"`
	nrow_std=`echo "$chr $merged_std" | sed "s/^ //" | tr -s " " | tr " " "\t"`
	nrow_ffs=`echo "$chr $merged_ffs" | sed "s/^ //" | tr -s " " | tr " " "\t"`
	nrow_cov=`echo "$chr $merged_cov" | sed "s/^ //" | tr -s " " | tr " " "\t"`

	# Adding row to matrix
	matrix_prs="$matrix_prs$nrow_prs\n"
	matrix_crs="$matrix_crs$nrow_crs\n"
	matrix_rcs="$matrix_rcs$nrow_rcs\n"
	matrix_std="$matrix_std$nrow_std\n"
	matrix_ffs="$matrix_ffs$nrow_ffs\n"
	matrix_cov="$matrix_cov$nrow_cov\n"
done

if $debug; then
	echo -e "$matrix_prs" > $outFile".prs.matrix.tmp.tsv"
	echo -e "$matrix_crs" > $outFile".crs.matrix.tmp.tsv"
	echo -e "$matrix_rcs" > $outFile".rcs.matrix.tmp.tsv"
	echo -e "$matrix_std" > $outFile".std.matrix.tmp.tsv"
	echo -e "$matrix_ffs" > $outFile".ffs.matrix.tmp.tsv"
	echo -e "$matrix_cov" > $outFile".cov.matrix.tmp.tsv"
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
		print $1 OFS $NF-$2;
	}'
	echo -e "$matrix" | awk "$awkprogram"
}
function ratmatrix_2points() {
	matrix=$1
	awkprogram='@include "join";
	NF {
		OFS=FS="\t";
		OFMT="%.4f"
		if ( "nan" == $2 ) {
			print $1 OFS "nan";
		} else {
			print $1 OFS $NF/$2;
		}
	}'
	echo -e "$matrix" | awk "$awkprogram"
}
function ratmatrix_first() {
	matrix=$1
	awkprogram='@include "join";
	NF {
		c=1;
		for (B = 3; B <= NF; B++) {
			if ( "nan" == $2 ) {
				a[c]="nan"
			} else {
				a[c]=$B / $2;
			}
			c++;
		}

		OFS=FS="\t";
		OFMT="%.4f"
		print $1 OFS join(a, 1, c, OFS);
	}'
	echo -e "$matrix" | awk "$awkprogram"
}

normatrix_prs_2p=`normatrix_2points "$matrix_prs"`
normatrix_prs=`normatrix_prev "$matrix_prs"`
normatrix_prs_fixed=`normatrix_first "$matrix_prs"`

normatrix_crs=`normatrix_prev "$matrix_crs"`
normatrix_crs_fixed=`normatrix_first "$matrix_crs"`
normatrix_crs_2p=`normatrix_2points "$matrix_crs"`

normatrix_rcs_2p=`normatrix_2points "$matrix_rcs"`
normatrix_rcs_fixed=`normatrix_first "$matrix_rcs"`

ratmatrix_std_2p=`ratmatrix_2points "$matrix_std"`
ratmatrix_std_fixed=`ratmatrix_first "$matrix_std"`

difmatrix_ffs_2p=`difmatrix_2points "$matrix_ffs"`
difmatrix_ffs_fixed=`difmatrix_first "$matrix_ffs"`

difmatrix_cov_2p=`difmatrix_2points "$matrix_cov"`
difmatrix_cov_fixed=`difmatrix_first "$matrix_cov"`

if $debug; then
	echo -e "$normatrix_prs_2p" > $outFile".prs.normatrix.2p"
	echo -e "$normatrix_prs" > $outFile".prs.normatrix"
	echo -e "$normatrix_prs_fixed" > $outFile".prs.normatrix.fixed"
	echo -e "$normatrix_crs" > $outFile".crs.normatrix"
	echo -e "$normatrix_crs_fixed" > $outFile".crs.normatrix.fixed"
	echo -e "$normatrix_crs_2p" > $outFile".crs.normatrix.2p"
	echo -e "$normatrix_rcs_2p" > $outFile".rcs.normatrix.2p"
	echo -e "$normatrix_rcs_fixed" > $outFile".rcs.normatrix.fixed"
	echo -e "$ratmatrix_std_2p" > $outFile".std.difmatrix.2p"
	echo -e "$ratmatrix_std_fixed" > $outFile".std.difmatrix.fixed"
	echo -e "$difmatrix_ffs_2p" > $outFile".ffs.difmatrix.2p"
	echo -e "$difmatrix_ffs_fixed" > $outFile".ffs.difmatrix.fixed"
	echo -e "$difmatrix_cov_2p" > $outFile".cov.difmatrix.2p"
	echo -e "$difmatrix_cov_fixed" > $outFile".cov.difmatrix.fixed"
fi

# Sum rankings & sort ----------------------------------------------------------
echo -e " · Summing and sorting..."

function sumsort_ranks() {
	ranked=$1

	awkprogram='
	{
		OFS=FS="\t";
		tot=0;
		nancount=0;
		for ( i = 2; i <= NF; i++ ) {
			if ( "nan" == $i ) {
				nancount++;
			} else {
				tot+=$i;
			}
		}

		if ( 0 != nancount ) {
			print "nan" OFS "nan";
		} else {
			print $1 OFS tot;
		}
	}'
	echo -e "$ranked" | tr " " "\t" | awk "$awkprogram" | sort -k2,2n
}

ranked_prs_2p=`sumsort_ranks "$normatrix_prs_2p"`
ranked_prs=`sumsort_ranks "$normatrix_prs"`
ranked_prs_fixed=`sumsort_ranks "$normatrix_prs_fixed"`
ranked_crs=`sumsort_ranks "$normatrix_crs"`
ranked_crs_fixed=`sumsort_ranks "$normatrix_crs_fixed"`
ranked_crs_2p=`sumsort_ranks "$normatrix_crs_2p"`
ranked_rcs_2p=`sumsort_ranks "$normatrix_rcs_2p"`
ranked_rcs_fixed=`sumsort_ranks "$normatrix_rcs_fixed"`
ranked_std_2p=`sumsort_ranks "$ratmatrix_std_2p"`
ranked_std_fixed=`sumsort_ranks "$ratmatrix_std_fixed"`
ranked_ffs_2p=`sumsort_ranks "$difmatrix_ffs_2p"`
ranked_ffs_fixed=`sumsort_ranks "$difmatrix_ffs_fixed"`
ranked_cov_2p=`sumsort_ranks "$difmatrix_cov_2p"`
ranked_cov_fixed=`sumsort_ranks "$difmatrix_cov_fixed"`

if $debug; then
	echo -e "$ranked_prs_2p" > $outFile".prs.2p.tmp.txt"
	echo -e "$ranked_prs" > $outFile".prs.tmp.txt"
	echo -e "$ranked_prs_fixed" > $outFile".prs.fixed.tmp.txt"
	echo -e "$ranked_crs" > $outFile".crs.tmp.txt"
	echo -e "$ranked_crs_fixed" > $outFile".crs.fixed.tmp.txt"
	echo -e "$ranked_crs_2p" > $outFile".crs.2p.tmp.txt"
	echo -e "$ranked_rcs_2p" > $outFile".rcs.2p.tmp.txt"
	echo -e "$ranked_rcs_fixed" > $outFile".rcs.fixed.tmp.txt"
	echo -e "$ranked_std_2p" > $outFile".std.2p.tmp.txt"
	echo -e "$ranked_std_fixed" > $outFile".std.fixed.tmp.txt"
	echo -e "$ranked_ffs_2p" > $outFile".ffs.2p.tmp.txt"
	echo -e "$ranked_ffs_fixed" > $outFile".ffs.fixed.tmp.txt"
	echo -e "$ranked_cov_2p" > $outFile".cov.2p.tmp.txt"
	echo -e "$ranked_cov_fixed" > $outFile".cov.fixed.tmp.txt"
fi

# Recap ------------------------------------------------------------------------
echo -e " · Recapping..."

# Print header
header="P2p\tPg\tPgf\tCR2p\tCRg\tCRgf\tRC2p\tRCgf\tV2p\tVgf"
header=$header"\tFF2p\tFFgf\tCV2p\tCVgf"
echo -e $header > $outFile".recap.txt"

# Merge and print rankings
paste \
	<(echo -e "$ranked_prs_2p" | cut -f 1) \
	<(echo -e "$ranked_prs" | cut -f 1) \
	<(echo -e "$ranked_prs_fixed" | cut -f 1) \
	<(echo -e "$ranked_crs" | cut -f 1) \
	<(echo -e "$ranked_crs_fixed" | cut -f 1) \
	<(echo -e "$ranked_crs_2p" | cut -f 1) \
	<(echo -e "$ranked_rcs_2p" | cut -f 1) \
	<(echo -e "$ranked_rcs_fixed" | cut -f 1) \
	<(echo -e "$ranked_std_2p" | cut -f 1) \
	<(echo -e "$ranked_std_fixed" | cut -f 1) \
	<(echo -e "$ranked_ffs_2p" | cut -f 1) \
	<(echo -e "$ranked_ffs_fixed" | cut -f 1) \
	<(echo -e "$ranked_cov_2p" | cut -f 1) \
	<(echo -e "$ranked_cov_fixed" | cut -f 1) \
	>> $outFile".recap.txt"

# End --------------------------------------------------------------------------

echo -e "\n        ~ FIN ~        "
echo -e "   └[∵┌]└[ ∵ ]┘[┐∵]┘   "

################################################################################
