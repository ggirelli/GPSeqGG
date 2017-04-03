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
IFS=', ' read -r -a array <<< "$string"

# INPUT ========================================================================

# Help string
helps="
 usage: ./global_centrality.sh [-hd] -c csList -i mergedBed -o outFile

 Description:
  Calculate global centrality metrics.

 Mandatory arguments:
  -c csList		Cutsite list file.
  -i mergedBed	Merged bed file(s). Expected to be ordered per condition.
  -o outFile	Output file prefix. (e.g., TK51.global_centrality)

 Optional arguments:
  -h	Show this help page.
  -d	Debug mode: write out intermediate results.
"

# Default options
debug=false

# Parse options
while getopts hdc:i:o: opt "${bedfiles[@]}"; do
	case $opt in
		h)
			echo -e "$helps\n"
			exit 0
		;;
		d)
			debug=true
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
		i)
			if [ -e $OPTARG ]; then
				inFile=$OPTARG
			else
				msg="!!! Invalid -i option, file not found.\n    File: $OPTARG"
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
if [ -z "$inFile" ]; then
	msg="!!! Missing mandatory -i option."
	echo -e "$helps\n$msg"
	exit
fi
if [ -z "$outFile" ]; then
	msg="!!! Missing mandatory -o option."
	echo -e "$helps\n$msg"
	exit
fi

# TEST =========================================================================

# RUN ==========================================================================

# Identify available chromosomes
IFS=' ' read -r -a chrlist <<< `cut -f 1 $inFile | sort | uniq`
echo -e " · Found ${#chrlist[@]} chromosomes."

# Default matrices
m_sums=""
m_tots=""
m_probs=""
m_cors=""
m_rocs=""
m_cons=""
m_bays=""
m_sigmas=""
m_ffs=""
m_cvs=""
m_exs=""

# Iterate over chromosomes
for chr in ${chrlist[@]}; do
	echo -e " >>> Working on $chr..."

	# Number of cutsite in the chromosome
	ncs=`cat $csList | awk -v chr=$chr '$1 == chr' | wc -l`

	# Sum columns awk program
	colSumPrg='@include "join"
	{
		for (i = 1; i <= NF; i++ )
			a[i] = a[i] + $i;
	}
	END {
		print join(a, 1, length(a));
	}'

	# Select and sum counts for current chromosome
	IFS=' ' read -r -a tots <<< \
		`cat $inFile | cut -f 5- | awk "$colSumPrg"`
	counts=`cat $inFile | awk -v chr=$chr '$1 == chr' | cut -f 5-`
	IFS=' ' read -r -a sums <<< \
		`echo -e "$counts" | awk "$colSumPrg"`

	# Calculate single-condition probability -----------------------------------
	probs=()
	for i in $(seq 0 `bc <<< "${#tots[@]}-1"`); do
		probs+=(`bc -l <<< "${sums[$i]}/${tots[$i]}/$ncs" | sed 's/^\./0./'`)
	done

	# Calculate cumulative of ratios (CoR) probability -------------------------
	cors=()
	for i in $(seq 0 `bc <<< "${#probs[@]}-1"`); do
		cors+=(`echo ${probs[@]:0:$i} | tr " " "+" | bc -l | sed 's/^\./0./'`)
	done

	# Calculate ratio of cumulatives (RoC) probability -------------------------
	rocs=()
	for i in $(seq 1 `bc <<< "${#probs[@]}"`); do
		rocn=`echo ${sums[@]:0:$i} | tr " " "+" | bc -l`
		roct=`echo ${tots[@]:0:$i} | tr " " "+" | bc -l`
		rocs+=(`bc -l <<< "$rocn/$roct/$ncs" | sed 's/^\./0./'`)
	done

	# Calculate conditional probability ----------------------------------------
	
	# Calculate conditional probability per cutsite
	conditional_p='{ OFS=FS=" "; OFMT="%.12f"

		split(Ts, T, ",");

		B=2;
		A=1;
		pconds=($B/T[B] + 1-$A/T[A] + ($B+T[A]-$A)/(T[B]+T[A]))/(1-$A/T[A]);

		for ( B=3; B <= NF; B++ ) {
			A=B-1;

			pcond=($B/T[B] + 1-$A/T[A] + ($B+T[A]-$A)/(T[B]+T[A]))/(1-$A/T[A]);
			pconds=pconds OFS pcond;
		}

		print pconds;
	}'
	cstots=`IFS=","; shift; echo "${tots[*]}";`
	pcons=`echo -e "$counts" | awk -v Ts="$cstots" "$conditional_p"`

	# Identify maximum per cutsite
	maxi_col='{ OFS=FS=" ";
		maxi=0;
		maxp=0;
		for ( i = 1; i <= NF; i++ ) {
			if ( $i > maxp ) {
				maxp=$i;
				maxi=i;
			}
		}
		print maxi;
	}'
	ccons=`echo -e "$pcons" | awk "$maxi_col"`

	# Average over region
	ccon=`echo -e "$ccons" | datamash mean 1`

	# Calculate Bayes probability ----------------------------------------------
	
	# Calculate Bayesian probability per cutsite
	bayes_p='{ OFS=FS=" "; OFMT="%.12f"

		sum=0;
		for ( i=1; i <= NF; i++ )
			sum=sum+$i;

		if ( 0 == sum ) next;

		p=$1/sum;
		for ( i=2; i <= NF; i++ )
			p=p OFS $i/sum;

		print p;
	}'
	pbays=`echo -e "$counts" | awk "$bayes_p"`
	
	# Identify maximum per cutsite
	cbays=`echo -e "$pbays" | awk "$maxi_col"`

	# Average over region
	cbay=`echo -e "$ccons" | datamash mean 1`

	# Calculate variability-based metrics --------------------------------------
	sigmas=()
	mus=()
	ffs=()
	cvs=()
	for i in $(seq 1 ${#tots[@]}); do
		sigma=`echo -e "$counts" | datamash sstdev $i`
		mu=`echo -e "$counts" | datamash mean $i`
		ffs+=(`bc -l <<< "$sigma^2 / $mu" | sed 's/^\./0./'`)
		cvs+=(`bc -l <<< "$sigma / $mu" | sed 's/^\./0./'`)
		sigmas+=($sigma)
		mus+=($mu)
	done

	# Calculate exclusivity ----------------------------------------------------
	
	# Calculate exclusivity per cutsite
	exclu_p='{ OFS=FS=" "; OFMT="%.12f"

		p=$1;

		for ( B=2; B <= NF; B++ ) {
			A=B-1;

			diff=$B-$A;

			if ( diff < 0 ) { diff=0; }

			p=p OFS diff;
		}

		print p;
	}'
	pexs=`echo -e "$probs" | awk "$exclu_p"`

	# Identify maximum per cutsite
	cexs=`echo -e "$pbays" | awk "$maxi_col"`

	# Average over region
	cex=`echo -e "$ccons" | datamash mean 1`

	# Make matrices ------------------------------------------------------------
	m_sums="$m_sums$chr `echo ${sums[@]}`\n"
	m_tots="$m_tots$chr `echo ${tots[@]}`\n"
	m_probs="$m_probs$chr `echo ${probs[@]}`\n"
	m_cors="$m_cors$chr `echo ${cors[@]}`\n"
	m_rocs="$m_rocs$chr `echo ${rocs[@]}`\n"
	m_cons="$m_cons$chr\t$ccon\n"
	m_bays="$m_bays$chr\t$cbay\n"
	m_sigmas="$m_sigmas$chr `echo ${sigmas[@]}`\n"
	m_mus="$m_mus$chr `echo ${mus[@]}`\n"
	m_ffs="$m_ffs$chr `echo ${ffs[@]}`\n"
	m_cvs="$m_cvs$chr `echo ${cvs[@]}`\n"
	m_exs="$m_exs$chr\t$cex\n"

	# Rank
	m_cons=`echo -e "$m_cons" | sort -k2`
	m_bays=`echo -e "$m_bays" | sort -k2`
	m_exs=`echo -e "$m_exs" | sort -k2`
done

if $debug; then
	echo -e "$m_sums" > $outFile".sums.matrix.tmp.txt"
	echo -e "$m_tots" > $outFile".tots.matrix.tmp.txt"
	echo -e "$m_probs" > $outFile".probs.matrix.tmp.txt"
	echo -e "$m_cors" > $outFile".cors.matrix.tmp.txt"
	echo -e "$m_rocs" > $outFile".rocs.matrix.tmp.txt"

	echo -e "$pcons" > $outFile".cons.$chr.matrix.tmp.txt"
	echo -e "$ccons" > $outFile".cons.$chr.normatrix.tmp.txt"
	echo -e "$m_cons" > $outFile".cons.ranked.tmp.txt"

	echo -e "$pbays" > $outFile".bays.$chr.matrix.tmp.txt"
	echo -e "$cbays" > $outFile".bays.$chr.normatrix.tmp.txt"
	echo -e "$m_bays" > $outFile".bays.ranked.tmp.txt"

	echo -e "$m_sigmas" > $outFile".sigmas.matrix.tmp.txt"
	echo -e "$m_mus" > $outFile".mus.matrix.tmp.txt"
	echo -e "$m_ffs" > $outFile".ffs.matrix.tmp.txt"
	echo -e "$m_cvs" > $outFile".cvs.matrix.tmp.txt"

	echo -e "$pexs" > $outFile".exs.matrix.$chr.tmp.txt"
	echo -e "$cexs" > $outFile".exs.normatrix.$chr.tmp.txt"
	echo -e "$m_exs" > $outFile".exs.ranked.tmp.txt"
fi

# Calculate centralities =======================================================

# awk programs
two_points='{ OFS=FS=" "; if ( 0 == NF ) next;
	print $1 "\t" $NF / $2;
}'
global='{ OFS=FS=" "; if ( 0 == NF ) next;
	sum=0;
	for ( B=3; B <= NF; B++) {
		A=B-1;
		sum=sum + $B / $A;
	}
	print $1 "\t" sum;
}'
gfixed='{ OFS=FS=" "; if ( 0 == NF ) next;
	sum=0;
	for ( B=3; B <= NF; B++)
		sum=sum + $B / $2;
	print $1 "\t" sum;
}'
log_two_points='{ OFS=FS=" "; if ( 0 == NF ) next;
	print $1 "\t" log($NF / $2);
}'
log_gfixed='{ OFS=FS=" "; if ( 0 == NF ) next;
	sum=0;
	for ( B=3; B <= NF; B++)
		sum=sum + log($B / $2);
	print $1 "\t" sum;
}'
diff_two_points='{ OFS=FS=" "; if ( 0 == NF ) next;
	print $1 "\t" $NF - $2;
}'
diff_gfixed='{ OFS=FS=" "; if ( 0 == NF ) next;
	sum=0;
	for ( B=3; B <= NF; B++)
		sum=sum + $B - $2;
	print $1 "\t" sum;
}'

# Probability-based centrality -------------------------------------------------

p_ranked_two_points=`echo -e "$m_probs" | awk "$two_points" | sort -k2`
p_ranked_global=`echo -e "$m_probs" | awk "$global" | sort -k2`
p_ranked_gfixed=`echo -e "$m_probs" | awk "$gfixed" | sort -k2`

# RoC centrality ---------------------------------------------------------------

roc_ranked_two_points=`echo -e "$m_rocs" | awk "$two_points" | sort -k2`
roc_ranked_global=`echo -e "$m_rocs" | awk "$global" | sort -k2`
roc_ranked_gfixed=`echo -e "$m_rocs" | awk "$gfixed" | sort -k2`

# CoR centrality ---------------------------------------------------------------

cor_ranked_two_points=`echo -e "$m_cors" | awk "$two_points" | sort -k2`
cor_ranked_global=`echo -e "$m_cors" | awk "$global" | sort -k2`
cor_ranked_gfixed=`echo -e "$m_cors" | awk "$gfixed" | sort -k2`

# Variance centrality ----------------------------------------------------------

var_ranked_two_points=`echo -e "$m_sigmas" | awk "$log_two_points" | sort -k2`
var_ranked_gfixed=`echo -e "$m_sigmas" | awk "$log_gfixed" | sort -k2`

# FF centrality ----------------------------------------------------------------

ff_ranked_two_points=`echo -e "$m_ffs" | awk "$diff_two_points" | sort -k2`
ff_ranked_gfixed=`echo -e "$m_ffs" | awk "$diff_gfixed" | sort -k2`

# CoV centrality ---------------------------------------------------------------

cv_ranked_two_points=`echo -e "$m_cvs" | awk "$diff_two_points" | sort -k2`
cv_ranked_gfixed=`echo -e "$m_cvs" | awk "$diff_gfixed" | sort -k2`

if $debug; then
	echo -e "$p_ranked_two_points" > $outFile".probs.ranked.2p.tmp.txt"
	echo -e "$p_ranked_global" > $outFile".probs.ranked.g.tmp.txt"
	echo -e "$p_ranked_gfixed" > $outFile".probs.ranked.gf.tmp.txt"

	echo -e "$cor_ranked_two_points" > $outFile".cors.ranked.2p.tmp.txt"
	echo -e "$cor_ranked_global" > $outFile".cors.ranked.g.tmp.txt"
	echo -e "$cor_ranked_gfixed" > $outFile".cors.ranked.gf.tmp.txt"

	echo -e "$roc_ranked_two_points" > $outFile".rocs.ranked.2p.tmp.txt"
	echo -e "$roc_ranked_global" > $outFile".rocs.ranked.g.tmp.txt"
	echo -e "$roc_ranked_gfixed" > $outFile".rocs.ranked.gf.tmp.txt"

	echo -e "$var_ranked_two_points" > $outFile".sigmas.ranked.2p.tmp.txt"
	echo -e "$var_ranked_gfixed" > $outFile".sigmas.ranked.gf.tmp.txt"

	echo -e "$ff_ranked_two_points" > $outFile".ffs.ranked.2p.tmp.txt"
	echo -e "$ff_ranked_gfixed" > $outFile".ffs.ranked.gf.tmp.txt"

	echo -e "$cv_ranked_two_points" > $outFile".cvs.ranked.2p.tmp.txt"
	echo -e "$cv_ranked_gfixed" > $outFile".cvs.ranked.gf.tmp.txt"
fi

# Rank centralities ============================================================

# Print header
out="P2p\tPg\tPgf\tCR2p\tCRg\tCRgf\tRC2p\tRCg\tRCgf\tcond\tbayes\tV2p\tVgf"
out=$out"\tF2p\tFgf\tCV2p\tCVgf\texc"
echo -e $out > $outFile".recap.txt"

# Merge and print rankings

paste \
	<(echo -e "$p_ranked_two_points" | cut -f 1) \
	<(echo -e "$p_ranked_global" | cut -f 1) \
	<(echo -e "$p_ranked_gfixed" | cut -f 1) \
	<(echo -e "$cor_ranked_two_points" | cut -f 1) \
	<(echo -e "$cor_ranked_global" | cut -f 1) \
	<(echo -e "$cor_ranked_gfixed" | cut -f 1) \
	<(echo -e "$roc_ranked_two_points" | cut -f 1) \
	<(echo -e "$roc_ranked_global" | cut -f 1) \
	<(echo -e "$roc_ranked_gfixed" | cut -f 1) \
	<(echo -e "$m_cons" | cut -f 1) \
	<(echo -e "$m_bays" | cut -f 1) \
	<(echo -e "$var_ranked_two_points" | cut -f 1) \
	<(echo -e "$var_ranked_gfixed" | cut -f 1) \
	<(echo -e "$ff_ranked_two_points" | cut -f 1) \
	<(echo -e "$ff_ranked_gfixed" | cut -f 1) \
	<(echo -e "$cv_ranked_two_points" | cut -f 1) \
	<(echo -e "$cv_ranked_gfixed" | cut -f 1) \
	<(echo -e "$m_exs" | cut -f 1) \
	>> $outFile".recap.txt"

# End --------------------------------------------------------------------------

echo -e "\n        ~ FIN ~        "
echo -e "   └[∵┌]└[ ∵ ]┘[┐∵]┘   "

################################################################################
