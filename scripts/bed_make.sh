#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 1.0.0
# Description: generates BED files
# 
# ------------------------------------------------------------------------------



# ENV VAR ======================================================================

export LC_ALL=C

# INPUT ========================================================================

# Help string
helps="
 usage: ./bed_make.sh [-h] -o outDir [-c csSeq][-f csList][-t headerLine]

 Description:
  Generates bed file from UMI txt file with chr-start-UMIs-counts columns.

 Mandatory arguments:
  -o outDir	Condition folder, containing UMI counts.

 Optional arguments:
  -h		Show this help page.
  -c csSeq	Cutsite sequence. Mandatory with -f.
  -f csList	Cutsite list file. Mandatory with -c.
  -t headerLine	Text header line.
"

# Parse options
while getopts ho:c:f:t: opt; do
	case $opt in
		h)
			echo -e "$helps\n"
			exit 0
		;;
		o)
			out_dir=$OPTARG
			if [ ! -d "$OPTARG" ]; then
				msg="Output folder not found, creating it."
				mkdir -p $out_dir
			fi
		;;
		c)
			csSeq=$OPTARG
		;;
		f)
			if [ -e "$OPTARG" ]; then
				csList=$OPTARG
			else
				msg="Invalid -f option, file not found.\nFile: $OPTARG"
				echo -e "$helps\n$msg"
				exit 1
			fi
		;;
		t)
			headerLine=$OPTARG
		;;
	esac
done

# Check mandatory options
if [ -z "$out_dir" ]; then
	echo -e "$helps\n!!! Missing mandatory -o option.\n"
	exit 1
fi

# RUN ==========================================================================

# Count cutsite size
csLen=${#csSeq}

if [ -n "$csList" ]; then
	# Print bedfile header line
	echo -e "$headerLine" >  "$out_dir/UMIpos.unique.atcs.bed"

	# Obtain cutsites from list and compare with uniqued UMIs
	awkprogram='
	("chr23"==$1){ $1="chrX" }
	("chr24"==$1){ $1="chrY" }
	(FNR==NR){
		FS=OFS="\t";

		a[$1"~"$2]=$0"\t"NR;
		next
	}

	(23==$1){ $1="X" }
	(24==$1){ $1="Y" }
	{
		FS=OFS="\t";
		gsub(/ /, "", $1)
		gsub(/ /, "", $2)
		k="chr"$1"~"$2;
	}

	(k in a) {
		split(a[k], cs, "\t")
		n=split($3, umi, " ")
		print "chr"$1 OFS $2 OFS $2+cslen OFS "cs_" cs[3] OFS n
	}'
	awk -v cslen=$csLen "$awkprogram" <(cat "$csList") \
		<(cat "$out_dir/UMIpos.unique.atcs.txt") \
		>> "$out_dir/UMIpos.unique.atcs.bed"
else
	# Print bedfile header line
	echo -e "$headerLine" >  "$out_dir/UMIpos.unique.bed"

	# Without cutsite assignment
	awkprogram='
	{
		split($0, r, "\t");
		n=split(r[3], u, " ");
	}

	(23==r[1]){ r[1]="X" }
	(24==r[1]){ r[1]="Y" }

	(0!=n){
		FS=OFS="";
		r[1]="chr"r[1];
		print r[1],"\t",r[2],"\t",r[2]+cslen-1,"\tloc_",NR,"\t",n
	}
	'
	awk -v cslen=$cslength "$awkprogram" \
		<(cat "$out_dir/UMIpos.unique.txt" | tr -d " " | sort -nk 1) \
		> "$out_dir/UMIpos.unique.bed"
fi

# End --------------------------------------------------------------------------

################################################################################
