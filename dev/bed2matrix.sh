#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 1.0.1
# Description: 	merges bedfiles into a matrix.
# 				The score column is merged based on the position or the name.
# 
# ------------------------------------------------------------------------------



# ENV VAR ======================================================================

export LC_ALL=C

# DEPENDENCIES =================================================================

function join_by { local IFS="$1"; shift; echo "$*"; }

# INPUT ========================================================================

# Help string
helps="
 usage: ./bed2matrix.sh [-hn] [BEDFILEs]...

 Description:
  Merge bedfiles into a matrix. The score column is merged based on the positon
  given by the chr+start+end columns (default) or by the name column (-n option)

 Mandatory arguments:
  BEDFILEs	Bed file(s). Expected to be ordered per condition.

 Optional arguments:
  -h		Show this help page.
  -n		Merge bedfiles based on name instead of location.
"

# Default options
byName=false

# Parse options
while getopts hn opt "${bedfiles[@]}"; do
	case $opt in
		h)
			echo -e "$helps\n"
			exit 0
		;;
		n)
			byName=true
		;;
	esac
done

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

if [ 0 -eq ${#bedfiles[@]} ]; then
	msg="!!! No bedfiles provided."
	echo -e "$helps\n$msg"
	exit 1
fi

# TEST =========================================================================

# RUN ==========================================================================

if $byName; then
	# Merge by name ------------------------------------------------------------

	# Merge files
	mergePrg='
	BEGIN{
		OFS=FS="\t";
	}

	(FNR==NR){
		a[$4]=$0;
		nFiles=NF-4;
		next;
	}

	{
		if ($4 in a){
			a[$4] = a[$4] OFS $5;
		} else {
			a[$4] = $1 OFS $2 OFS $3 OFS $4;
			i = 0;
			while ( i < nFiles ) {
				a[$4] = a[$4] OFS 0;
				i++;
			}
			a[$4] = a[$4] OFS $5;
		}
	}

	END{
		for ( k in a ) {
			n=split(a[k], s, OFS);
			for ( i = 0; i < nFiles+5-n; i++ )
				a[k] = a[k] OFS 0;
			print a[k];
		}
	}'

	# Cycle through files and merge them
	merged=""
	for bf in ${bedfiles[@]}; do
		if [ -z "$merged" ]; then
			merged=`cat $bf | sed 1d`
		else
			bf=`cat $bf | sed 1d`
			merged=`awk "$mergePrg" <(echo -e "$merged") <(echo -e "$bf")`
		fi
	done

	# Output
	echo -e "$merged" | sort -k1.4,2
else
	# Merge by location --------------------------------------------------------

	# Add chr~start~end column
	addIDprg='{
		OFS=FS="\t";
		print $1"~"$2"~"$3 OFS $0
	}'

	# Merge files
	mergePrg='
	BEGIN{
		OFS=FS="\t";
	}

	(FNR==NR){
		a[$1]=$0;
		next;
	}

	{
		if ($1 in a){
			a[$1] = a[$1] OFS $6;
		} else {
			a[$1] = $1 OFS $2 OFS $3 OFS $4 OFS $5;
			i = 0;
			while ( i < nFiles ) {
				a[$1] = a[$1] OFS 0;
				i++;
			}
			a[$1] = a[$1] OFS $6;
		}
	}

	END{
		for ( k in a ) {
			n=split(a[k], s, OFS);
			for ( i = 0; i < nFiles+5+1-n; i++ )
				a[k] = a[k] OFS 0;
			print a[k];
		}
	}'

	# Cycle through files and merge them
	merged=""
	for bfi in $(seq 0 `bc <<< "${#bedfiles[@]}-1"`); do
		bf=${bedfiles[$bfi]}
		if [ -z "$merged" ]; then
			merged=`cat $bf | sed 1d | awk "$addIDprg"`
		else
			bf=`cat $bf | sed 1d | awk "$addIDprg"`
			merged=`awk -v nFiles=$bfi "$mergePrg" <(echo -e "$merged") <(echo -e "$bf")`
		fi
	done

	# Output
	echo -e "$merged" | cut -f 2- | sort -k1.4,2
fi

# End --------------------------------------------------------------------------

################################################################################
