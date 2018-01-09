#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 1.0.0
# Description: filter alignment output (SAM file)
# 
# ------------------------------------------------------------------------------



# ENV VAR ======================================================================

export LC_ALL=C

# INPUT ========================================================================

# Help string
helps="
 usage: ./sam_filter.sh [-h] -i samfile

 Description:
  Filter alignment output. Allows to apply a number of filters, depending on the
  selected options. Use samtools view to perform filters. Use awk for line-by-
  -line file manipulations. Only chromosomes 1-24 (23 = X, 24 = Y) are kept.

 Mandatory arguments:
  -i samfile	Input SAM file to be filterd.

 Optional arguments:
  -h		Show this help page.
  -p 		Keep only primary alignments.
  -1 		Keep only R1, for PE-seq.
  -c 		Remove chimeric reads, for PE-seq.
  -m 		Keep only mapped reads.
  -q qual 	Keep only alignments with MAPQ >= qual.
  -C chrlist	Comma-separated list of chromosomes to remove. E.g., 1,3,X.
  -l cslength	Cutsite length for alignment shift.
  -t threads	Number of threads for parallelization.
"

# Default values
keep_mapped=false
keep_r1=false
keep_primary=false
rm_chimeric=false
min_mapq=0
threads=1
cslength=0

# Parse options
while getopts hp1cmi:q:C:l:t: opt; do
	case $opt in
		h)
			echo -e "$helps\n"
			exit 1
		;;
		p) # Keep primary alignments only
			keep_primary=true
		;;
		1) # Keep R1 only
			keep_r1=true
		;;
		c) # Remove chimeras
			rm_chimeric=true
		;;
		m) # Keep only mapped
			keep_mapped=true
		;;
		i) # Input SAM file
			if [ -f "$OPTARG" ]; then
				samfile="$OPTARG"
			else
				msg="!ERROR! Invalid value for option -f."
				msg="$msg\n        File not found: $OPTARG."
				echo -e "$helps\n $msg"
				exit 0
			fi
		;;
		q) # MAPQ filter
			if [ $OPTARG -ge 0 ]; then
				min_mapq=$OPTARG
			else
				msg="!ERROR! Invalid value for option -q."
				msg="$msg\n        qual must be a positive integer."
				echo -e "$helps\n $msg"
				exit 0
			fi
		;;
		C) # Comma-separated list of chromosomes
			chrlist=$OPTARG
		;;
		l)
			if [ 0 -gt "$OPTARG" ]; then
				echo -e "Enforcing minimum cutsite length of 0.\n"
			else
				cslength=$OPTARG
			fi
		;;
		t) # Threads
			if [ 0 -ge "$OPTARG" ]; then
				echo -e "Enforcing a minimum of 1 thread.\n"
			else
				threads=$OPTARG
			fi
		;;
		?) # Default
			echo -e "$helps\n !ERROR! Unrecognized option ${@:$OPTIND-1}."
			exit 0
		;;
	esac
done

# Check mandatory options
if [ -z "$samfile" ]; then
	echo -e "$helps\n !ERROR! Missing mandatory -i option."
	exit 0
fi

# Check that some software is installed
if [ ! -x "$(command -v samtools)" ]; then
	echo -e "gpseq-seq-gg requires samtools to work."
	exit 1
fi
if [ ! -x "$(command -v awk)" ]; then
	echo -e "gpseq-seq-gg requires awk to work."
	exit 1
fi
if [ ! -x "$(command -v datamash)" ]; then
	echo -e "gpseq-seq-gg requires datamash to work."
	exit 1
fi

# RUN ==========================================================================

# Identify input/output directory
outdir=$(dirname "$samfile")
fname=$(basename "$samfile" .sam)

tspath="$outdir/$fname.tmp.sam"
fspath="$outdir/$fname.filtered.sam"
hspath="$outdir/$fname.sam.header"
samtools view -H "$samfile" > "$hspath"

# Count input ------------------------------------------------------------------
acount=$(samtools view -c "$samfile")
echo " · Found $acount records."

# Primary alignments -----------------------------------------------------------
if ( $keep_primary ); then
	echo " · Keeping only primary alignments..."
	samtools view -F 256 -h "$samfile" -@ $threads \
		> "$tspath"
	wait
	mv "$tspath" $fspath
	
	tcount=$(samtools view -c "$fspath" -@ $threads)
	echo -e "   > Left with $tcount records."
	
	fcount=$(bc <<< "$acount - $tcount")
	>&2 echo "$fcount secondary alignments."
	acount=$tcount
fi

# Filter chimeras --------------------------------------------------------------
if ( $rm_chimeric ); then
	echo " · Removing chimeras..."
	awkprg='
	BEGIN { OFS=FS="\t"; }
	{
		# Setup read ID counter
		if ( $1 in c ) {
			c[$1]++;
		} else {
			c[$1] = 1;
		}

		# Store read alignments and print them when coupled and non-chimeric
		if ( $1 in a ) {
			split(a[$1], tt, OFS);
			if ( $3 == tt[3] ) {
				print a[$1] "\n" $0
			}
		} else {
			a[$1] = $0;
		}
	}
	END {
		# Print unpaired reads
		for ( k in c ) {
			if ( 1 == c[k] ) {
				print a[k];
			}
		}
	}'
	cp "$hspath" "$tspath"
	cat "$fspath" | grep -v "^\@" \
		| awk "$awkprg" >> "$tspath"
	wait
	mv "$tspath" "$fspath"

	tcount=$(samtools view -c "$fspath" -@ $threads)
	echo -e "   > Left with $tcount records."
	
	fcount=$(bc <<< "$acount - $tcount")
	>&2 echo "$fcount chimeric reads."
	acount=$tcount
fi

# Convert to BAM ---------------------------------------------------------------
echo " · Converting to BAM..."
bpath="$outdir/$fname.bam"
samtools view -b "$fspath" -@ $threads > "$bpath"
rm "$fspath"

# R1 ---------------------------------------------------------------------------
if ( $keep_mapped ); then
	echo " · Keeping only R1..."
	samtools view -f 64 -h "$bpath" -@ $threads \
		> "$outdir/$fname.tmp.bam"
	wait
	mv "$outdir/$fname.tmp.bam" "$bpath"

	tcount=$(samtools view -c "$bpath" -@ $threads)
	echo -e "   > Left with $tcount records."

	fcount=$(bc <<< "$acount - $tcount")
	>&2 echo "$fcount R2 reads."
	acount=$tcount
fi

# Mapped reads -----------------------------------------------------------------
if ( $keep_r1 ); then
	echo " · Keeping only mapped reads..."
	samtools view -F 4 -h "$bpath" -@ $threads \
		> "$outdir/$fname.tmp.bam"
	wait
	mv "$outdir/$fname.tmp.bam" "$bpath"

	tcount=$(samtools view -c "$bpath" -@ $threads)
	echo -e "   > Left with $tcount records."

	fcount=$(bc <<< "$acount - $tcount")
	>&2 echo "$fcount unmapped reads."
	acount=$tcount
fi

# MAPQ threshold ---------------------------------------------------------------
if ( $keep_r1 ); then
	echo " · Applying MAPQ threshold of $min_mapq..."
	samtools view -q $min_mapq -h "$bpath" -@ $threads \
		> "$outdir/$fname.tmp.bam"
	wait
	mv "$outdir/$fname.tmp.bam" "$bpath"

	tcount=$(samtools view -c "$bpath" -@ $threads)
	echo -e "   > Left with $tcount records."

	fcount=$(bc <<< "$acount - $tcount")
	>&2 echo "$fcount reads with MAPQ < $min_mapq."
	acount=$tcount
fi

# Convert back to sam ----------------------------------------------------------
echo " · Converting to SAM..."
samtools view "$bpath" -@ $threads > "$tspath"
rm "$bpath"

# Chromosome filter ------------------------------------------------------------
if [ -n "$chrlist" ]; then
	echo " · Discarding chromosomes [$chrlist]..."
	awkprg='
	BEGIN { OFS=FS="\t"; }
	{
		# Identify chromosomes to discard
		split(cl, ca, ",");

		# Print if not in chromosome list.
		if ( !($3 in ca) ) {
			# Output
			if ( $3+0 > 0 || $3 == "X" || $3 == "Y" ) { print $0; }
		}
	}'
	cat "$hspath" > "$fspath"
	cat "$tspath" | awk -v cl="$chrlist" "$awkprg" >> "$fspath"
	wait
	mv "$fspath" "$tspath"


	tcount=$(samtools view -c "$tspath" -@ $threads)
	echo -e "   > Left with $tcount records."

	fcount=$(bc <<< "$acount - $tcount")
	>&2 echo "$fcount reads from removed chromosomes."
	acount=$tcount
fi

# Count number of reads after filtering ----------------------------------------
nreads=$(samtools view -c "$tspath")
>&2 echo "$nreads reads left after filtering."

# Perform read shift -----------------------------------------------------------

# Print SAM file header
cat "$hspath" >> "$fspath"

# Identify reverse-strand reads with samtools -f 16 and apply CIGAR shift
echo " · Shifting reverse-strand alignments..."
awkprg='
BEGIN { OFS = FS = "\t"; }
{
	# Identify CIGAR bits
	split($6, cv, /[MIDNSHP=X]/, cs);

	# Setup default shift
	shift = 0;

	# Cycle through CIGAR bits and update shift
	for ( i = 1; i <= length(cv); i++ ) {
		if ( cs[i] == "M" || cs[i] == "D" ) {
			shift += cv[i];
		}
		if ( cs[i] == "I" ) {
			shift -= cv[i];
		}
	}

	# Shift position
	$4 = $4 + shift;
	print shift | "cat 1>&2"

	# Output
	print $0;
}'
samtools view -f 16 "$tspath" | awk "$awkprg" \
	>> "$fspath" 2> "$outdir/shifts.txt"
wait

# Average shift
aveshift=$(cat "$outdir/shifts.txt" | datamash mean 1)
rm "$outdir/shifts.txt"
msg="$aveshift bp of average position correction for reverse strand alignments."
>&2 echo "$msg"

# Identify non-reverse-strand reads with samtools -f 16 and apply CIGAR shift
if [[ 0 -ne $cslength ]]; then
	echo " · Shifting positive strand alignments..."
	awkprg='
	BEGIN { OFS = FS = "\t"; }
	{
		# Apply shift
		$4 = $4 - cslength;

		# Output
		print $0;
	}'
	samtools view -F 16 "$tspath" | \
		awk -v cslength=$cslength "$awkprg" >> "$fspath"
	wait
fi

# Remove temporary files
rm "$hspath" "$tspath"

# END --------------------------------------------------------------------------

################################################################################
