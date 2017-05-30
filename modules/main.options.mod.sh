#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 1.0.0
# Description: module for main script input option definition.
# 
# ------------------------------------------------------------------------------



# START ========================================================================

# Help string
helps="
usage: ./main.sh [-h][-w][-t threads] -i inDir -o outDir -e expID
 [-n][-a aligner][-g refGenome][-d bwaIndex][-x][-y][-q mapqThr]
 [-p platform][-u umiLength][-r csRange][-j emax][-k eperc][-z binSize]
 [-b binStep][-l csList][-m maskFile][-s chrLengths]

 Description:
  Run a step-by-step interactive GPSeq sequencing data analysis.

 Required files:
  Requires R1 (and R2 if paired-end sequencing) and a pattern files in the
  input directory (inDir). The patterns file should have a condition per row
  with condition name, pattern (scan_for_mateches format), cutsite sequence and
  non-genomic portion length, separated by tabulations.

 Mandatory arguments:
  -i indir	Input directory.
  -o outdir	Output directory. Created if not found.
  -e expID	Experiment ID.

 Optional arguments:
  -h	Show this help page.
  -w	Perform every step of the pipeline without asking.
  -n	Negative present. Expected label: neg;
  -x	Remove X chromosome after alignment.
  -y	Remove Y chromosome after alignment.
  -t threads	Number of threads for parallelization.
  -a aligner	Aligner. Either 'bwa' (default) or 'bowtie2'.
  -g refGenome	Path to reference genome file. Default: 'hg19'.
  -d bwaIndex	Path to BWA index file. Required if BWA is the selected aligner.
  -q mapqThr	Mapping quality threshold. Default: 30.
  -p platform	Sequencing platform. Default: 'L'.
  -u umilength	UMI sequence length. Default: 8.
  -r csRange	Range around cutsite for UMI assignment. Default: 40.
  -j emax	Maximum error probability for read quality filtering. Default: 1e-3.
  -k eperc	Maximum % of bases with emax error probability. Default: 20.
  -z binSize	Bin size. Default: 1e6.
  -b binStep	Bin step. Default: 1e5.
  -l csList	File with cutsite list. Columns: chr|pos. No header.
  -m maskFile	File with masked regions. Columns: id|chr|start|end. No header.
  -s chrLengths	File with chromosome lengths. chr|len. No header.
"

# Default values
dontask=0
threads=1
rmX=false
rmY=false
aligner='bwa'
refGenome='hg19'
mapqThr=30
platform='L'
umiLength=8
csRange=40
binSize=1e6
binStep=1e5
pthr=0
emax=1e-3
eperc=20

# Parse options
while getopts hwt:i:o:e:ng:a:d:xyq:p:u:r:z:b:j:k:l:m:s: opt; do
	case $opt in
		h)
			# Help
			echo -e "$helps"
			exit 0
		;;
		w)
			# Automatic
			dontask=1
		;;
		t)
			# Threads
			if [ 0 -ge "$OPTARG" ]; then
				echo -e "Enforcing a minimum of 1 thread.\n"
			else
				threads=$OPTARG
			fi
		;;
		i)
			# Input directory
			if [ -d "$OPTARG" ]; then
				indir=$OPTARG
			else
				msg="Invalid -i option, folder not found.\nFolder: $OPTARG"
				echo -e "$helps\n!!! ERROR! $msg"
				exit 1
			fi
		;;
		o)
			# Output directory
			outdir=$OPTARG
			if [ ! -d "$OPTARG" ]; then
				msg="Output folder not found, creating it."
				mkdir -p $outdir
			fi
		;;
		e)
			# Experiment ID
			expID=$OPTARG
		;;
		n)
			# Negative condition label
			neg='neg'
		;;
		g)
			# Reference genome
			refGenome=$OPTARG
		;;
		a)
			# Aligner
			if [ 'bwa' == "$OPTARG" -o 'bowtie2' == "$OPTARG" ]; then
				aligner=$OPTARG
			else
				msg="Invalid -a option. Available values: 'bwa', 'bowtie2'."
				echo -e "$helps\n!!! ERROR! $msg"
				exit 1
			fi
		;;
		d)
			# File with cutsite list
			if [ -e $OPTARG ]; then
				bwaIndex=$OPTARG
			else
				msg="Invalid -d option, file not found.\nFile: $OPTARG"
				echo -e "$helps\n!!! ERROR! $msg"
				exi 1
			fi
		;;
		x)
			# Remove chrX after alignment
			rmX=true
		;;
		y)
			# Remove chrY after alignment
			rmY=true
		;;
		q)
			# Mapping quality threshold
			mapqThr=$OPTARG
		;;
		p)
			# Sequencing platform
			platform=$OPTARG
		;;
		u)
			# UMI length in nt
			umiLength=$OPTARG
		;;
		r)
			# Range around cutsite
			csRange=$OPTARG
		;;
		z)
			# Bin size in nt
			binSize=$OPTARG
		;;
		b)
			# Bin step nin nt
			binStep=$OPTARG
		;;
		j)
			emax=$OPTARG
		;;
		k)
			eperc=$OPTARG
		;;
		l)
			# File with cutsite list
			if [ -e $OPTARG ]; then
				csList=$OPTARG
			else
				msg="Invalid -l option, file not found.\nFile: $OPTARG"
				echo -e "$helps\n!!! ERROR! $msg"
				exit 1
			fi
		;;
		m)
			# File with masked regions
			if [ -e $OPTARG ]; then
				maskFile=$OPTARG
			else
				msg="Invalid -m option, file not found.\nFile: $OPTARG"
				echo -e "$helps\n!!! ERROR! $msg"
				exit 1
			fi
		;;
		s)
			# File with chromosome lengths
			if [ -e $OPTARG ]; then
				chrLengths=$OPTARG
			else
				msg="Invalid -s option, file not found.\nFile: $OPTARG"
				echo -e "$helps\n!!! ERROR! $msg"
				exit 1
			fi
		;;
	esac
done

# Check mandatory options
if [ -z "$indir" ]; then
	echo -e "$helps\n!!! ERROR! Missing mandatory -i option.\n"
	exit 1
fi
if [ -z "$outdir" ]; then
	echo -e "$helps\n!!! ERROR! Missing mandatory -o option.\n"
	exit 1
fi
if [ -z "$expID" ]; then
	echo -e "$helps\n!!! ERROR! Missing mandatory -e option.\n"
	exit 1
fi

# Additional checks
if [ "bwa" == "$aligner" -a -z "$bwaIndex" ]; then
	echo -e "$helps\n!!! ERROR! Missing mandatory -d option.\n"
	exit 1
fi

# Check pat_files --------------------------------------------------------------

if [ -e "$indir/pat_files" ]; then

	# Count columns
	for i in $(seq 1 `wc -l $indir/pat_files | cut -d " " -f 1`); do
		# Retrieve pattern row
		row=`cat $indir/pat_files | head -n $i | tail -n 1`

		# Count fields
		nc=`echo "$row" | awk '{ n=split($0, t, "\t"); print n; }' | uniq`

		# Check field number
		if [ 4 -ne $nc ]; then
			msg="!!! Missing columns in pat_files, row $i.\n"
			msg="$msg    Expected 4 columns, found $nc."
			echo -e "$helps\n$msg"
			exit 1
		fi
	done

	conds=`cut -f 1 $indir/pat_files | tr '\n' ',' | sed -r 's/^(.*),$/\1/'`
	IFS=',' read -r -a condv <<< "$conds"
else
	msg="!!! ERROR! Missing pat_files.\n"
	msg=$msg"    Standard pat_files row:\n"
	msg=$msg"    conditionLabel	linkerPattern	cutsiteSeq	linkerSize\n"
	echo -e "$helps\n$msg"
	exit 1
fi

# Check settings ---------------------------------------------------------------

# Settings header
settings="
# SETTINGS
"
if [ 1 -eq $dontask ]; then
	settings="$settings
 Perform EVERY step of the pipeline.
"
fi

# Settings mandatory options
settings="$settings
 Input directory:\t$indir
 Output directory:\t$outdir
 
 Experiment ID: $expID
 Conditions: $conds"
if [ -n "$neg" ]; then
	settings="$settings\n Negative: $neg"
else
	settings="$settings\n!!! No negative condition."
fi

settings=$settings"\n\n Pattern instructions:\n"
patterns=" conditionLabel\tlinkerPattern\tcutsiteSeq\tlinkerSize\n"
patterns=$patterns" "$(cat -t $indir/pat_files | sed "s/\^I/\\\t/g")
patterns="$(echo -e ' '$patterns | sed ':a;N;$!ba;s/\n/\n /g')"
patterns="$(column -c 4 -s '~' -t <(echo "$patterns" | sed 's/\t/~/g'))"
settings=$settings"${patterns}"

# Other settings
settings="$settings

 Threads: $threads

 Reference genome: $refGenome
 Aligner: $aligner"
if [ -n $bwaIndex ]; then
	settings="$settings\n BWA index: $bwaIndex"
fi

settings="$settings
 MAPQ threshold: $mapqThr
 Platforhm: $platform

 UMI length: $umiLength nt
 Cutsite range: $csRange nt
 Max error probability: $emax
 Max emax bases percentage: $eperc%

 Bin size: $binSize nt
 Bin step: $binStep nt

"
# Chromosome removal
if $rmX; then
	settings="$settings!!! Removing chrX after alignment.\n"
fi
if $rmY; then
	settings="$settings!!! Removing chrY after alignment.\n"
fi

# Additional files
if [ -n "$csList" ]; then
	settings="$settings\n Cutsites:\n $csList\n"
fi
if [ -n "$maskFile" ]; then
	settings="$settings\n Masked regions:\n $maskFile\n"
fi
if [ -n "$chrLengths" ]; then
	settings="$settings\n Chromosome lengths:\n $chrLengths\n"
fi

# Print settings
clear
echo -e "$settings"

################################################################################
