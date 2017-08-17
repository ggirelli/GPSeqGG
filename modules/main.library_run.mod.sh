#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 1.1.0
# Description: module for single library run. Run after options.mod.
# 
# ------------------------------------------------------------------------------



# START ========================================================================

function library_run() {

	# Identify experiment ID and log it.
	expID=$1

	# PREPARE DIRECTORY STRUCTURE ==============================================

	# Input folders
	in=$indir/indata_$expID && mkdir -p $in

	# Output folders
	out=$outdir/$expID && mkdir -p $out
	cout=$out/conditions && mkdir -p $cout
	pout=$out/plots && mkdir -p $pout
	xout=$out/aux && mkdir -p $xout
	log=$out/log && mkdir -p $log

	# Additional outputs
	outcontrol=$out/tmp && mkdir -p $outcontrol
	logpath="$log/$expID.log"

	# START PIPELINE ===========================================================
	{

	# START LOG ----------------------------------------------------------------
	line=""
	for (( i = 1; i < $(echo "# $expID #" | wc -c); i++ )); do
		line=$line"#"; done
	echo -e "\n"$line"\n# $expID #\n"$line"\n"

	# IDENTIFY CONDITIONS ------------------------------------------------------
	condv=$(cat $indir/patterns.tsv | grep "^$expID" | cut -f 2)

	# IDENTIFY DATA FILES ------------------------------------------------------

	find $indir -maxdepth 1 -type f -iname "*$expID*R[12]*" | sort > filelist
	numb_of_files=`cat filelist | wc -l`

	if [ 0 -eq $numb_of_files ]; then
		echo -e "!!! ERROR! No R1/R2 files found for $expID.\n"
		exit 1
	fi

	r1=`cat filelist | head -n1`
	echo "R1 is " $r1
	if [ $numb_of_files -eq 2 ]; then
	    r2=`cat filelist | tail -n1`
	    echo "R2 is " $r2
	fi
	rm filelist

	# QUALITY CONTROL ----------------------------------------------------------

	# Load quality control module
	source $moddir/main.qc.mod.sh
	execute_step $dontask 'quality control' quality_control

	# FILE GENERATION ----------------------------------------------------------

	# Load file generation module
	source $moddir/main.file_gen.mod.sh
	execute_step $dontask 'file generation' file_generation

	# PATTERN FILTERING --------------------------------------------------------

	# Load pattern filtering module
	source $moddir/main.pattern_filtering.mod.sh
	execute_step $dontask 'pattern_filtering' pattern_filtering

	# ALIGNMENT ----------------------------------------------------------------

	# Load alignment module
	source $moddir/main.alignment.mod.sh
	execute_step $dontask 'alignment' alignment

	# Filter SAM ---------------------------------------------------------------

	# Load sam filter module
	source $moddir/main.sam_filter.mod.sh
	execute_step $dontask 'SAM filtering' filter_sam

	# PREPARE UMI --------------------------------------------------------------

	# Load UMI preparation module
	source $moddir/main.umi_prep.mod.sh
	execute_step $dontask 'UMI preparation' prepare_umi

	# LIBRARY COMPLEXITY -------------------------------------------------------
	
	# Load library complexity module
	source $moddir/main.library_complexity.mod.sh
	execute_step $dontask 'Library complexity estimation' library_complexity

	# BIN ----------------------------------------------------------------------

	# Load binning module
	source $moddir/main.binning.mod.sh
	execute_step $dontask 'binning' bin_step

	# PLOT UMI -----------------------------------------------------------------

	# Load plot module
	source $moddir/main.plot.mod.sh
	time execute_step $dontask 'UMI plot' plot_umi

	# Clean --------------------------------------------------------------------

	if [ 2 -le $neatness ]; then
		echo -e "\n~ Cleaning..."
		rm -r $outcontrol $in
	fi

	} &> >(tee $log/$timestamp.$expID.log)

}

################################################################################
