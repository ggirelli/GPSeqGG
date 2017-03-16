#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 0.1.1
# Description: analyze GPSeq sequencing data
# 
# Help page: ./main.sh -h
# 
# ------------------------------------------------------------------------------



# ENV VAR ======================================================================
export LC_ALL=C

# DEPENDENCIES =================================================================

# Code folder
scriptdir="`dirname ${BASH_SOURCE}`/scripts/"
moddir="`dirname ${BASH_SOURCE}`/modules/"

# Load functions
source $moddir/main.functions.sh

# CHECK OPTIONS ================================================================

# Load input options module
source $moddir/main.options.mod.sh

# Ask the user to double-check everything
check_settings

# PREPARE DIRECTORY STRUCTURE ==================================================

# Input folders
in=$indir/indata && mkdir -p $in

# Output folders
out=$outdir/ && mkdir -p $out
cout=$outdir/conditions && mkdir -p $cout
pout=$outdir/plot && mkdir -p $pout
xout=$outdir/aux && mkdir -p $xout

# Additional outputs
outcontrol=$outdir/tmp && mkdir -p $outcontrol
logpath="$out/$expID.log"

# Save command line
echo "$0 $*" > "$xout/CMD"

# START LOG ====================================================================
clear
{

# Start
echo -e "$settings
START\n=====================\n"

# IDENTIFY DATA FILES ==========================================================

find $indir -maxdepth 1 -type f -iname "*$expID*R[12]*" | sort > filelist
numb_of_files=`cat filelist | wc -l`

if [ 0 -eq $numb_of_files ]; then
	echo -e "!!! ERROR. No R1/R2 files found for $expID.\n"
	exit 1
fi

r1=`cat filelist | head -n1`
echo "R1 is " $r1
if [ $numb_of_files -eq 2 ]; then
    r2=`cat filelist | tail -n1`
    echo "R2 is " $r2
fi
rm filelist

# START PIPELINE ===============================================================

# QUALITY CONTROL --------------------------------------------------------------

# Load quality control module
source $moddir/main.qc.mod.sh
execute_step $dontask 'quality control' quality_control

# FILE GENERATION --------------------------------------------------------------

# Load file generation module
source $moddir/main.file_gen.mod.sh
execute_step $dontask 'file generation' file_generation

# PATTERN FILTERING ------------------------------------------------------------

# Load pattern filtering module
source $moddir/main.pattern_filtering.mod.sh
execute_step $dontask 'pattern_filtering' pattern_filtering

# ALIGNMENT --------------------------------------------------------------------

# Load alignment module
source $moddir/main.alignment.mod.sh
execute_step $dontask 'alignment' alignment

# Filter SAM -------------------------------------------------------------------

# Load sam filter module
source $moddir/main.sam_filter.mod.sh
execute_step $dontask 'SAM filtering' filter_sam


# PREPARE UMI ------------------------------------------------------------------

# Load UMI preparation module
source $moddir/main.umi_prep.mod.sh
execute_step $dontask 'UMI preparation' prepare_umi

# BIN --------------------------------------------------------------------------

# Load binning module
source $moddir/main.binning.mod.sh
execute_step $dontask 'binning' bin_step

# PLOT UMI ---------------------------------------------------------------------

# Load plot module
source $moddir/main.plot.mod.sh
time execute_step $dontask 'UMI plot' plot_umi

# END --------------------------------------------------------------------------

echo -e "\n\n~~ fin ~~"

# END LOG ======================================================================

} &> >(tee $logpath)

################################################################################
