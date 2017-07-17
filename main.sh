#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 1.1.0
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

# Save command line
echo "$0 $*" > "$outdir/CMD"
mkdir -p $outdir/log

# GLOBAL VARS ==================================================================

timestamp=`date +"%Y-%m-%e.%H-%M-%S"`
main_logpath=$outdir/log/$timestamp".main.log"

# START LOG ====================================================================
clear
{

# Start
echo -e "$settings
START\n====================="

source $moddir/main.library_run.mod.sh

# RUN BY LIBRARY ---------------------------------------------------------------

for exp in ${expv[@]}; do library_run $exp; done

# Merge library summary
exp0=`echo ${expv[0]} | cut -d ' ' -f 1`
echo -e "expID\t`head -n 1 "$outdir/$exp0/summary"`" > $outdir/summary
for exp in ${expv[@]}; do
	echo -e "$exp\t$(cat $outdir/$exp/summary | sed 1d)" >> $outdir/summary
done

# END --------------------------------------------------------------------------

echo -e "\n\n~~ fin ~~"

# END LOG ======================================================================

} &> >(tee $main_logpath)

################################################################################
