#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 2.0.1
# Description: analyze GPSeq sequencing data
# 
# Help page: ./main.sh -h
# 
# ------------------------------------------------------------------------------



# ENV VAR ======================================================================
export LC_ALL=C.UTF-8
export LANG=C.UTF-8
version="2.0.1"

# DEPENDENCIES =================================================================

# Code folder
libdir="`dirname ${BASH_SOURCE}`/lib/"
scriptdir="`dirname ${BASH_SOURCE}`/scripts/"
moddir="`dirname ${BASH_SOURCE}`/modules/"

# Load functions
source $moddir/main.functions.sh

# CHECK OPTIONS ================================================================

# Load input options module
source $moddir/main.options.mod.sh

# Ask the user to double-check everything
check_settings "$settings"

# Save command line
echo "$0 $*" > "$outdir/CMD"
mkdir -p $outdir/log

# GLOBAL VARS ==================================================================

timestamp=`date +"%Y-%m-%e.%H-%M-%S" | tr -d ' '`
main_logpath=$outdir/log/$timestamp".main.log"
echo -e "$main_logpath"

# START LOG ====================================================================
clear
{

# Start
echo -e "
#----------#
# SETTINGS #
#----------#
$settings

START\n====================="

source $moddir/main.library_run.mod.sh

# RUN BY LIBRARY ---------------------------------------------------------------

for exp in ${expv[@]}; do library_run $exp; done

# Merge library summary
exp0=`echo ${expv[0]} | cut -d ' ' -f 1`
echo -e "expID\t`head -n 1 "$outdir/$exp0/summary"`" > $outdir/summary
for exp in ${expv[@]}; do
	cat $outdir/$exp/summary | sed 1d | \
		awk -v expID="$exp" '{ print expID "\t" $0; }' >> $outdir/summary
done

# Recap qc reports with multiqc ------------------------------------------------

source $moddir/main.qc_recap.mod.sh
execute_step $dontask 'QC recap' qc_recap

# END --------------------------------------------------------------------------

echo -e "\n\n~~ fin ~~"

# END LOG ======================================================================

} &> >(tee $main_logpath)

################################################################################
