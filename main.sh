#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 1.0.1
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
main_logpath="$outdir/main.log"

# START LOG ====================================================================
clear
{

# Start
echo -e "$settings
START\n=====================\n"

source $moddir/main.library_run.mod.sh

# RUN BY LIBRARY ---------------------------------------------------------------

expIDs=$(cut -f 1 $indir/patterns.tsv | sort | uniq)
for expID in ${expIDs[@]}; do
	{
		library_run $expID
	} &> >(tee $outdir/$expID/$expID.log)
done

# END --------------------------------------------------------------------------

echo -e "\n\n~~ fin ~~"

# END LOG ======================================================================

} &> >(tee $main_logpath)

################################################################################
