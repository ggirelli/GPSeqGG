#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 1.0.1
# Description: 
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
 usage: ./reads_shuffle.sh [-h]

 Description:
  

 Mandatory arguments:
  

 Optional arguments:
  -h	Show this help page.
"

# Parse options
while getopts h opt "${bedfiles[@]}"; do
	case $opt in
		h)
			echo -e "$helps\n"
			exit 0
		;;
	esac
done

# Check mandatory options

# TEST =========================================================================

# RUN ==========================================================================

# End --------------------------------------------------------------------------

################################################################################
