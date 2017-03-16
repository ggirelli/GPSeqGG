#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 0.1.0
# Description: calculates global centrality metrics
# 
# ------------------------------------------------------------------------------



# ENV VAR ======================================================================

export LC_ALL=C

# INPUT ========================================================================

# Help string
helps="
 usage: ./global_centrality.sh [-h] [BEDFILE]...

 Description:
  Calculate global centrality metrics.

 Mandatory arguments:
  BEDFILE	Bed file(s).

 Optional arguments:
  -h		Show this help page.
"

# Parse options
while getopts ho:c:f:t: opt; do
	case $opt in
		h)
			echo -e "$helps\n"
			exit 0
		;;
	esac
done

# Read list of bedfile paths


# RUN ==========================================================================


# End --------------------------------------------------------------------------

################################################################################
