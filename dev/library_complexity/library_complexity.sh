#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 1.0.0
# Date: 20170816
# Project: GPSeq-library-complexity
# 
# Description: estimate GPSeq library complexity and species richness using the
#              preseq tool [1].
# 
# References:
#  [1]: Daley & Smith, Apr 2013, Nat Meth 10(4):325 (DOI: 10.1038/NMETH.2375)
# 
# ------------------------------------------------------------------------------



# ENV VAR ======================================================================

export LC_ALL=C

# PARAMS =======================================================================


# Help string
helps="
usage: ./library_complexity.sh [-h] -i input [-d default][-f]

 Description:
  Estimate GPSeq library complexity and species richness using the preseq tool.

 Mandatory arguments:
  -i indir	Condition folder, containing the UMIpos.unique file.

 Optional arguments:
  -h	Show this help page.
  -f	Flag option.
  -d default	Param with default value. Default: 10
"

# Default values
default=10
flag=false

# Parse options
while getopts hfi:d: opt; do
	case $opt in
		h)
			# Help page
			echo -e "$helps"
			exit 0
		;;
		f)
			# Flag
			flag=true
		;;
		i)
			# Input
			input=$OPTARG
		;;
		d)
			# Default
			default=$OPTARG
		;;
	esac
done

# Check mandatory options
if [ -z "$input" ]; then
	echo -e "$helps\n!!! ERROR! Missing mandatory -i option.\n"
	exit 1
fi

# Additional checks
# ...

# RUN ==========================================================================

# END ==========================================================================

################################################################################
