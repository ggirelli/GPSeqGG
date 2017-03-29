#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 1.0.0
# Description: module for file generation.
# 
# ------------------------------------------------------------------------------



# START ========================================================================

function file_generation() {
	echo -e 'File generation\n=====================\n'
	# Generate necessary files
	if [ -z "$r2" ]; then
		time $scriptdir/files_prepare.sh -t $threads -o $in -1 $r1
	else
		time $scriptdir/files_prepare.sh -t $threads -o $in -1 $r1 -2 $r2
	fi
}

################################################################################
