#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 1.0.0
# Description: module for quality control.
# 
# ------------------------------------------------------------------------------



# START ========================================================================

function quality_control() {
	echo -e 'Quality control\n=====================\n'
	# Produce quality control summarie(s)
	if [ -z "$r2" ]; then
		time $scriptdir/quality_control.sh -t $threads -o $xout -1 $r1
	else
		time $scriptdir/quality_control.sh -t $threads -o $xout -1 $r1 -2 $r2
	fi
}

################################################################################
