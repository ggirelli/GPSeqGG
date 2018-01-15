#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 1.0.0
# Description: module for quality control recap with multiqc.
# 
# ------------------------------------------------------------------------------



# START ========================================================================

function qc_recap() {
	echo -e '\nQC recap\n====================='
	multiqc -o "$outdir" "$outdir"
}

################################################################################
