#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 1.0.1
# Description: runs the global_centrality.sh script on shuffled bed files.
# 
# ------------------------------------------------------------------------------



# ENV VAR ======================================================================

export LC_ALL=C

# INPUT ========================================================================

# Cutsite file
csList="/media/gire/Data/BiCro-Resources/hg19.HindIII.txt"

# Directory with shuffled bed files
inDir="/home/gire/Desktop/BiCro/Code/src/"
inDir=$inDir"170424_global_centrality_chr_robustness/TK57_61_n100/"

# Bed files prefix
bedfiles=()
bedfiles+=($inDir"TK57_1min_GG__cutsiteLoc-umiCount")
bedfiles+=($inDir"TK58_5min_GG__cutsiteLoc-umiCount")
bedfiles+=($inDir"TK59_10min_GG__cutsiteLoc-umiCount")
bedfiles+=($inDir"TK60_15min_GG__cutsiteLoc-umiCount")
bedfiles+=($inDir"TK61_30min_GG__cutsiteLoc-umiCount")
# The shuffled bed files should be generated with beds_shuffle.py
# And have the following naming format: prefix.iterX.Yperc.bed

# Number of iterations
nIter=100

# Percentages of shuffled reads
IFS=',' read -r -a percs <<< "10,20,30,40"

# Number of cores
ncores=25

# Output directory
outDir="/home/gire/Desktop/BiCro/Code/src/"
outDir=$outDir"170424_global_centrality_chr_robustness/"

# Iterate on percentages
c=0
for perc in ${percs[@]}; do

	# Iterate on shufflings
	for i in $(seq 1 $nIter); do

		# Prepare script call
		script="../../dev/global_centrality.sh -c $csList"
		script=$script" -o "$outDir"TK57_61_gc_n"$nIter
		script=$script"/TK57_61.global_centrality.iter$i."$perc"perc"
		for bf in ${bedfiles[@]}; do
			script=$script" $bf.iter$i."$perc"perc.bed"
		done

		# Call script
		$script &

		# Increase job counter
		c=`bc <<< "$c+1"`
		if [ $c -ge $ncores ]; then
			wait
			c=0
		fi
	done
done

# END --------------------------------------------------------------------------

################################################################################
