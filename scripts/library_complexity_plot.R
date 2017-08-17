#!/usr/bin/Rscript

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 1.0.0
# Description: produces plots for library complexity estimation
# 
# ------------------------------------------------------------------------------



# DEPENDENCIES =================================================================

suppressMessages(library(argparser))
suppressMessages(library(ggplot2))

# INPUT ========================================================================

# Create arguent parser
parser = arg_parser('Generate plot for library complexity estimation.',
	name = 'library_complexity_plot.R')

# Define mandatory arguments
parser = add_argument(parser, arg = 'outDir',
	help = 'GPSeq dataset output folder.')
parser = add_argument(parser, arg = 'expID',
	help = 'Experiment ID.')

# Parse arguments
p = parse_args(parser)

# Attach argument values to variables
attach(p['' != names(p)])

# RUN ==========================================================================

# Aux directory
auxDir = paste0(outDir, '/aux/')

# Start plotting
pdf(paste0(outDir, '/plots/library_complexity.pdf'), width = 10, height = 10)

# Plot observed curve
data = read.delim(paste0(auxDir, 'lc.c_curve.txt'))
colnames(data) = toupper(colnames(data))
ggplot(data, aes(x = TOTAL_READS, y = DISTINCT_READS,
	color = CONDITION)) + geom_line()

# Plot estimation
data = read.delim(paste0(auxDir, 'lc.lc_extrap.txt'))
ggplot(data, aes(x = TOTAL_READS, y = EXPECTED_DISTINCT,
	color = CONDITION)) + geom_line()

# # Plot enrichment
# data = read.delim(paste0(auxDir, 'lc.bound_pop.txt'))
# ggplot(data, aes(x = CONDITION, y = median_estimated_unobs)
# 	) + geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci))

# Write output
graphics.off()

# END --------------------------------------------------------------------------

################################################################################
