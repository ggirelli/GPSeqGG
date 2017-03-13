#!/usr/bin/Rscript

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 0.1.1
# Description: generates GPSeq final plots
# 
# ------------------------------------------------------------------------------



# DEPENDENCIES =================================================================

library(argparser)
source('umi_plot.functions.R')

# INPUT ========================================================================

# Create arguent parser
parser = arg_parser('Generate UMI plots.', name = 'umi_plot.R')

# Define mandatory arguments
parser = add_argument(parser, arg = 'dirpath',
	help = 'Experiment condition directory, contains UMI counts.')
parser = add_argument(parser, arg = 'experiment',
	help = 'Experiment ID (e.g., TK26).')
parser = add_argument(parser, arg = 'conditions',
	help = 'Comma-separated conditions (e.g., 400U2h,200U1h).')
parser = add_argument(parser, arg = 'cutsites',
	help = 'File containing the cutsite positions. (chr|pos)')
parser = add_argument(parser, arg = 'chrlengths',
	help = 'File containing the chromosome lengths. (chr|len)')
parser = add_argument(parser, arg = 'mask',
	help = 'File containing the masked regions. (num|chr|start|end|type)')

# Define elective arguments
parser = add_argument(parser, arg = '--neg', short = '-n',
	help = 'Negative condition label.',
	default = '', nargs = 1)
parser = add_argument(parser, arg = '--bin-size', short = '-i',
	help = 'Bin size in bp.',
	default = 1e6, type = class(0))
parser = add_argument(parser, arg = '--bin-step', short = '-t',
	help = 'Distance between the starting point of consecutive bins in bp.',
	default = 1e5, type = class(0))
parser = add_argument(parser, arg = '--num-proc', short = '-c',
	help = 'Number of cores for parallel computation.',
	default = 1, type = class(0))
parser = add_argument(parser, arg = '--suffix',
	help = 'Suffix to be added to output files.',
	default = '', nargs = 1)
parser = add_argument(parser, arg = '--rmChrX', flag = T,
	help = 'Remove ChrX from the analysis')
parser = add_argument(parser, arg = '--rmChrY', flag = T,
	help = 'Remove ChrY from the analysis')

# Parse arguments
p = parse_args(parser)

# Attach argument values to variables
attach(p['' != names(p)])

# Split comma-separated values
conditions <- unlist(strsplit(conditions, ',', fixed = T))

# CONSTANTS ====================================================================

# Color palette
col <- c('#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99',
	'#e31a1c', '#fdbf6f', '#ff7f00', '#cab2d6', '#6a3d9a')

# Bin definition string
binstring <- paste0('.bsize', bin_size, '.bstep', bin_step)

# RUN ==========================================================================

# File input -------------------------------------------------------------------

# Load UMItable
l <- load(paste0(dirpath, '/aux/',
	experiment, '.umi_table', binstring, '.RData'))

# Get list of present chromosomes
chr_list <- unique(unlist(lapply(umi_tab, names)))

# From chrN to N
chr_id_list <- unlist(lapply(chr_list,
	FUN = function(x) substr(x, 4, nchar(x))))
chr_id_list[chr_id_list == 'X'] <- 23
chr_id_list[chr_id_list == 'Y'] <- 24
chr_id_list <- as.numeric(chr_id_list)

# Remove chromosomes if requested
toRM = c()
if ( rmChrX ) toRM = which(23 == chr_id_list)
if ( rmChrY ) toRM = c(toRM, which(23 == chr_id_list))
if ( 0 != length(toRM) ) chr_id_list = chr_id_list[-toRM]

# Re-order chromosome list
chr_list <- chr_list[order(chr_id_list)]

# Negative condition id
if ( !exists('neg') ) neg = ''
neg_id = which(neg == names(l))
neg_id = ifelse(0 == length(neg_id), 0, neg_id)

# Read maskfile
n0n = numeric(0)
mf_ori = NULL
mf = NULL
if ( 0 != file.info(mask)$size ) {
	mf_ori = read.delim(mask, header = F, as.is = T)
	colnames(mf_ori) = c('num', 'chr', 'start', 'end', 'type')
	mf = mf_ori[, 2:5]
	mf$type = 'masked'
}

dt = mf
dt_field = 'type'
dt_levels = c('masked')
dt_title = 'Mask'

# PLOT =========================================================================

cat('Plotting ...\n')

# UMI thr study ----------------------------------------------------------------
cat(' >>> UMI distribution across condition ...\n')

plotUMIdistrib(dirpath, conditions, experiment, cutsites,
	ncores = num_proc)

# Sensed CSs -------------------------------------------------------------------
cat(' >>> Sensed cutsites profile ...\n')

pdf(paste0(dirpath, "plots/",
	experiment, '.sensed_cutsites', binstring, '.pdf'),
	width = 15, height = 15)
plotField(chr_list, conditions, 'ncs',
	'Genomic coordinates [nt]', 'Number of cut cutsites', default = 'cs',
	dt = dt, dt_field = dt_field, dt_levels = dt_levels, dt_title = dt_title)
graphics.off()

# UMI count profiles -----------------------------------------------------------
cat(' >>> UMI count profiles ...\n')

pdf(paste0(dirpath, "plots/", experiment, '.umi_absolute', binstring, '.pdf'),
	width = 15, height = 15)
plotField(chr_list, 
	conditions, 'n.sum',
	'Genomic coordinates [nt]', 'Number of unique UMIs',
	dt = dt, dt_field = dt_field, dt_levels = dt_levels, dt_title = dt_title)
graphics.off()

pdf(paste0(dirpath, "plots/", experiment, '.umi_mean', binstring, '.pdf'),
	width = 15, height = 15)
plotField(chr_list, conditions, 'n.mean',
	'Genomic coordinates [nt]', 'Number of unique UMIs per cutsite',
	dt = dt, dt_field = dt_field, dt_levels = dt_levels, dt_title = dt_title)
graphics.off()

pdf(paste0(dirpath, "plots/", experiment, '.umi_sd', binstring, '.pdf'),
	width = 15, height = 15)
plotField(chr_list, conditions, 'n.sd',
	'Genomic coordinates [nt]', 'SD of unique UMIs per cutsite',
	dt = dt, dt_field = dt_field, dt_levels = dt_levels, dt_title = dt_title)
graphics.off()

pdf(paste0(dirpath, "plots/", experiment, '.umi_CV', binstring, '.pdf'),
	width = 15, height = 15)
plotField(chr_list, conditions, 'n.sd/n.mean',
	'Genomic coordinates [nt]', 'CV of unique UMIs per cutsite',
	dt = dt, dt_field = dt_field, dt_levels = dt_levels, dt_title = dt_title)
graphics.off()

pdf(paste0(dirpath, "plots/", experiment, '.umi_FF', binstring, '.pdf'),
	width = 15, height = 15)
plotField(chr_list, conditions, 'n.sd^2/n.mean',
	'Genomic coordinates [nt]', 'FF of unique UMIs per cutsite',
	dt = dt, dt_field = dt_field, dt_levels = dt_levels, dt_title = dt_title)
graphics.off()

# UMI probability profiles -----------------------------------------------------
cat(' >>> UMI probability profiles ...\n')

pdf(paste0(dirpath, "plots/", experiment, '.p_umi_absolute', binstring, '.pdf'),
	width = 15, height = 15)
plotField(chr_list, conditions, 'p.sum',
	'Genomic coordinates [nt]', 'P(UMI)',
	dt = dt, dt_field = dt_field, dt_levels = dt_levels, dt_title = dt_title)
graphics.off()

pdf(paste0(dirpath, "plots/", experiment, '.p_umi_mean', binstring, '.pdf'),
	width = 15, height = 15)
plotField(chr_list, conditions, 'p.mean',
	'Genomic coordinates [nt]', '<P(UMI)>',
	dt = dt, dt_field = dt_field, dt_levels = dt_levels, dt_title = dt_title)
graphics.off()

# Same plots without the negative
if ( 0 != nchar(neg) & 0 != neg_id) {
	pdf(paste0(dirpath, "plots/",
		experiment, '.p_umi_absolute.no_neg', binstring,
		'.pdf'), width = 15, height = 15)
	plotField(chr_list, conditions[-neg_id], 'p.sum',
		'Genomic coordinates [nt]', 'P(UMI)',
		col = col[-1], dt = dt, dt_field = dt_field,
		dt_levels = dt_levels, dt_title = dt_title)
	graphics.off()

	pdf(paste0(dirpath, "plots/",
		experiment, '.p_umi_mean.no_neg', binstring, '.pdf'),
		width = 15, height = 15)
	plotField(chr_list, conditions[-1], 'p.mean',
		'Genomic coordinates [nt]', '<P(UMI)>',
		col = col[-1], dt = dt, dt_field = dt_field,
		dt_levels = dt_levels, dt_title = dt_title)
	graphics.off()
}

pdf(paste0(dirpath, "plots/", experiment, '.p_umi_sd', binstring, '.pdf'),
	width = 15, height = 15)
plotField(chr_list, conditions, 'p.sd',
	'Genomic coordinates [nt]', 'SD(P(UMI))',
	dt = dt, dt_field = dt_field, dt_levels = dt_levels, dt_title = dt_title)
graphics.off()

pdf(paste0(dirpath, "plots/", experiment, '.p_umi_CV', binstring, '.pdf'),
	width = 15, height = 15)
plotField(chr_list, conditions, 'p.sd/p.mean',
	'Genomic coordinates [nt]', 'CV(P(UMI))',
	dt = dt, dt_field = dt_field, dt_levels = dt_levels, dt_title = dt_title)
graphics.off()

pdf(paste0(dirpath, "plots/", experiment, '.p_umi_FF', binstring, '.pdf'),
	width = 15, height = 15)
plotField(chr_list, conditions, 'p.sd^2/p.mean',
	'Genomic coordinates [nt]', 'FF(P(UMI))',
	dt = dt, dt_field = dt_field, dt_levels = dt_levels, dt_title = dt_title)
graphics.off()

# Difference profiles ----------------------------------------------------------
cat(' >>> Difference profiles ...\n')

pdf(paste0(dirpath, "plots/", experiment, '.deltaCV', binstring, '.pdf'),
	width = 15, height = 15)
plotDiffField(chr_list, conditions, 'p.sd/p.mean',
	'Genomic coordinates [nt]', 'deltaCV(P(UMI))',
	col = col, neg_id = neg_id, rm_neg = T, y_hline = 0,
	dt = dt, dt_field = dt_field, dt_levels = dt_levels, dt_title = dt_title)
graphics.off()

pdf(paste0(dirpath, "plots/", experiment, '.deltaCV.extr', binstring, '.pdf'),
	width = 15, height = 15)
plotDiffField(chr_list, conditions, 'p.sd/p.mean',
	'Genomic coordinates [nt]', 'deltaCV(P(UMI))',
	col = col, neg_id = neg_id, rm_neg = T, only_extr = T, y_hline = 0,
	dt = dt, dt_field = dt_field, dt_levels = dt_levels, dt_title = dt_title)
graphics.off()

pdf(paste0(dirpath, "plots/", experiment, '.deltaFF', binstring, '.pdf'),
	width = 15, height = 15)
plotDiffField(chr_list, conditions, 'p.sd^2/p.mean',
	'Genomic coordinates [nt]', 'deltaFF(P(UMI))',
	col = col, neg_id = neg_id, rm_neg = T, y_hline = 0,
	dt = dt, dt_field = dt_field, dt_levels = dt_levels, dt_title = dt_title)
graphics.off()

pdf(paste0(dirpath, "plots/", experiment, '.deltaFF.extr', binstring, '.pdf'),
	width = 15, height = 15)
plotDiffField(chr_list, conditions, 'p.sd^2/p.mean',
	'Genomic coordinates [nt]', 'deltaFF(P(UMI))',
	col = col, neg_id = neg_id, rm_neg = T, only_extr = T, y_hline = 0,
	dt = dt, dt_field = dt_field, dt_levels = dt_levels, dt_title = dt_title)
graphics.off()

# Plot chr-wide FF -------------------------------------------------------------
cat(' >>> Chr-wide FF ...\n')

pdf(paste0(dirpath, "plots/", experiment, '.meanFF.chr_wide_ff.pdf'),
	width = 16, height = 6)
plotChrwideVar(dirpath, conditions, experiment, bin_size, cutsites,
	var_type = 'FF', col = col, neg_id = neg_id,
	rm.X = rmChrX, rm.Y = rmChrY, ncores = num_proc)
graphics.off()

# Plot chr-wide CV -------------------------------------------------------------
cat(' >>> Chr-wide CV ...\n')

pdf(paste0(dirpath, "plots/", experiment, '.meanCV.chr_wide_cv.pdf'),
	width = 16, height = 6)
plotChrwideVar(dirpath, conditions, experiment, bin_size, cutsites,
	var_type = 'CV', col = col, neg_id = neg_id,
	rm.X = rmChrX, rm.Y = rmChrY, ncores = num_proc)
graphics.off()

# Plot genomic overview plot ---------------------------------------------------
cat(' >>> Genomic overview ...\n')

pdf(paste0(dirpath, "plots/",
	experiment, '.genome_view.chromo.bsize', bin_size, '.pdf'),
	width = 15, height = 15)
plotGenome(dirpath, experiment, bin_size, bin_step, maskfile = mf_ori,
	rm.X = rmChrX, rm.Y = rmChrY, global = F
)
graphics.off()

pdf(paste0(dirpath, "plots/",
	experiment, '.genome_view.global.bsize', bin_size, '.pdf'),
	width = 15, height = 15)
plotGenome(dirpath, experiment, bin_size, bin_step, maskfile = mf_ori,
	rm.X = rmChrX, rm.Y = rmChrY, global = T
)
graphics.off()

# END --------------------------------------------------------------------------

################################################################################
