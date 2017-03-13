#!/usr/bin/Rscript

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 0.1.1
# Description: bin cutsites.
# 
# ------------------------------------------------------------------------------



# DEPENDENCIES =================================================================

suppressMessages(library(argparser))
suppressMessages(library(parallel))
suppressMessages(library(readr))

# INPUT ========================================================================

# Create arguent parser
parser = arg_parser('Bin cutsites.', name = 'cs_bin.R')

# Define mandatory arguments
parser = add_argument(parser, arg = 'dirpath',
	help = 'Directory to which the binned cutsites will be saved.')
parser = add_argument(parser, arg = 'cutsites',
	help = 'File containing the cutsite positions. (chr | pos)')
parser = add_argument(parser, arg = 'chrlengths',
	help = 'File containing the chromosome lengths. (chr | len)')

# Define elective arguments
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

# Parse arguments
p = parse_args(parser)

# Attach argument values to variables
attach(p['' != names(p)])

# RUN ==========================================================================

# Count cutsites per bin, to normalize later on
if ( file.exists(cutsites) ) {

	# Read files ---------------------------------------------------------------

	# Retrieve chromosome lengths
	c <- read.delim(chrlengths, as.is = T, header = F)
	colnames(c) <- c('chr', 'len')

	# Retrieve cutsite list
	el <- read.delim(cutsites, as.is = T, header = F)
	colnames(el) <- c('chr', 'pos')

	# Count cutsites per bin ---------------------------------------------------
	cat('Counting cutsites per bin, bin_size:', bin_size, ' ...\n')

	# Per chromosome
	csbin <- by(el, el$chr,
		FUN = function(st, c, bin_size) {
			# Working on chr
			chr <- st$chr[1]

			# Retrieve chromosome size
			size <- c$len[c$chr == chr]

			# Calculate size breaks
			bins <- seq(0, size - bin_size, by = bin_step)

			# Per bin
			t <- unlist(mclapply(1:(length(bins)-1),
				FUN = function(i, bins, st) {
					# Select bin
					bin <- c(bins[i], bins[i] + bin_size)

					# Identify cutsites in the bin
					ids <- which(st$pos >= bin[1] & st$pos < bin[2])

					# Count
					return(length(ids))
				}, bins, st
				, mc.cores = num_proc
			))

			# Output
			return(t)
		}, c, bin_size
	)

	# From chrN to N
	chr_id <- unlist(lapply(unique(el$chr),
		FUN = function(x) { substr(x, 4, nchar(x)) } ))
	chr_id[chr_id == 'X'] <- 23
	chr_id[chr_id == 'Y'] <- 24
	chr_id <- as.numeric(chr_id)

	# Reorder chromosomes
	csbin <- csbin[order(chr_id)]

	# Output -------------------------------------------------------------------

	# Save binning
	save('csbin', file = paste0(dirpath, '/cs_per_bin',
		'.bsize', bin_size, '.bstep', bin_step, '.RData'))
}

# END --------------------------------------------------------------------------

################################################################################
