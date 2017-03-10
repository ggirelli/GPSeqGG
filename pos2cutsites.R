#!/usr/bin/Rscript

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 0.1.1
# Description: convert genomic coordinates to cutsite position (in a window)
# 
# ------------------------------------------------------------------------------



# DEPENDENCIES =================================================================

library(argparser)
library(data.table)
library(ggplot2)
library(parallel)
library(readr)

# INPUT ========================================================================

# Create arguent parser
parser = arg_parser('Convert genomic coordinate to cutsite.',
	name = 'pos2cutsite.R')

# Define mandatory arguments
parser = add_argument(parser, arg = 'dirpath',
	help = 'Experiment ccondition directory, contains UMI counts.')
parser = add_argument(parser, arg = 'experiment',
	help = 'Experiment ID (e.g., TK26).')
parser = add_argument(parser, arg = 'condition',
	help = 'Condition folder name (e.g., 400U2h).')
parser = add_argument(parser, arg = 'cutsites',
	help = 'File containing the cutsite positions. (chr | pos)')

# Define elective arguments
parser = add_argument(parser, arg = '--bin-size', short = '-i',
	help = 'Range size in bp around the cutsite.',
	default = 40, type = class(0))
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

cat('\nGrouping locations per cutsite:\n')

# Read files -------------------------------------------------------------------

cat(' · Retrieving cutsites from the provided list ...\n')

# Read cutstite list
cs <- read.delim(cutsites, as.is = T, header = F)
colnames(cs) <- c('chr', 'pos')
cs$pos <- cs$pos

# Read grouped UMIs
cat(' · Retrieving grouped UMIs ...\n')
umi <- read.delim(paste0(dirpath, 'UMIpos.txt'), as.is = T, header = F)
colnames(umi) <- c('chr', 'pos', 'seq', 'qual')
umi = umi[order(as.numeric(umi$chr)),]
umi_num = umi$chr
umi$chr[umi$chr == 23] <- 'X'
umi$chr[umi$chr == 24] <- 'Y'
umi$chr <- paste0('chr', umi$chr)

cat(paste0(' · Range: +-', bin_size / 2, '\n'))

# Associate positions to cutsites ----------------------------------------------
cat(' · Associating locations to cutsites ...\n')
bt <- rbindlist(lapply(split(umi, umi_num),
	FUN = function(t, cs, bs) {
		chr <- t$chr[1]

		# Select cutsites on the current chromosome
		cs <- cs[cs$chr == chr,]

		# Calculate distance from closest cutsite
		ds <- unlist(mclapply(t$pos,
			FUN = function(x, cs) { min(abs(x - cs)) }, cs$pos
			, mc.cores = num_proc))

		# Plot
		# plot(ds, pch = 16, col = rgb(0, 0, 0, .3), cex = .5, log = 'y',
		# 	xlab = 'Distance from closest cutsite [nt]', ylab = 'Occurrences',
		# 	main = names(table(ds[ds != 0]))[which.max(table(ds[ds != 0]))])
		# abline(h = 70)
		# abline(h = 76, col=2)

		# Identify closest cutsite
		dsid <- unlist(mclapply(t$pos,
			FUN = function(x, cs) { which.min(abs(x - cs)) }, cs$pos
			, mc.cores = num_proc))

		# Log orphan reads percentage
		cat(paste0(' >>> ',
			chr, ': ', length(which(ds > bs/2)), '/', length(ds),
			' (', round(length(which(ds > bs/2)) / length(ds) *100, 2),
			'%) orphan reads found.\n'))

		# Save assignment and remove orphan reads
		t$mind <- dsid
		torm <- which(ds > bs / 2)
		if ( 0 != length(torm) ) t <- t[-torm,]

		# Move the UMIs together at the cutsite position
		c <- rbindlist(mclapply(split(t, t$mind),
			FUN = function(st, cs) {
				id <- st$mind[1]
				data.frame(
					chr = cs$chr[id],
					pos = cs$pos[id],
					seq = paste(st$seq, collapse = ' '),
					qual = paste(st$qual, collapse = ' '),
					stringsAsFactors = F
				)
			}, cs
			, mc.cores = num_proc
		))

		# Output
		return(c)

	}, cs, bin_size
))

# Re-format table
bt$pos <- as.numeric(bt$pos)
colnames(bt) <- c('chr', 'pos', 'seq', 'qual')
rownames(bt) <- NULL

# Update chromosome signatures
bt$chr[bt$chr == 'chrX'] <- 'chr23'
bt$chr[bt$chr == 'chrY'] <- 'chr24'
bt$chr <- as.numeric(lapply(bt$chr,
	FUN = function(x) { substr(x, 4, nchar(x)) }))

# Save output
cat(' · Writing output ...\n')
write.table(bt, file = paste0(dirpath, 'UMIpos.atcs.txt'),
	quote = F, row.names = F, col.names = F, sep = '\t')


# END --------------------------------------------------------------------------

################################################################################
