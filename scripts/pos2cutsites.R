#!/usr/bin/Rscript

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 1.0.0
# Description: convert genomic coordinates to cutsite position (in a window)
# 
# ------------------------------------------------------------------------------



# DEPENDENCIES =================================================================

suppressMessages(library(argparser))
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
suppressMessages(library(parallel))
suppressMessages(library(readr))

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

# Reset orphan positions
orph_fname = paste0(dirpath, 'orphans.txt')
if ( file.exists(orph_fname) ) {
	tmp = file.remove(orph_fname)
}

# Parallelize
btt <- mclapply(split(umi, umi_num),
	FUN = function(t, cs, bs) {
		chr <- t$chr[1]

		# Select cutsites on the current chromosome
		cs <- cs[cs$chr == chr,]

		# Calculate distance from closest cutsite and identify it
		ds <- rbindlist(mclapply(t$pos,
			FUN = function(x, cs) {
				if ( x > max(cs) ) {
					csid = which.max(cs)
					darr = x - cs[csid]
					return(data.frame(d = darr, id = csid))
				}

				csid = which(cs >= x)[1]
				darr = c(cs[csid] - x, x - cs[min(1, csid - 1)])
				if ( 2 == which.min(darr) ) csid = min(1, csid -1)
				return(data.frame(d = min(darr), id = csid))
			}, cs$pos
			, mc.cores = 1
		))

		# Count orphan locations
		n_locs = length(ds$d)
		n_orphs_loc = length(which(ds$d > bs/2))

		# Count orphan reads
		n_reads = sum(unlist(mclapply(t$seq,
			FUN = function(x) { length(unlist(strsplit(x, ' ', fixed = T))) }
			, mc.cores = 1)))
		n_orphs = sum(unlist(mclapply(t$seq[which(ds$d > bs/2)],
			FUN = function(x) { length(unlist(strsplit(x, ' ', fixed = T))) }
			, mc.cores = 1)))

		# Log orphan reads percentage
		cat(paste0(' >>> ', chr, ': ', n_orphs_loc, '/', n_locs,
			' (', round(n_orphs_loc / n_locs *100, 2),
			'%) orphan locations found.\n'))
		cat(paste0('     ', paste(rep(' ', nchar(chr)), collapse = ''), '  ',
			n_orphs, '/', n_reads,
			' (', round(n_orphs / n_reads *100, 2), '%) orphan reads found.\n'))

		# Save assignment and remove orphan reads
		t$mind <- ds$id
		torm <- which(ds$d > bs / 2)
		if ( 0 != length(torm) ) {
			# Write orphan read positions
			write.table(t[torm,], file = paste0(dirpath, 'orphans.txt'),
				append = T, col.names = F, row.names = F, quote = F, sep = "\t")
			t <- t[-torm,]
		}

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
			, mc.cores = 1
		))

		# Output table and orphan/total read counts
		return(list(c,
			data.frame(
				lorp = n_orphs_loc,
				ltot = n_locs,
				rorp = n_orphs,
				rtot = n_reads
			)))

	}, cs, bin_size
	, mc.cores = num_proc
)

# Prepare orphan reads log -----------------------------------------------------
or <- rbindlist(lapply(btt, FUN = function(x) { x[[2]] }))

# Calculate orphan location absolute count and ratio
n_orphan_loc = sum(or$lorp)
p_orphan_loc = round(n_orphan_loc / sum(or$ltot) * 100, 2)
cat(paste0(' · Found ', n_orphan_loc,
	' (', p_orphan_loc, '%) orphan locations in total.\n'))

# Calculate orphan read absolute count and ratio
n_orphan = sum(or$rorp)
p_orphan = round(n_orphan / sum(or$rtot) * 100, 2)
cat(paste0(' · Found ', n_orphan,
	' (', p_orphan, '%) orphan reads in total.\n'))

# Write to logfile
logfile = paste0(dirpath, condition, '.umi_prep_notes.txt')
log = c(
	paste0(n_orphan_loc, ' orphan locations (', p_orphan_loc, '%).'),
	paste0(n_orphan, ' orphan reads (', p_orphan, '%).'),
	paste0(bin_size / 2, ' nt as maximum distance from closest cutsite.')
)
write.table(log, logfile, row.names = F, col.names = F, quote = F, append = T)

# Retrieve and merge tables
bt <- rbindlist(lapply(btt, FUN = function(x) { x[[1]] }))

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
