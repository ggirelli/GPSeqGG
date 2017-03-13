#!/usr/bin/Rscript

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 0.1.1
# Description: bin UMIs
# 
# ------------------------------------------------------------------------------



# DEPENDENCIES =================================================================

suppressMessages(library(argparser))
suppressMessages(library(data.table))
suppressMessages(library(parallel))
suppressMessages(library(readr))

# INPUT ========================================================================

# Create arguent parser
parser = arg_parser('Bin UMIs.', name = 'umi_bin.R')

# Define mandatory arguments
parser = add_argument(parser, arg = 'dirpath',
	help = 'Experiment main directory, contains binned cutsites if available.')
parser = add_argument(parser, arg = 'experiment',
	help = 'Experiment ID (e.g., TK26).')
parser = add_argument(parser, arg = 'conditions',
	help = 'Comma-separated conditions (e.g., 400U2h,200U1h).')
parser = add_argument(parser, arg = 'chrlengths',
	help = 'File containing the chromosome lengths. (chr | len)')
parser = add_argument(parser, arg = '--cutsites', flag = T,
	help = 'Whether cutsites where used.')

# Define elective arguments
parser = add_argument(parser, arg = '--bin-size', short = '-i',
	help = 'Bin size in bp.',
	default = 1e6, type = class(0))
parser = add_argument(parser, arg = '--bin-step', short = '-t',
	help = 'Distance between the starting point of consecutive bins in bp.',
	default = 1e5, type = class(0))
parser = add_argument(parser, arg = '--num-proc', short = '-p',
	help = 'Number of cores for parallel computation.',
	default = 1, type = class(0))
parser = add_argument(parser, arg = '--suffix',
	help = 'Suffix to be added to output files.',
	default = '', nargs = 1)

# Parse arguments
p = parse_args(parser)

# Attach argument values to variables
attach(p['' != names(p)])

# Split comma-separated values
conditions = unlist(strsplit(conditions, ','))

# Set variables ----------------------------------------------------------------
suff = ''

# FUNCTIONS ====================================================================

bin_summary = function(i, bins, st) {
	# Statistic summary of a chromosome single bin UMIs.
	# 
	# Args:
	# 	i (int): bin index
	# 	bins (data.frame): bins
	# 	st (data.frame): UMI table
	# 
	
	# Identify current chromosome
	chr <- st$chr[1]

	# Identify current bin
	bin <- c(bins[i], bins[i]+bin_size)

	# Identify in-bin cut-cutsite
	ids <- which(st$pos >= bin[1] & st$pos < bin[2])

	# Select in-bin cut-cutsite
	bindata <- st[ids,]

	# Count in-bin cutsites
	ncs <- nrow(bindata)

	# Set default bindata variable
	if ( 0 == ncs )
		bindata <- data.frame(count = 0, p = 0)

	# Count expected in-bin cutsites
	cs <- 0
	if ( cutsites ) {
		cs <- csbin[[paste0('chr', chr)]][i]
	}
	if ( 0 == length(cs) ) cs <- 0

	# Absolute statistics summary for in-bin counts
	nsummary <- data.frame(
		t(c(summary(bindata$count))),
		sum = sum(bindata$count),
		sd = sd(bindata$count)
	)
	colnames(nsummary) <- c('min', 'q1', 'median',
		'mean', 'q3', 'max', 'sum', 'sd')
	if ( is.na(nsummary$sd) ) nsummary$sd <- 0

	# Relative statistics summary for in-bin counts
	psummary <- data.frame(
		t(c(summary(bindata$p))),
		sum = sum(bindata$p),
		sd = sd(bindata$p)
	)
	colnames(psummary) <- c('min', 'q1', 'median',
		'mean', 'q3', 'max', 'sum', 'sd')
	if ( is.na(psummary$sd) ) psummary$sd <- 0

	# Output
	return(data.frame(
		start = bin[1],
		mid = (bin[2] + bin[1]) / 2,
		end = bin[2],
		cs = cs,
		ncs = ncs,
		n = nsummary,
		p = psummary,
		stringsAsFactors = F
	))
}

chr_summary = function(st, c, bin_size, num_proc) {
	# Statistic summary of a chromosome binned UMIs.
	# 
	# Args:
	# 	st (data.frame): UMI table
	# 	c (data.frame): chromosome size table
	# 	bin_size (int): bin size in bp
	# 	num_proc (int): number of processors for parallel computation
	# 

	# Identify current chromosome
	chr <- st$chr[1]

	# Retrieve chromosome size
	size <- c$len[c$chr == chr]
	
	# Calculate bins
	if ( bin_step > bin_size ) bin_step <- bin_size
	bins <- seq(0, size - bin_size, by = bin_step)

	# Rename heterologous chromosomes
	if (chr == 23) chr <- 'X'
	if (chr == 24) chr <- 'Y'

	# Log
	cat(paste0(' >>> Working on chr', chr,
		' [', length(bins)-1, ' bins]\n'))

	# Per bin
	t <- rbindlist(lapply(1:(length(bins) - 1),
		FUN = bin_summary, bins, st
	))

	# Output
	return(data.frame(t, stringsAsFactors = F))
}

# RUN ==========================================================================

# Read files -------------------------------------------------------------------

# Retrieve chr length
cat(' >>> Retrieve chromosome lengths ...')
c <- read_delim(chrlengths,	'\t', col_names = c('chr', 'len'), col_types = 'ci')

# From chrN to N
c$chr[c$chr == 'chrX'] <- 'chr23'
c$chr[c$chr == 'chrY'] <- 'chr24'
c$chr <- as.numeric(unlist(lapply(c$chr,
	FUN = function(x) { substr(x, 4, nchar(x)) })))

# Retrieve binned cutsites, to normalize later on
if ( cutsites ) {
	rdata_path = paste0(dirpath, '/aux/cs_per_bin.bsize', bin_size,
		'.bstep', bin_step, '.RData')

	if ( file.exists(rdata_path) ) {
		load(rdata_path)
	} else {
		stop(paste0('ERROR: missing ', rdata_path, ' file.\n'))
	}
}

# Produce UMI probability table ------------------------------------------------

# Per condition
umi_tab <- lapply(conditions,
	FUN = function(condition) {

		# Retrieve UMIs
		cat(paste0('\n · Retrieve UMIs for "', condition , '"...\n'))
		fname <- 'UMIpos.unique'
		if ( cutsites ) fname <- paste0(fname, '.atcs')
		u <- read.delim(paste0(dirpath, '/conditions/', condition, '/',
			fname, suff, '.txt'), as.is = T, header = F)
		colnames(u) <- c('chr', 'pos', 'seq')

		# Add UMI count
		cat(paste0(' · Add UMI counts to "', condition , '"...\n'))
		u$count <- unlist(lapply(strsplit(u$seq, ' '), length))

		# Add UMI probability
		cat(paste0(' · Add CS probabilities to "', condition , '"...\n'))
		u$p <- u$count / sum(u$count)

		# Re-order chromosomes
		u <- u[order(u$chr),]

		# Group UMIs per bin ---------------------------------------------------
		cat(' · Group UMIs per bin ...\n')

		# Per chromosome
		r <- by(u, u$chr, FUN = chr_summary, c, bin_size, num_proc)
		names(r) <- paste0('chr', sort(unique(u$chr)))
		names(r)[names(r) == 'chr23'] <- 'chrX'
		names(r)[names(r) == 'chr24'] <- 'chrY'

		# Add chromosome column
		rt <- rbindlist(lapply(1:length(r),
			FUN = function(i) {
				rsub <- r[[i]]
				rsub$chr <- names(r)[i]
				return(as.data.frame(rsub, stringsAsFactors = F))
			}
		))

		# Output 1 -------------------------------------------------------------

		# Export multi-chromosome condition-specific UMI table
		fname = paste0(dirpath, '/conditions/', condition, '/', experiment, '.',
			condition, '.umi_table.bsize', bin_size, '.bstep', bin_step, '.tsv')
		write.table(rt, fname, quote = F, row.names = F, sep = '\t')

		return(r)
	}
)
names(umi_tab) <- conditions

# Output 2 ---------------------------------------------------------------------

# Save multi-condition umi_table
save('umi_tab', file = paste0(dirpath, '/aux/', experiment, '.umi_table',
	'.bsize', bin_size, '.bstep', bin_step, '.RData'))

# END --------------------------------------------------------------------------

################################################################################
