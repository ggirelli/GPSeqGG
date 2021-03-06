#!/usr/bin/Rscript

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 2.0.0
# Description: filter and deduplicate UMIs
# 
# ------------------------------------------------------------------------------



# DEPENDENCIES =================================================================

suppressMessages(require(argparser))
suppressMessages(require(cowplot))
suppressMessages(require(data.table))
suppressMessages(require(ggplot2))
suppressMessages(require(parallel))
theme_set(theme_cowplot())

# INPUT ========================================================================

# Create arguent parser
parser = arg_parser('Deduplicate UMIs.', name = 'umi_dedupl.R')

# Define mandatory arguments
parser = add_argument(parser, arg = 'dirpath',
	help = 'Experiment condition directory, contains UMI counts.')
parser = add_argument(parser, arg = 'experiment',
	help = 'Experiment ID (e.g., TK26).')
parser = add_argument(parser, arg = 'condition',
	help = 'Condition folder name (e.g., 400U2h).')

# Define elective arguments
parser = add_argument(parser, arg = '--cutsites', short = '--cs',
	default = 0, nargs = 1,
	help = 'Binary flag for cutsite assignment. 1 for, and 0 for no cutsites')
parser = add_argument(parser, arg = '--platform', short = '-p',
	default = 'L', nargs = 1,
	help = 'Sequencing platform identifier.')
parser = add_argument(parser, arg = '--cutoff', short = '--co',
	default = 1, nargs = 1,
	help = 'Probability cutoff, compared to the automatic for filtering.')
parser = add_argument(parser, arg = '--emax', short = '--em',
	default = 1e-3, nargs = 1,
	help = 'Maximum error probability for filtering.')
parser = add_argument(parser, arg = '--eperc', short = '--ep',
	default = 20, nargs = 1,
	help = 'Maximum percentage of bases with emax error probability.')
parser = add_argument(parser, arg = '--num-proc', short = '-c',
	help = 'Number of cores for parallel computation.',
	default = 1, type = class(0))
parser = add_argument(parser, arg = '--num-reg', short = '-r',
	help = 'Number of regions per job during parallel computation.',
	default = 1000, type = class(0))
parser = add_argument(parser, arg = '--suffix',
	help = 'Suffix to be added to output files.',
	default = '', nargs = 1)

# Parse arguments
p = parse_args(parser)

# Attach argument values to variables
attach(p['' != names(p)])

setDTthreads(num_proc)

# CONSTANTS ====================================================================

# dirpath = "/mnt/data/preProcess_GPSeq/B48_1min/TK75/conditions/1min"
# experiment = "TK75"
# condition = "1min"
# platform = "L"
# cutoff = 0
# num_reg = 1000
# num_proc = 10
# cutsites = 1
# emax = 1e-3
# eperc = 20
# suffix = ""

qabs = list(
  S = list(
    qab = '!"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHI',
    min = 0,
    sep = '~'
  ),
  X = list(
    qab = ';<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefgh',
    min = -5,
    sep = '~'
  ),
  I = list(
    qab = '@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefgh',
    min = 0,
    sep = '~'
  ),
  J = list(
    qab = 'DEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefgh',
    min = 3,
    sep = '~'
  ),
  L = list(
    qab = '!"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJ',
    min = 0,
    sep = '~'
  )
)
qab <- qabs[[platform]]$qab
qab_min <- qabs[[platform]]$min
qab_sep <- qabs[[platform]]$sep

# FUNCTIONS ====================================================================

Qchar_to_Perr = function(c, qab, qab_min) {
	# Convert a quality char into probability of error.
	# 
	# Args:
	# 	c (char): quality character
	# 	qab (string): platform-specific quality alphabet
	# 	qab_min (int): Phred value of the first char in qab
	# 	
	
	10 ** (-(which(c == unlist(strsplit(qab, ''))) - 1 + qab_min)/10)
}

mk_qab_df = function(qab, qab_min) {
	# Build quality alphabet data frame.
	# 
	# Args:
	# 	qab (string): platform-specific quality alphabet string
	# 	qab_min (int): Phred value of the first char in the qab string
	# 
	# Returns:
	# 	The probability of error.
	# 
	
	# Make alphabet vector
	qabc = unlist(strsplit(qab, '', fixed = T))

	# Calculate single qab char error probability
	qs = unlist(lapply(qabc, FUN = Qchar_to_Perr, qab, qab_min))

	# Make data-frame
	qs = t(as.data.frame(qs, stringsAsFactors = F))
	colnames(qs) = qabc

	# Output
	return(qs)
}

Qstring_to_Perrs = function(s, qab = NULL, qab_min = NULL, qabd = NULL) {
	# Convert a quality string into a vector of error probabilities.
	# 
	# Notes:
	# 	Require either qab/qab_min or qabd.
	# 
	# Args:
	# 	s (string): quality string
	# 	qab (string): platform-specific quality alphabet
	# 	qab_min (int): Phred value of the first char in qab
	# 	qabd (data.frame): quality alphabet data.frame with Perrs
	# 	
	# Returns:
	# 	The vector of error probabilities.
	# 
	
	if ( !is.null(qabd) ) {
		# Use qab dictionary
		return(as.numeric(qabd[,unlist(strsplit(s, '', fixed = T))]))
	} else {
		# Calculate every Perr
		return(unlist(lapply(unlist(strsplit(s, '', fixed = T)),
			FUN = Qchar_to_Perr, qab, qab_min)))
	}
}

p_true_match = function(p1, p2) {
	# Calculate the probability that two matching bases are actually matching.
	# 
	# Args:
	# 	p1 (float): probability of base 1 being wrong
	# 	p2 (float): probability of base 2 being wrong
	# 
	# Returns:
	# 	The probability that a mismatch is actually a match, based on the two
	# 	single-based read quality values.
	# 
	
	return((1 - p1) * (1 - p2) + (p1 * p2 / 3))
}

p_false_mismatch = function(p1, p2) {
	# Calculate the probability that two different bases are actually matching.
	# 
	# Args:
	# 	p1 (float): probability of base 1 being wrong
	# 	p2 (float): probability of base 2 being wrong
	# 
	# Returns:
	# 	The probability that a mismatch is actually a match, based on the two
	# 	single-base read quality values.
	# 
	
	return(((1 - p1) * p2 / 3) + (p1 * (1 - p2) / 3) + (p1 * p2 / 2))
}

p_seq_match = function(s1, s2, q1, q2,
	qab = NULL, qab_min = NULL, qabd = NULL) {
	# Calculate the prob. that two sequences of equal length are identical.
	# 
	# Args:
	# 	s1 (string): first sequence
	# 	s2 (string): second sequence
	# 	q1 (string): first sequence quality string
	# 	q2 (string): second sequence quality string
	# 	qab (string): platform-specific quality alphabet
	# 	qab_min (int): Phred value of the first char in qab
	# 
	# Returns:
	# 	Probability that two sequences of equal length are identical, treating
	# 	separately matches and mismatches, based on the single-base read quality
	# 	values.
	# 
	
	# Check sequences length
	if ( nchar(s1) != nchar(s2) ) return(0)
	
	# Check quality string length
	if ( nchar(s1) != nchar(q1) || nchar(s2) != nchar(q2) ) {
		msg = 'The provided quality string(s) don\'t'
		msg = paste0(msg, ' match the corresponding sequence.\n')
		cat(msg)
		return(NULL)
	}

	if ( is.null(qabd) ) {
		# Build quality data-frame
		qabd = mk_qab_df(qab, qab_min)
	}
	
	# Calculate every single-base read error probability
	p1 <- Qstring_to_Perrs(q1, qabd = qabd)
	p2 <- Qstring_to_Perrs(q2, qabd = qabd)
	
	# Split sequences in bases
	ss1 <- unlist(strsplit(s1, ''))
	ss2 <- unlist(strsplit(s2, ''))
	
	# Identify matches and mismatches
	mid <- which(ss1 == ss2)
	eid <- which(ss1 != ss2)

	# Calculate single-base match probability
	ps <- c(
		p_true_match(p1[mid], p2[mid]),
		p_false_mismatch(p1[eid], p2[eid])
	)

	# Calculate string match probability
	p <- prod(ps)

	# Output
	return(p)
}

# RUN ==========================================================================

# Read UMI file ----------------------------------------------------------------

fname <- 'UMIpos'
if ( 1 == cutsites ) fname <- sprintf('%s.atcs', fname)

cat(' · Reading UMIs ...\n')
u <- fread(file.path(dirpath, sprintf('%s.txt', fname)),
	col.names=c('chr', 'pos', 'seq', 'qual'), nThread=num_proc)

cat(' >>> Pre-processing ...\n')
u$seq = mclapply(u$seq, function(x) unlist(strsplit(x, " ", fixed = T)),
	mc.cores = num_proc)
u$qual = mclapply(u$qual, function(x) unlist(strsplit(x, " ", fixed = T)),
	mc.cores = num_proc)
u$n = unlist(mclapply(u$seq, length, mc.cores = num_proc))

# Initialize -------------------------------------------------------------------

# Check UMI length
cat(' · Checking UMI length ...\n')

ulen = unique(nchar(unique(unlist(u$seq))))

if ( 1 < length(ulen) ) {
	cat(paste0('  >> Multiple UMI length detected: ',
		paste(ulen, collapse = ' '), ' [nt]\n'))
	cat(paste0('  >> Using the average to calculate the threshold: ',
		mean(ulen), ' [nt]\n'))
} else {
	cat(paste0('  >> UMI length is consistently ', ulen, ' nt.\n'))
}

# Count UMIs
cat(' · Counting UMIs ...\n')
log = paste0(u[, sum(n)], ' non-orphan reads.\n')

# Build quality data.frame
qabd = mk_qab_df(qab, qab_min)

# Filter based on single-base quality ------------------------------------------

cat(paste0(' · Filtering UMIs based on provided parameters.\n'))

# Check the parameters based on the selected quality alphabet
emin = Qchar_to_Perr(substr(qab, nchar(qab), nchar(qab)), qab, qab_min)
emax = max(emin, emax)
eperc = max(0, min(eperc, 100)) / 100
cat(paste0(' · Parameters used:\n'))
cat(paste0(' >>> emin: ', round(emin, 6), '\n'))
cat(paste0(' >>> emax: ', round(emax, 6), '\n'))
cat(paste0(' >>> eperc: ', eperc, '\n'))

# Calculate quality threshold
ethr = emin * ulen * eperc + emax * ulen * (1 - eperc)
cat(paste0(' · Calculated error probability threshold: ', round(ethr, 6), '\n'))

qData = data.table(q = unlist(u$qual)
	)[, .(n = .N), by = q]
qData$ps = unlist(mclapply(qData$q, function(q) {
	sum(Qstring_to_Perrs(q, qabd=qabd)) }))

qData[, toRM := F]
qData[ps > ethr, toRM := T]

nq = qData[, sum(n)]
nqkept = qData[(!toRM), sum(n)]
cat(sprintf(' >>> %d/%d (%.2f%%) UMIs pass the filter.\n',
	nqkept, nq, nqkept/nq*100))
log = paste0(log, sprintf('%d reads pass the read quality filter (%.2f%%).\n',
	nqkept, nqkept/nq*100))

# Remove those that do not pass the threshold by checking from the overall index
cat(' · Removing UMIs ...\n')
q2rm = qData[(toRM), q]

u[, toRM := F]
u[1 == n, toRM := qual %in% q2rm]
u = u[(!toRM)]
u[, toRM := NULL]

u$group = findInterval(seq(1, nrow(u)), seq(1, nrow(u), by = num_reg))
u = rbindlist(mclapply(split(u, u$group),  function(rowList, q2rm, maxGroup) {
	cat(sprintf("%d/%d (%.2f%%)\n", rowList[1, group], maxGroup,
		rowList[1, group]/maxGroup*100))
	rowList[, group := NULL]
	out = rbindlist(lapply(split(rowList, 1:nrow(rowList)), function(row) {
		if ( 1 == row$n ) return(row)
		toRM = unlist(row$qual) %in% q2rm
		if ( all(toRM) ) return(NULL)
		out = data.table(
			chr = row$chr,
			pos = row$pos,
			seq = list(unlist(row$seq)[!toRM]),
			qual = list(unlist(row$qual)[!toRM]),
			n = row$n - sum(toRM)
		)
		unlink(toRM)
		unlink(row)
		return(out)
	}))
	unlink(q2rm)
	unlink(rowList)
	return(out)
}, q2rm, u[, max(group)], mc.cores = num_proc, mc.preschedule = F))
cat(' >>> Removed.\n')

# Strict unique ----------------------------------------------------------------

# Perform strict unique
cat(' · Performing strict UMI deduplication...\n')
u$group = findInterval(seq(1, nrow(u)), seq(1, nrow(u), by = num_reg))
u = rbindlist(mclapply(split(u, u$group), function(rowList, maxGroup) {
	cat(sprintf("%d/%d (%.2f%%)\n", rowList[1, group], maxGroup,
		rowList[1, group]/maxGroup*100))
	rowList[, group := NULL]
	out = rbindlist(lapply(split(rowList, 1:nrow(rowList)), function(row) {
		uniqSeq = unique(unlist(row$seq))
		out = data.table(
			chr = row$chr,
			pos = row$pos,
			seq = list(uniqSeq),
			preN = row$n,
			postN = length(uniqSeq)
		)
		unlink(uniqSeq)
		unlink(row)
		return(out)
	}))
	unlink(rowList)
	return(out)
}, u[, max(group)], mc.cores = num_proc, mc.preschedule = F))

nu = u[, sum(preN)]
nurm = u[, sum(preN-postN)]
nuk = u[, sum(postN)]
cat(sprintf(' >>> %d (%.2f%%) UMIs identified as duplicates and removed.\n',
	nurm, nurm/nu*100))
log=sprintf('%s%d (%.2f%%) duplicated UMIs.\n%d UMIs left after deduplication.\n',
	log, nurm, nurm/nu*100, nuk)
cat(sprintf(' >>> Remaining UMIs: %d\n', nuk))

# Plot -------------------------------------------------------------------------

cat(' · Generating UMI deduplication report ...\n')

png(file.path(dirpath, 'umi_dedup.report.png'), width = 1200, height = 1200)
p1 = ggplot(u, aes(x = postN / preN)
	) + geom_histogram(breaks = seq(0, 1, by = 1/20),
		aes(y = ..density..), fill = NA, color = "black"
	) + geom_density(color = "blue", linetype = "dotted"
	) + xlab("Unique/Total UMIs per cutsite"
	) + xlim(0, 1
	) + ylab("Density"
	) + ggtitle(sprintf('%s ~ %s', experiment, condition))

p2 = ggplot(u, aes(y = preN, x = factor(1))
	) + geom_boxplot(outlier.shape = NA
	) + geom_jitter(width = .2, alpha = .5, size = .2, color = "#009687"
	) + xlab("") + ylab("UMIs per cutsite"
	) + scale_y_log10() + theme(axis.text.x = element_blank(),
		axis.ticks.x = element_blank())

p3 = ggplot(u, aes(y = postN, x = factor(1))
	) + geom_boxplot(outlier.shape = NA
	) + geom_jitter(width = .2, alpha = .5, size = .2, color = "#009687"
	) + xlab("") + ylab("Unique UMIs per cutsite"
	) + scale_y_log10() + theme(axis.text.x = element_blank(),
		axis.ticks.x = element_blank())

plot_grid(p1, plot_grid(p2, p3, ncol = 2), nrow = 2)
graphics.off()

# Export deduplicated UMIs -----------------------------------------------------

setnames(u, "postN", "counts")
u[, preN := NULL]
u$seq = unlist(lapply(u$seq, paste, collapse = " "))
cat(' · Saving de-duplicated UMI list ...\n')
fname <- 'UMIpos.unique'
if ( 1 == cutsites ) fname <- sprintf('%s.atcs', fname)
write.table(as.matrix(u),
	file.path(dirpath, sprintf('%s%s.txt', fname, suffix)),
	row.names = F, col.names = F, sep = '\t', quote = F)

# Write log --------------------------------------------------------------------

logfile = file.path(dirpath, sprintf('%s.umi_prep_notes.txt', condition))
log = unlist(strsplit(log, '\n', fixed = T))
write.table(log, logfile, row.names = F, col.names = F, quote = F, append = T)

# END --------------------------------------------------------------------------

################################################################################
