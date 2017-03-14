#!/usr/bin/Rscript

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 0.1.1
# Description: filter and deduplicate UMIs
# 
# ------------------------------------------------------------------------------



# DEPENDENCIES =================================================================

suppressMessages(library(argparser))
suppressMessages(library(data.table))
suppressMessages(library(parallel))
suppressMessages(library(readr))

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
parser = add_argument(parser, arg = '--cutsites', short = '-cs',
	default = 0, nargs = 1,
	help = 'Binary flag for cutsite assignment. 1 for, and 0 for no cutsites')
parser = add_argument(parser, arg = '--platform', short = '-p',
	default = 'L', nargs = 1,
	help = 'Sequencing platform identifier.')
parser = add_argument(parser, arg = '--cutoff', short = '-co',
	default = 1, nargs = 1,
	help = 'Probability cutoff, compared to the automatic for filtering.')
parser = add_argument(parser, arg = '--emax', short = '-em',
	default = 1e-3, nargs = 1,
	help = 'Maximum error probability for filtering.')
parser = add_argument(parser, arg = '--eperc', short = '-ep',
	default = 20, nargs = 1,
	help = 'Maximum percentage of bases with emax error probability.')
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

# CONSTANTS ====================================================================

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

# Prepare filename
fname <- 'UMIpos'
if ( 1 == cutsites ) fname <- paste0(fname, '.atcs')

# Read UMI table
u <- read_delim(paste0(dirpath, fname, '.txt'),
	'\t', col_names=c('chr', 'pos', 'seq', 'qual'), col_types='iccc')
u$pos <- as.numeric(u$pos)

# Initialize -------------------------------------------------------------------

# Check UMI length
cat(' · Checking UMI length ...\n')
ulen = unique(unlist(mclapply(unlist(strsplit(u$seq, ' ', fixed = T)),
	FUN = nchar, mc.cores = num_proc)))
if ( 1 < length(ulen) ) {
	cat(paste0(' · Multiple UMI length detected: ',
		paste(ulen, collapse = ' '), ' [nt]\n'))
	cat(paste0(' · Using the average to calculate the threshold: ',
		mean(ulen), ' [nt]\n'))
} else {
	cat(paste0(' · UMI length is consistently ', ulen, ' nt.\n'))
}

# Count UMIs
cat(' · Counting UMIs ...\n')
n = unlist(mclapply(strsplit(u$seq, ' ', fixed = T),
	FUN = length, mc.cores = num_proc))

# Build quality data.frame
qabd = mk_qab_df(qab, qab_min)

# Filter based on self-match probability ---------------------------------------

if ( 0 < cutoff ) {
	cat(paste0(' · Filtering UMIs based on self-match probability.\n'))

	# Identify unique quals
	quals = unlist(mclapply(u$qual, FUN = strsplit, ' ', fixed = T
		, mc.cores = num_proc))
	unique_quals = unique(quals)

	# Calculate self-match probability for unique quals
	uq_ps = unlist(mclapply(unique_quals,
		FUN = function(x) {
			p_seq_match(x, x, x, x, qabd = qabd)
		}
		, mc.cores = num_proc))

	# Get self-match probability distribution
	qual_ps = uq_ps[match(quals, unique_quals)]

	# Calculate automatic cutoff
	auto_cutoff = quantile(qual_ps, .25) - 1.5 * IQR(qual_ps)

	# Compare with manual cutoff
	if ( 1 == which.min(c(auto_cutoff, cutoff)) ) {
		cat(paste0(' ·̣ Using automatic cutoff: ', auto_cutoff, '\n'))
	} else {
		cat(paste0(' ·̣ Automatic cutoff: ', auto_cutoff, '\n'))
		cat(paste0(' ·̣ Using manual cutoff: ', cutoff, '\n'))
	}
	cutoff = min(auto_cutoff, cutoff)

	# Count outliers
	perc_outliers = round(sum(qual_ps < cutoff) / sum(n) * 100, 2)
	cat(paste0(' >>> ', perc_outliers,
		'% of the UMIs are being discarded...\n'))

	# Identify outlier quals
	out_quals = unique_quals[which(uq_ps < cutoff)]

	# Discard outliers
	u = rbindlist(mclapply(1:nrow(u),
		FUN = function(i, u, out_quals) {
			# Select UMI table row
			row = u[i,]

			# Retrieve single sequences and quality strings
			ss = unlist(strsplit(u$seq[i], ' ', fixed = T))
			qs = unlist(strsplit(u$qual[i], ' ', fixed = T))

			# Identify non-outliers
			toRemove = which(qs %in% out_quals)

			if ( 0 != length(toRemove) ) {
				# Update row
				row$seq = paste(ss[-toRemove], collapse = ' ')
				row$qual = paste(qs[-toRemove], collapse = ' ')
			}

			# Output
			return(row)
		}, u, out_quals
		, mc.cores = num_proc
	))

	# Reset counters
	n = unlist(mclapply(strsplit(u$seq, ' ', fixed = T),
		FUN = length, mc.cores = num_proc))

	# Log
	cat(paste0(' >>> Left with ', sum(n), ' UMIs.\n'))
} else {
	cat(paste0(' · Skipped UMI filtering based on self-match probability.\n'))
}

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

# Retrieve qualities
quals = unlist(strsplit(u$qual, ' ', fixed = T))
uquals = unique(quals)
cat(paste0(' · Found ', length(uquals), ' unique quality strings (/',
	sum(n), ').\n'))

# Calculate quality of every uniqued quality string
pqs = lapply(uquals,
	FUN = function(q) {
		sum(Qstring_to_Perrs(q, qabd = qabd))
	}
)
rmq = uquals[which(pqs > ethr)]
nqkept = length(which(! quals %in% rmq))
cat(paste0(' >>> ', nqkept, '/', length(quals),
	' (', round(nqkept/length(quals)*100, 2), '%) UMIs pass the filter.', '\n'))

# Remove those that do not pass the threshold by checking from the overall index
cat(' · Removing UMIs ...\n')
u = as.data.frame(rbindlist(mclapply(1:nrow(u),
	FUN = function(i) {
		ss = unlist(strsplit(u$seq[i], ' ', fixed = T))
		qs = unlist(strsplit(u$qual[i], ' ', fixed = T))

		return(data.frame(
			chr = u$chr[i],
			pos = u$pos[i],
			seq = paste(ss[!qs %in% rmq], collapse = ' '),
			qual = paste(qs[!qs %in% rmq], collapse = ' '),
			stringsAsFactors = F
		))
	}
	, mc.cores = num_proc
)), stringsAsFactors = F)
cat(' >>> Removed.\n')

# Strict unique ----------------------------------------------------------------

# Perform strict unique
cat(' · Performing strict UMI deduplication...\n')
uniqued_seq = unlist(mclapply(u$seq,
	FUN = function(ss) {
		paste(unique(unlist(strsplit(ss, ' ', fixed = T))), collapse = ' ')
	}
	, mc.cores = num_proc
))

# Count delta-N
deltan = unlist(mclapply(strsplit(uniqued_seq, ' ', fixed = T),
	FUN = length, mc.cores = num_proc))

# Log
cat(paste0(' >>> ', round(sum(n - deltan) / sum(n) * 100, 2),
	'% UMIs identified as duplicates and removed.\n'))
cat(paste0(' >>> Remaining UMIs: ', sum(deltan), '\n'))

# Unique UMI list --------------------------------------------------------------

uu = data.frame(
	chr = u$chr,
	pos = u$pos,
	seq = uniqued_seq,
	stringsAsFactors = F
)

cat(' · Saving de-duplicated UMI list ...\n')
fname <- 'UMIpos.unique'
if ( 1 == cutsites ) fname <- paste0(fname, '.atcs')
write.table(as.matrix(uu), paste0(dirpath, fname, suffix, '.txt'),
	row.names = F, col.names = F, sep = '\t', quote = F)

# Plot -------------------------------------------------------------------------

cat(' · Generating UMI deduplication report ...\n')

# Output to umi_dedup.report.png
png(paste0(dirpath, 'umi_dedup.report.png'), width = 1200, height = 1200)
layout(matrix(c(1, 1, 2, 3), ncol = 2, byrow = T))

# Unique/total UMIs distribution
hist(deltan / n, prob = T, xlab = 'Unique/Total UMIs per cutsite',
	main = paste0(experiment, ' ~ ', condition), cex.lab = 1.5,
	cex.axis = 1.5, cex.main = 2)
lines(density(deltan / n, na.rm = T), col = 4, lty = 3)

# Pre-deduplication boxplot
boxplot(n, outline = F, ylim = c(1, max(n)), log = 'y',
	ylab = 'Number of UMIs per cutsite', cex.lab = 1.5,
	cex.axis = 1.5, cex.main = 2)
stripchart(n, vertical = T, add = T, method = 'j', pch = 20, cex = .5,
	col = rgb(0, 150, 135, 50, maxColorValue = 255))

# Post-deduplication boxplot
boxplot(deltan, outline = F, ylim = c(1, max(deltan)), log = 'y',
	ylab = 'Unique UMIs per cutsite', cex.lab = 1.5,
	cex.axis = 1.5, cex.main = 2)
stripchart(deltan, vertical = T, add = T, method = 'j', pch = 20, cex = .5,
	col = rgb(0, 150, 135, 50, maxColorValue = 255))

# Save to file
graphics.off()

# END --------------------------------------------------------------------------

################################################################################
