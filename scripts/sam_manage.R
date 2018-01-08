#!/usr/bin/Rscript

# ------------------------------------------------------------------------------
#
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 1.1.0
# Description: script to manipulate alignment output (SAM file)
# 
# ------------------------------------------------------------------------------



# DEPENDENCIES =================================================================

suppressMessages(library(argparser))
suppressMessages(library(data.table))
suppressMessages(library(parallel))
suppressMessages(library(readr))

# PARAMS =======================================================================

# Create arguent parser
parser = arg_parser('Filter alignment output.', name = 'sam_filter.R')

# Define mandatory arguments
parser = add_argument(parser, arg = 'dirpath',
	help = 'Experiment condition directory, contains UMI counts.')
parser = add_argument(parser, arg = 'experiment',
	help = 'Experiment ID (e.g., TK26).')
parser = add_argument(parser, arg = 'condition',
	help = 'Condition folder name (e.g., 400U2h).')

# Define elective arguments
parser = add_argument(parser, arg = '--no-cutsite', flag = T,
	help = 'Used blunt linker, thus not cutsite sequence.')
parser = add_argument(parser, arg = '--cutsite', short = '-cs',
	default = 'AAGCTT',
	help = 'Cutsite sequence (e.g., for HindIII AAGCTT).')
parser = add_argument(parser, arg = '--mapq-thr', short = '-mt',
	help = 'Mapping quality threshold.',
	default = 30, type = class(0))
parser = add_argument(parser, arg = '--num-proc', short = '-c',
	help = 'Number of cores for parallel computation.',
	default = 1, type = class(0))
parser = add_argument(parser, arg = '--suffix',
	help = 'Suffix to be added to output files.',
	default = '', nargs = 1)
parser = add_argument(parser, arg = '--rmChr', short = '-r',
	default = '', nargs = 1,
	help = 'Comma separated chromosomes to remove.')

# Parse arguments
p = parse_args(parser)

# Attach argument values to variables
attach(p['' != names(p)])

# Additional manipulations
rmChr = unlist(strsplit(rmChr, ',', fixed = T))

# Output notes
notes = ""

# RUN ==========================================================================

setwd(dirpath)

# Read input -------------------------------------------------------------------

cat(paste0('\nFiltering SAM file:\n'))
cat(paste0(' · Collecting aligned reads.\n'))
t <- read_delim(paste0(condition, '.linkers.tmp.sam'), '\t',
	col_names = c('qname', 'flag', 'chr', 'pos', 'mapq', 'cigar', 'rnext',
		'pnext', 'tlen', 'lseq', 'seq'),
	col_types = 'ciciicciiccc')

# Reset position of reverse complement strand ----------------------------------

# Update binary flags
uf2 <- sort(unique(mt$flag))
ubf2 <- ubf[which(uf %in% uf2)]

# Identify rcs alignments
cat(' · Retrieving 3\' position ...\n')

# Identify RC reads from the flag (shift based on CIGAR)
rc <- unlist(mclapply(ubf2,
	FUN = function(x) { as.numeric(substr(x, 8, 8)) }
	, mc.cores = num_proc))
rc <- uf2[rc == 1]

# Calculate shift per CIGAR
move <- rbindlist(mclapply(unique(mt$cigar[mt$flag %in% rc]),
	FUN = function(x) {

		# Retrieve CIGAR letters
		clet <- unlist(strsplit(
			paste(unlist(strsplit(x, '[0-9]')), collapse=''),
			'', fixed = T))

		# Retrieve CIGAR numbers
		cnum <- as.numeric(unlist(strsplit(x, '[MIDNSHP=X]')))

		# Couple letters and numbers
		tab <- unlist(lapply(unique(clet),
			FUN = function(char) sum(cnum[clet == char])))
		names(tab) <- unique(clet)

		# Output CIGAR shift
		out = data.frame(
			cigar = x,
			move = sum(c(tab['M'], tab['D'], - tab['I']), na.rm = T)
		)

		return(out)
	}
	, mc.cores = num_proc
))
cat(' · Average position correction for reverse strand alignments:',
	round(mean(move$move, na.rm = T), 2), 'bp\n')
notes = paste0(notes, round(mean(move$move, na.rm = T), 2),
	" bp of average position correction for reverse strand alignments.\n")

# Apply CIGAR-based shift
mt$pos[mt$flag %in% rc] <- mt$pos[mt$flag %in% rc] + move$move[
	match(mt$cigar[mt$flag %in% rc], move$cigar)]

if ( !no_cutsite ) {
	# Apply CS-based shift
	mt$pos[!mt$flag %in% rc] <- mt$pos[!mt$flag %in% rc] - nchar(cutsite)
}

# Output 2 ---------------------------------------------------------------------

cat(' · Writing output 3\' position ...\n')
write.table(mt[, c('chr', 'pos', 'mapq', 'lseq', 'lqual', 'qname')],
	paste0(condition, '.filtered.umi.pos.txt'),
	row.names = F, col.names = F, quote = F, sep = '\t')

# Save notes -------------------------------------------------------------------

write.table(notes, paste0(condition, '.sam_filter_notes.txt'),
	quote = F, col.names = F, row.names = F, append = T)

# END --------------------------------------------------------------------------

################################################################################
