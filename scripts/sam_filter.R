#!/usr/bin/Rscript

# ------------------------------------------------------------------------------
#
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 1.0.0
# Description: script to filter alignment output (SAM file)
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
t <- read_delim(paste0(condition, '.linkers.sam'), '\t',
	col_names = c('qname', 'flag', 'chr', 'pos', 'mapq', 'cigar', 'rnext',
		'pnext', 'tlen', 'lseq', 'seq', 'lqual', 'qual'),
	col_types = 'ciciicciicccc')

# Binarize flags ---------------------------------------------------------------

cat(paste0(' · Binarizing flags.\n'))

# Unique flags
uf <- sort(unique(t$flag))

# Unique binary flags
ubf <- unlist(mclapply(uf,
	FUN = function(x) {
		rev(substr(paste(rev(as.integer(intToBits(x))), collapse = ""), 21, 32))
	}
	, mc.cores = num_proc
))

# Remove supplementary alignments ----------------------------------------------

cat(paste0(' · Removing secondary alignments...\n'))

# Identify secondary alignments from the flag
sec <- unlist(mclapply(ubf,
	FUN = function(x) { as.numeric(substr(x, 1, 1)) == 0 }
	, mc.cores = num_proc))

# Keep only primary alignments
mt <- t[t$flag %in% uf[sec],]

cat(paste0(' >>> Kept ', nrow(mt), ' rows out of ', nrow(t),
	' (', round(nrow(mt) / nrow(t) * 100, 2), '%)\n'))
notes = paste0(notes, nrow(t) - nrow(mt), " secondary alignments.\n")

# Filter chimeras --------------------------------------------------------------

# Update unique binary flags
uf2 <- sort(unique(mt$flag))
ubf2 <- ubf[which(uf %in% uf2)]

# Check for R2 presence from the flags
r2 <- unlist(mclapply(ubf2,
	FUN = function(x) { as.numeric(substr(x, 5, 5)) }
	, mc.cores = num_proc))

if ( 0 != length(which(1 == r2)) ) {
	cat(paste0(' · Removing chimeras...\n'))
	
	# Identify chimeras (R1 and R2 on different chromosomes)
	tmp <- unique(data.frame(qname = mt$qname, chr = mt$chr))
	tmp <- table(tmp$qname)
	chimeras <- names(tmp)[tmp > 1]
	names(chimeras) <- NULL

	if ( 0 != length(chimeras) ) {
		# Remove chimeras
		old_nrow = nrow(mt)
		mt <- mt[-which(mt$qname %in% chimeras),]

		cat(paste0(' >>> Kept ', nrow(mt), ' rows out of ', old_nrow,
			' (', round(nrow(mt) / old_nrow * 100, 2), '%)\n'))
		notes = paste0(notes, old_nrow - nrow(mt), " chimeric reads.\n")
	} else {
		cat(paste0(' >>> No chimeras found.\n'))
		notes = paste0(notes, 0, " chimeric reads.\n")
	}
} else {
	cat(paste0(' >>> No R2 found, thus no chimeras.\n'))
	notes = paste0(notes, 0, " chimeric reads.\n")
}

# Filter unmapped reads --------------------------------------------------------

cat(paste0(' · Removing unmapped reads.\n'))

old_nrow = nrow(mt)
mt <- mt[mt$chr %in% c(1:22, 'X', 'Y'),]

cat(paste0(' >>> Kept ', nrow(mt), ' rows out of ', old_nrow,
	' (', round(nrow(mt) / old_nrow * 100, 2), '%)\n'))
notes = paste0(notes, old_nrow - nrow(mt), " unmapped reads.\n")

# Filter R2s -------------------------------------------------------------------

# Update binary flags
uf2 <- sort(unique(mt$flag))
ubf2 <- ubf[which(uf %in% uf2)]

# Identify R2s from flag
r2 <- unlist(mclapply(ubf2,
	FUN = function(x) { as.numeric(substr(x, 5, 5)) }
	, mc.cores = num_proc))

cat(paste0(' · Removing R2...\n'))
if ( 0 != length(which(1 == r2)) ) {

	old_nrow = nrow(mt)
	mt <- mt[mt$flag %in% uf[r2 == 0],]

	cat(paste0(' >>> Kept ', nrow(mt), ' out of ', old_nrow,
		' (', round(nrow(mt) / old_nrow * 100, 2), '%)\n'))
	notes = paste0(notes, old_nrow - nrow(mt), " R2 reads.\n")
} else {
	cat(paste0(' >>> No R2 found...\n'))
	notes = paste0(notes, 0, " R2 reads.\n")
}

# Filter MAPQ ------------------------------------------------------------------

cat(paste0(' · Filtering MAPQ...\n'))

old_nrow = nrow(mt)
mt <- mt[mt$mapq >= mapq_thr,]

cat(paste0(' >>> MAPQ threshold [', mapq_thr, ']\n'))
cat(paste0(' >>> Selected ', nrow(mt), ' out of ', old_nrow,
	' (', round(nrow(mt) / old_nrow * 100, 2), '%)\n'))
notes = paste0(notes, old_nrow - nrow(mt),
	" reads with MAPQ < ", mapq_thr, ".\n")

# Filtering chromosomes --------------------------------------------------------

cat(paste0(' · Removing chromosomes...\n'))
if ( 0 != length(rmChr) ) {

	# Removed reads counter
	nrm = 0

	for (chr in rmChr) {
		cat(paste0(' >>> Removing chr', chr, '...\n'))

		# Current count of removed reads
		cnrm = sum(chr == mt$chr)
		cat(paste0(' >>> Removing ', cnrm, ' reads...\n'))

		# Remove reads
		if ( 0 != cnrm )
			mt = mt[-which(chr == mt$chr),]

		# Update general count
		nrm = nrm + cnrm
	}

	notes=paste0(notes, paste0(nrm, ' reads from removed chromosomes (',
		paste(paste0('chr', rmChr), collapse = ', '), ').\n'))
} else {
	cat(paste0(' >>> No chromosomes to be removed.\n'))
	notes=paste0(notes, paste0(0, ' reads from removed chromosomes (',
		paste(paste0('chr', rmChr), collapse = ', '), ').\n'))
}


# Output 1 ---------------------------------------------------------------------

#cat(' · Writing output ...\n')
notes = paste0(notes, nrow(mt), " reads left after filtering.\n")
write.table(mt, paste0(condition, '.filtered.sam'),
	row.names = F, col.names = F, quote = F, sep = '\t')
#save(mt, file = paste0(condition, '.filtered.sam.RData'))

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
	quote = F, col.names = F, row.names = F)

# END --------------------------------------------------------------------------

################################################################################
