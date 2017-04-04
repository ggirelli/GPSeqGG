#!/usr/bin/Rscript

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 1.0.0
# Description: Shuffle a certain percentage of reads in a bed file.
# 
# ------------------------------------------------------------------------------



# DEPENDENCIES =================================================================

suppressMessages(library(argparser))
suppressMessages(library(data.table))
suppressMessages(library(parallel))
suppressMessages(library(readr))

# INPUT ========================================================================

# Create arguent parser
parser = arg_parser('Shuffle bed file reads.', name = 'bed_shuffle.R')

# Define mandatory arguments
parser = add_argument(parser, arg = 'seed', short = '-s',
	help = 'Seed for random number generation.', type = class(0))
parser = add_argument(parser, arg = 'bedfile', short = '-b',
	help = 'path to bedfile.', type = class(''))

# Define elective arguments
parser = add_argument(parser, arg = '--nIter', short = '-n',
	help = 'Number of iterations.',
	default = 100, nargs = 1, type = class(0))
parser = add_argument(parser, arg = '--perc', short = '-p',
	help = 'Percentage of reads to shuffle.',
	default = 10, nargs = 1, type = class(0))
parser = add_argument(parser, arg = '--outDir', short = '-o',
	help = 'Output directory.',
	default = "./shuffled/", nargs = 1, type = class(''))
parser = add_argument(parser, arg = '--keepSeed', short = '-k',
	help = 'Reload previous seed state. Use on subsequent runs.',
	default = TRUE, nargs = 1, type = class(TRUE))
parser = add_argument(parser, arg = '--threads', short = '-t',
	help = 'Number of threads for parallelization.',
	default = 1, nargs = 1, type = class(0))

# Parse arguments
p = parse_args(parser)

# Attach argument values to variables
attach(p['' != names(p)])

fname = unlist(strsplit(bedfile, '/', fixed = TRUE))
fname = fname[length(fname)]
fname = unlist(strsplit(fname, '.', fixed = TRUE))
fname = fname[-length(fname)]
fname = paste(fname, collapse = '.')

# FUNCTIONS ====================================================================

# RUN ==========================================================================

# Log info
cat(paste0(' · Shuffling x', nIter, ' ', perc, '% of ', bedfile, '\n'))

# Set seed
set.seed(seed)

# Load seed if available
fseed = paste0(outDir, '/.Random.seed.RData')
if ( file.exists(fseed) && keepSeed ) {
	cat(paste0(' >> Loading previous seed status...\n'))
	load(fseed)
}

# Read bedfile -----------------------------------------------------------------

bedcols = c('chr', 'start', 'end', 'name', 'score')
bed = read_delim(bedfile, '\t', skip = 1,
	col_names = bedcols, col_types = 'ciici')

# Shuffle ----------------------------------------------------------------------

# Count reads
nreads = sum(bed$score)
toShuffle = round(nreads * perc / 100)
cat(paste0(' >> Found ', nreads, ' reads.\n'))

cat(paste0(' >>> Pre-shuffling...\n'))
preshuffle = factor(unlist(mclapply(seq(nrow(bed)),
	FUN = function(nr) {
		row = bed[nr,]
		return(rep(gsub(' ', '', paste(row[1:4], collapse = '~')), row[5]))
	}
	, mc.cores = threads
)))
cat(paste0(' >>> Identifying locations...\n'))
positions = levels(preshuffle)

l = lapply(seq(nIter),
	function(i) {
		cat(paste0(' >>>> Iteration #', i, '\n'))

		shuffled = preshuffle

		cat(paste0(' >>>> Identifying origins...\n'))
		from = sample(seq(length(preshuffle)), toShuffle)
		cat(paste0(' >>>> Identifying destinations...\n'))
		to = sample(seq(length(positions)), toShuffle, replace = TRUE)

		cat(paste0(' >>>> Moving reads...\n'))
		shuffled[from] = positions[to]

		cat(paste0(' >>>> Counting shuffled reads...\n'))
		counts = table(shuffled)

		cat(paste0(' >>>> Preparing output...\n'))
		out = rbindlist(mclapply(seq(length(counts)),
			FUN = function(j) {
				loc = names(counts)[j]
				loc = unlist(strsplit(loc, '~', fixed = T))
				data = data.frame(
					chr = loc[1],
					start = loc[2],
					end = loc[3],
					name = loc[4],
					score = counts[j]
				)
				return(data)
			}
			, mc.cores = threads
		))
		#cat(paste0(' >> ', sum(out$score), ' reads in the shuffled file.\n'))

		outfname = paste0(outDir, '/', fname, '.iter', i, '.', perc, 'perc.bed')
		write.table(out, outfname, row.names = F, quote = F, sep = '\t')
	}
)

# Save current seed status
cat(paste0(' >> Saving seed status...\n'))
save('.Random.seed', file=paste0(outDir, '/.Random.seed.RData'))

# END --------------------------------------------------------------------------

################################################################################
