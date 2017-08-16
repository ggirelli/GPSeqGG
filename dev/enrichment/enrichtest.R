#!/usr/bin/Rscript

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 0.1.1
# Description: enrichment analysis of GPSeq data.
# Notes: based on regioneR [http://goo.gl/01Wjvc] with statistics from
# 		 Phipson&Smith, 2010 [http://goo.gl/ZHnvY7].
# 
# @TODO:
# 	- Make permutation generation more efficient (time and memory).
# 	- Non-ranked data-tracks compatibility.
# 	- No-cutsite list scenario.
# 
# ------------------------------------------------------------------------------



# DEPENDENCIES =================================================================


suppressMessages(library(argparser))
suppressMessages(library(data.table))
suppressMessages(library(GenomicRanges))
suppressMessages(library(parallel))
suppressMessages(library(readr))
suppressMessages(library(regioneR))

# INPUT ========================================================================

# Manual arguments -------------------------------------------------------------
# help = FALSE
# reuse_perm = TRUE
# opts = NA
# cutsite_range = 20
# min_perm = 500
# max_perm = 2000
# num_proc = 4
# seed = 1657984
# palette = "#74C957,#0068CB,#CC2343,#009688,#FF9800"
# dirpath = "/home/gire/Desktop/BiCro-Analysis/Sequencing/TK22/outdata/"
# experiment = "TK22"
# conditions = 'neg,2min,30min,on'
# data_track = paste0('/media/gire/Data2/BiCro-Resources/Data_tracks_byQ/',
# 	'Timing/0wgEncodeUwRepliSeqHelas3AllAlnRep1.bed')
# data_track_field = "phase"
# label = 'timing'
# cutsites = "/media/gire/Data2/BiCro-Resources/hg19.HindIII.txt"
# mask_file = paste0('/media/gire/Data2/BiCro-Data/Sequencing/TK22/',
# 	'consensusBlacklist.centro.telo.HeLaCNV.tsv')
# plot = T
# sig_thr = .05
# notsig_color = '#FFFFFF'
# dep_color = '#CC2343'
# enr_color = '#009688'

# Create arguent parser
parser = arg_parser('Perform enrichment analysis.
The scripts expect dirpath to contain a folder for every condition, whose names
are provided as a comma-separated values string.', name = 'enrichtest')

# Dataset arguments ------------------------------------------------------------
parser = add_argument(parser, arg = 'dirpath',
	help = 'Directory to which the binned cutsites will be saved.')
parser = add_argument(parser, arg = 'experiment',
	help = 'Experiment ID (e.g., TK26).')
parser = add_argument(parser, arg = 'conditions',
	help = 'Comma-separated conditions (e.g., 400U2h,200U1h).')

# Data track arguments ---------------------------------------------------------
parser = add_argument(parser, arg = 'data-track',
	help = 'Data-track file. Structure: (chr|start|end|...) with header.')
parser = add_argument(parser, arg = 'data-track-field',
	help = 'Data-track column to be tested for enrichment.')

parser = add_argument(parser, arg = 'cutsites',
	help = paste0('File containing the cutsite positions. ',
	'Structure: (chr|pos) with no header.'))
parser = add_argument(parser, arg = '--cutsite-range', short = '-u',
	default = 20, help = 'Cutsite window half-size.')

# Mask arguments ---------------------------------------------------------------
parser = add_argument(parser, arg = '--mask-file', default = '',
	help = paste0('File with regions to be masked. ',
	'Structure: (id|chr|start|end|status) with no header.'))

# Permutation arguments --------------------------------------------------------
parser = add_argument(parser, arg = '--min-perm', short = '-i',
	help = 'Minimum number of permutations to be used.',
	default = 500, type = class(0))
parser = add_argument(parser, arg = '--max-perm', short = '-a',
	help = 'Maximum number of permutations to be used.',
	default = 2e3, type = class(0))
parser = add_argument(parser, arg = '--reuse-perm', short = '-r',
	help = 'Use previous permutations if available? (y/n)', flag = T)

# General arguments ------------------------------------------------------------
parser = add_argument(parser, arg = '--num-proc', short = '-c',
	help = 'Number of cores for parallel computation.',
	default = 1, type = class(0))
parser = add_argument(parser, arg = '--seed', default = 1657984, short = '-s',
	type = class(0), help = 'Seed for random number generator.')
parser = add_argument(parser, arg = '--plot', short = '-p',
	flag = T, help = 'Enable plot to enrichtest.pdf')
parser = add_argument(parser, arg = '--palette',
	default = '#74C957,#0068CB,#CC2343,#009688,#FF9800',
	help = 'Comma-separated hexadecimal color codes.')
parser = add_argument(parser, arg = '--label', short = '-l',
	default = '', help = 'Data track label for output.')

# Plot arguments ---------------------------------------------------------------
parser = add_argument(parser, arg = '--sig-thr', short = '-t',
	default = .05, help = 'Significance threshold for p-value interpretation.')
parser = add_argument(parser, arg = '--dep-color',
	default = '#CC2343', help = 'Significant depletion color label.')
parser = add_argument(parser, arg = '--enr-color',
	default = '#009688', help = 'Significant enrichment color label.')
parser = add_argument(parser, arg = '--notsig-color',
	default = '#FFFFFF', help = 'Not significant color label.')

# Parse arguments --------------------------------------------------------------
p = parse_args(parser)

# Attach argument values to variables
attach(p['' != names(p)])

# Split comma-separated variables
conditions = unlist(strsplit(conditions, ',', fixed = T))

# Add trailing slash to paths
if ( '/'!= substr(dirpath, nchar(dirpath), nchar(dirpath)) )
	dirpath = paste0(dirpath, '/')

# Fix input strings
label = tolower(gsub(' ', '_', label, fixed = T))

# Fix colors
sigcol = data.frame(
	not_sig = notsig_color,
	dep_sig = dep_color,
	enr_sig = enr_color,
	stringsAsFactors = F
)

# FUNCTIONS ====================================================================

get_condition_grange = function(condition) {
	# Retreive the data of the specified condition.
	# 
	# Args:
	# 	condition (string): condition folder name
	# 
	# Returns:
	# 	GRange object with condition UMI at cutsite.
	# 

	# Read UMIpos.unique file
	umic <- read.delim( paste0(dirpath, condition, '/UMIpos.unique.atcs.txt'),
		as.is = T, header = F)
	colnames(umic) <- c('id', 'pos', 'seq')

	# Update format ------------------------------------------------------------
	umic$chr <- umic$id
	umic$chr[umic$id == 23] <- 'X'
	umic$chr[umic$id == 24] <- 'Y'

	# Add chr-name
	umic$chr <- paste0('chr', umic$chr)

	# Make range around CS
	umic$start <- umic$pos - cutsite_range
	umic$end <- umic$pos + cutsite_range

	# Add UMI count ------------------------------------------------------------
	umic$count <- unlist(lapply(strsplit(umic$seq, ' '), FUN = length))

	# Add UMI probability ------------------------------------------------------
	umic$p <- umic$count / sum(umic$count)

	# Output -------------------------------------------------------------------

	# Re-order columns
	umic <- umic[, c(4, 7, 8, 5:6, 1:3)]
	return(GRanges(umic))
}

mask_region_set = function(A, M,
	condition = '', verbose = 0, num_proc = 1, ...) {
	# Remove regions or portions of regions that overlap with regions in M.
	# 
	# Args:
	# 	A (GRange): the input region set
	# 	M (GRange): the mask set
	# 	condition (string): condition description
	# 	num_proc (int): number of thread for parallelization
	# 
	# Returns:
	# 	A without regions or portions of regions overlapping with M.
	# 
	
	# Dependencies -------------------------------------------------------------
	library(GenomicRanges)
	library(parallel)

	# Format -------------------------------------------------------------------
	A <- GRanges(A)
	M <- GRanges(M)

	# Keep only common seqnames
	com_seqnames = intersect(levels(seqnames(A)), levels(seqnames(M)))
	A <- keepSeqlevels(A, com_seqnames)
	M <- keepSeqlevels(M, com_seqnames)

	# Remove 'within' overlaps -------------------------------------------------
	
	ovls_within <- as.data.frame(findOverlaps(A, M, type = 'within'))	
	if ( 0 != nrow(ovls_within) )
		A <- A[-ovls_within$queryHits,]
	
	# Log
	if ( verbose > 0 )
		cat(paste0(' [', condition, '] Removed ', nrow(ovls_within),
			' "within" overlaps.\n'))

	# Check other types of overlaps --------------------------------------------
	ovls_any <- as.data.frame(findOverlaps(A, M, type = 'any'))

	# Remove already filtered 'within' overlaps
	ovls_any <- ovls_any[-which(ovls_any$queryHits %in% ovls_within$queryHits),]

	# Remove region portions
	L <- mclapply(unique(ovls_any$subjectHits),
		FUN=function(i) {

			# Get masked region borders
			start <- start(M)[i]
			end <- end(M)[i]

			# Identify overlaps
			queryHits <- ovls_any$queryHits[which(ovls_any$subjectHits == i)]

			# Update overlap borders
			start(A)[queryHits][which(start(A)[queryHits] >= start)] = start - 1
			end(A)[queryHits][which(end(A)[queryHits] <= end)] = end + 1

			# Empty output
			return(NULL)
		}
		, mc.cores = num_proc
	)

	# Log
	if ( verbose > 0 )
		cat(paste0(' [', condition, '] Partially masked ', nrow(ovls_any),
			' "non-within" overlaps.\n'))

	# Masked Region Set
	return(A)
}

gen_perm = function(A, universe, field, ...) {
	# Permutate maintaining the original UMI distribution.
	#
	# Args:
	# 	A (GRange): the input region set
	# 	universe (GRange): A is a subset of universe
	# 	field (string): the name of the column in A with the field of interest
	# 
	# Returns:
	# 	A sample of the universe keeping the UMI counts constant.
	# 	In other words, the number of times that the cutsites are cut is
	# 	maintained as in the provided A set (the overall distribution).
	# 
	
	# Check format
	A <- toDataframe(A)
	universe <- toDataframe(universe)

	# Sample the universe
	B <- universe[sample(1:nrow(universe), nrow(A)),]

	# Add permutated UMI counts
	B[, field] <- A[sample(1:nrow(A), nrow(A)), field]

	# Reorder rows
	B <- B[order(B$start),]
	B <- B[order(B$chr),]
	rownames(B) <- 1:nrow(B)

	# Output
	return(GRanges(B))
}

eval_region_set = function(A, B, field, ...) {
	# Count the unique overlaps of a region set
	# extended based on the number of UMIs.
	# 
	# 
	# Args:
	# 	A (GRange): the input region set
	# 	B (GRange): the region set to be compared (e.g., data-track)
	# 	field (string): the name of the column in A with the field of interest
	# 
	# Returns:
	# 	The total number of overlaps.
	# 

	# DEPENDENCIES -------------------------------------------------------------
	
	library(GenomicRanges)
	library(regioneR)

	# FORMAT -------------------------------------------------------------------
	
	A <- toDataframe(A)

	# RUN ----------------------------------------------------------------------

	# Duplicate based on umi counts
	extA <- GRanges(A[rep(1:nrow(A), A[, field]),])

	# Return uniqued overlaps with extended A set
	return(sum(countOverlaps(extA, B) > 1))
}

enrichment_test = function(
	A, U,
	Bs, B_field, B_level,
	condition = '',
	A_dup_field = 'count',
	min_perm = 500,
	max_perm = 2000,
	rand_set = NULL,
	sd_thr = 1e-4,
	sd_window = 100,
	num_proc = 1,
	plot = F
) {
	# Run condition-specific rank-specific enrichment test.
	# 
	# Args:
	# 	A (GRange): input region set
	# 	A_dup_field (string): mcols(A) header with the duplicate counts
	# 	U (GRange): universe (A is a subset of U)
	# 	Bs (GRange): ranked data track
	# 	B_field (string): data track rank column
	# 	B_level (string): rank to test for enrichment
	# 	condition (string): condition description for log
	# 	min_perm (int): minimum number of permutations
	# 	max_perm (int): maximum number of permutations
	# 	rand_set (list): list of permutations of A (if available)
	# 	sd_thr (float): sd threshold for the selection of the number of perms
	# 	sd_window (int): upper limit of range of number of perms to calculate sd
	# 	num_proc (int): number of clusters for parallelization
	# 	plot (bool): whether to plot
	# 
	# Returns:
	# 	List containing enrichment test results.
	# 

	# Subset Bs based on specified field and level -----------------------------
	B = Bs[mcols(Bs)[, B_field] == B_level]

	# Merge overlapping regions
	B = mergeRegions(B, B)

	# Use previous permutations ------------------------------------------------
	if ( !is.null(rand_set) ) {
		prev_iter = length(rand_set)
		cat(paste0(' · Found ', prev_iter, ' permutations ...\n'))

		# Generate more permutations if not enough
		mperms = min_perm - prev_iter
		if ( mperms > 0 ) {
			cat(paste0(' [', condition, '] ',
				'Generating more permutations [', mperms, '] ...\n'))

			pb = txtProgressBar(min = 0, max = mperms, style = 3)
			rand_set = c(rand_set, mclapply(1:mperms,
				FUN = function(i) {
					#cat(paste0(' · Generating permutation #', i, ' ...\n'))
					setTxtProgressBar(pb, i)
					gen_perm(A, U, A_dup_field)
				}
				, mc.cores = num_proc
			))
			setTxtProgressBar(pb, mperms)
			cat('\n')
		}
	} else {
		cat(paste0(' [', condition, '] ',
			'Generating permutations [',min_perm, '] ...\n'))
		prev_iter = 0

		pb = txtProgressBar(min = 0, max = min_perm, style = 3)
		rand_set = mclapply(1:min_perm,
			FUN = function(i) {
				# cat(paste0(' [', condition, '] ',
				# 	'Generating permutation #', i, ' ...\n'))
				setTxtProgressBar(pb, i)
				gen_perm(A, U, A_dup_field)
			}
			, mc.cores = num_proc
		)
		setTxtProgressBar(pb, min_perm)
		cat('\n')
	}

	# <--- ADD HERE AUTOMATIC NUMBER OF PERMUTATION CALCULATION !!!

	# Save current permutation set if bigger than the input --------------------
	if ( prev_iter < length(rand_set) ) {
		rand_set_path = paste0(dirpath, condition, '/',
			experiment, '.', condition, '.enrichment_analysis.rand_set.RData')
		cat(paste0(' >>> Dumping permutations ...\n'))
		cat(paste0(' >>> ', rand_set_path, '\n'))
		save('rand_set', file = rand_set_path)
	}

	# Evaluate permutations ----------------------------------------------------
	cat(' · Evaluating permutations ...\n')
	pb = txtProgressBar(min = 0, max = length(rand_set), style = 3)
	rand_ev = unlist(mclapply(1:length(rand_set),
		FUN = function(i) {
			setTxtProgressBar(pb, i)
			eval_region_set(rand_set[[i]], B, A_dup_field)
		}
		, mc.cores = num_proc
	))
	setTxtProgressBar(pb, length(rand_set))
	cat('\n')

	# Evaluate original --------------------------------------------------------
	orig_ev = eval_region_set(A, B, A_dup_field)

	# Calculate Z-score and p-value --------------------------------------------
	ntimes <- length(rand_set)
	num_nas <- length(which(is.na(rand_ev)))
	num_valid_values <- ntimes - num_nas
	rand_mean <- mean(rand_ev, na.rm = T)

	alt = ifelse(orig_ev < rand_mean, "less", "greater")

	if ( "less" == alt ) {
		pval = (sum(orig_ev > rand_ev, na.rm = TRUE) + 1)
		pval = pval / (num_valid_values + 1)
	} else {
		pval = (sum(orig_ev < rand_ev, na.rm = TRUE) + 1)
		pval = pval / (num_valid_values + 1)
	}

	if (orig_ev == 0 & all(rand_ev == 0)) {
		pval = 1
		zscore = NA
	} else {
		zscore = (orig_ev - mean(rand_ev, na.rm = TRUE))
		zscore = zscore  / sd(rand_ev, na.rm = TRUE)
		zscore = round(zscore, 4)
	}

	# Plot ---------------------------------------------------------------------
	
	if ( plot ) {
		# Plot histogram
		xlim = c(min(c(orig_ev, rand_ev)), max(c(orig_ev, rand_ev)))
		h = hist(rand_ev, breaks = length(rand_set), xlim = xlim,
			xlab = 'numOverlaps', ylab = '', main = '', prob = T,
			col = 'gray50', border = 'gray50')
		lines(density(rand_ev))

		# Add random mean and sample value
		abline(v = orig_ev, col = 4, lwd = 3)
		abline(v = rand_mean, col = 3, lwd = 3)

		# Add header
		mtext(paste0('[', experiment, '~', condition, ' | ', B_level,
			' | ', alt, '] green: perms; blue: orig'))
		mtext(paste0('n_perm: ', length(rand_set)), line = 1)
		mtext(paste0('Z-score: ', round(zscore, 2)), line = 2)
		mtext(paste0('p-value: ', round(pval, 6)), line = 3)

		# Add alternative region
		if ( 'less' == alt ) {
			q = qnorm(0.05, rand_mean, sd(rand_ev))
			rect(xleft = xlim[1], xright = q, ybottom = 0, ytop = max(h$counts),
				col = rgb(1, 0, 0, .2), border = 0)
		} else {
			q = qnorm(0.95, rand_mean, sd(rand_ev))
			rect(xleft = q, xright = xlim[2], ybottom = 0, ytop = max(h$counts),
				col = rgb(1, 0, 0, .2), border = 0)
		}
	}


	# Output -------------------------------------------------------------------

	return(list(
		level = B_level,
		rand_ev = rand_ev,
		rand_mean = rand_mean,
		ntimes = ntimes,
		num_nas = num_nas,
		num_valid_values = num_valid_values,
		orig_ev = orig_ev,
		alt = alt,
		pval = pval,
		zscore = zscore
	))
}

# INPUT ========================================================================

# Condition data ---------------------------------------------------------------

cat('Reading and formatting condition data...\n')
As = mclapply(conditions, FUN = get_condition_grange, mc.cores = num_proc)
names(As) = conditions

# Mask file --------------------------------------------------------------------

cat('Reading and formatting mask file...\n')
if ( 0 != nchar(mask_file) ) {
	M = suppressMessages(read_delim(mask_file, '\t',
		col_names = c('id', 'chr', 'start', 'end', 'status'), progress = F))
	M = M[, c(2:4, 1, 5)]
	M = GRanges(M)
} else {
	seqnams = unique(unlist(lapply(As, FUN = function(x) levels(seqnames(x)))))
	M = GRanges(paste0(seqnams, ':0-0.'))
}

# Mask condition data ----------------------------------------------------------

cat('Masking condition data...\n')
As = lapply(1:length(As),
	FUN = function(i) {
		mask_region_set(As[[i]], M, names(As)[i],
			verbose = 1, num_proc = num_proc)
	}
)

# Universal Region Set (+- range nt around known cutsite) ----------------------

cat('Reading and formatting universal set...\n')
U = suppressMessages(read_delim(cutsites, '\t',
	col_names = c('chr', 'pos'), progress = F))

# Make range around CS
U$start = U$pos - cutsite_range
U$end = U$pos + cutsite_range

# Select columns and make into a GRange object
U = U[, c('chr', 'start', 'end')]
U = GRanges(U)

# Apply mask
U = mask_region_set(U, M, 'universe', verbose = 1, num_proc = num_proc)

# Data-track -------------------------------------------------------------------

cat('Reading and formatting datatrack...\n')
Bs <- suppressMessages(read_delim(data_track, '\t', progress = F))
colnames(Bs)[1:3] = c('chr', 'start', 'end')

# Make into a GRance object
Bs <- GRanges(Bs)

# Apply mask
Bs <- mask_region_set(Bs, M, 'datatrack', verbose = 1, num_proc = num_proc)

# Keep only common seqnames
A_levels = unique(unlist(lapply(As, FUN = function(A) levels(seqnames(A)))))
com_seqnames = intersect(A_levels, levels(seqnames(Bs)))
As <- lapply(As,
	FUN = function(A, com_seqnames) {
		keepSeqlevels(A, intersect(levels(seqnames(A)), com_seqnames))
	}, com_seqnames
)
Bs <- keepSeqlevels(Bs, com_seqnames)

# RUN ==========================================================================

# Select data track ranks
B_levels <- unlist(unique(mcols(Bs)[data_track_field]))
names(B_levels) <- NULL

# Test enrichment for each condition (one at a time) ---------------------------
if ( plot ) pdf(paste0(experiment, '.', label, '.enrichtest.pdf'))
etr <- lapply(1:length(As), FUN = function(cond_id) {
	cat(paste0('  [ ', conditions[cond_id], ' ]\n'))

	rand_set_path = paste0(dirpath, conditions[cond_id], '/', experiment,
		'.', conditions[cond_id], '.enrichment_analysis.rand_set.RData')
	if ( reuse_perm ) {
		if ( file.exists(rand_set_path) ) {
			# Load previous permutations
			cat(paste0(' · Loading previous permutations ...\n'))
			load(rand_set_path)
		} else {
			# New permutations
			cat(paste0(' [', conditions[cond_id], '] ',
				'Generating permutations [', min_perm, '] ...\n'))
			prev_iter = 0

			# Generate
			pb = txtProgressBar(min = 0, max = min_perm, style = 3)
			rand_set = mclapply(1:min_perm,
				FUN = function(i) {
					setTxtProgressBar(pb, i)
					gen_perm(As[[cond_id]], U, 'count')
				}
				, mc.cores = num_proc
			)
			setTxtProgressBar(pb, min_perm)
			cat('\n')

			# Save
			cat(paste0(' >>> Dumping permutations ...\n'))
			cat(paste0(' >>> ', rand_set_path, '\n'))
			save('rand_set', file = rand_set_path)
		}
	} else { rand_set = NULL }


	lapply(B_levels, FUN = function(B_level, cond_id) {
		cat(paste0('  > ', B_level, ' <\n'))

		# Run enrichment test
		enrichment_test(As[[cond_id]], U, Bs, data_track_field, B_level,
			conditions[cond_id], rand_set = rand_set, num_proc = num_proc,
			plot = plot)

	}, cond_id)
})
if ( plot ) graphics.off()

# Output -----------------------------------------------------------------------

# Save output as RData
save(etr, file = paste0(experiment, '.', label, '.enrichtest.RData'))

# Save output as table
outdf = rbindlist(lapply(1:length(etr),
	FUN=function(cond_id) {
		rbindlist(lapply(1:length(etr[[cond_id]]),
			FUN = function(lev_id) {
				B_level = B_levels[lev_id]
				condition = conditions[cond_id]
				ltmp = etr[[cond_id]][[lev_id]]
				
				return(data.frame(
					condition = conditions[cond_id],
					level = B_level,
					ntimes = ltmp$ntimes,
					num_valid_values = ltmp$num_valid_values,
					rand_mean = ltmp$rand_mean,
					rand_sd = sd(ltmp$rand_ev),
					orig_ev = ltmp$orig_ev,
					zscore = ltmp$zscore,
					alt = ltmp$alt,
					pval = ltmp$pval
				))

			}
		))
	}
))
write.table(outdf, paste0(experiment, '.', label, '.enrichtest.txt'),
	quote = F, sep = '\t', row.names = F)

# Other plots ------------------------------------------------------------------

if ( plot ) {

	# A barplot per condition, with data-track levels on the X
	pdf(paste0(experiment, '.', label, '.enrichtest.per_condition.pdf'))
	l = by(outdf, outdf$condition,
		FUN = function(subdf) {

			# Set colors
			color = rep(0, nrow(subdf))
			color[subdf$pval <= sig_thr & subdf$zscore >= 0] = sigcol$enr_sig
			color[subdf$pval <= sig_thr & subdf$zscore < 0] = sigcol$dep_sig
			color[subdf$pval > sig_thr] = sigcol$not_sig

			# Plot
			barplot(subdf$zscore, names = subdf$level,
				col = color, ylab = 'z-score', xlab = 'Data-track value',
				main = paste0(experiment, ' ~ ', unique(subdf$condition))
			)

		}
	)
	graphics.off()

	# A barplot per data-track level, with conditions on the X
	pdf(paste0(experiment, '.', label, '.enrichtest.per_level.pdf'))
	l = by(outdf, outdf$level,
		FUN = function(subdf) {

			# Set colors
			color = rep(0, nrow(subdf))
			color[subdf$pval <= sig_thr & subdf$zscore >= 0] = sigcol$enr_sig
			color[subdf$pval <= sig_thr & subdf$zscore < 0] = sigcol$dep_sig
			color[subdf$pval > sig_thr] = sigcol$not_sig

			# Plot
			barplot(subdf$zscore, names = subdf$condition,
				col = color, ylab = 'z-score', xlab = 'condition',
				main = paste0(experiment, '~', unique(subdf$level))
			)

		}
	)
	graphics.off()

}

# END ==========================================================================

################################################################################
