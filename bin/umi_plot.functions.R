#!/usr/bin/Rscript

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 0.1.1
# Description: functions used by the `umi_plot` script
# 
# ------------------------------------------------------------------------------



# DEPENDENCIES =================================================================

suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
suppressMessages(library(grid))
suppressMessages(library(gridExtra))
suppressMessages(library(parallel))
suppressMessages(library(readr))
suppressMessages(library(reshape2))

# FUNCTIONS ====================================================================

plotUMIdistrib = function(
	dirpath, conditions, experiment, cutsites,
	ncores = 1, suff = ''
) {
	# Plot the distribution of UMIs across conditions.
	# 
	# Args:
	# 	dirpath (string): base directory path (outdata)
	# 	conditions (list): list of condition-specific data-frames
	# 	experiment (string): experiment ID
	# 	cutsites (string): path to file with cutsite list
	# 	ncores (int): number of threads for parallelization
	# 	suff (string): suffix used for input/output operations
	# 	
	
	# Calculate the cumulative umi distribution for every condition
	cum_umi_count <- lapply(conditions,
		FUN=function(condition) {

			# Read unique UMI file
			fname <- 'UMIpos.unique'
			if ( file.exists(cutsites) ) fname <- paste0(fname, '.atcs')
			fname = paste0(dirpath, condition, '/', fname, suff, '.txt')
			t <- read.delim(fname, header = F)
			colnames(t) <- c('chr', 'pos', 'seq')

			# Count UMIs at every cutsite
			t$count <- unlist(mclapply(strsplit(as.character(t$seq), ' '),
				FUN = length, mc.cores = ncores))

			# Prepare UMI frequency table
			t <- table(t$count)
			k <- as.numeric(names(t))

			# Make it cumulative
			unlist(lapply(1:max(k),
				FUN = function(i) {
					sum(t[which(k >= i)])
				}
			))
		}
	)

	# Save the data and plot
	save('cum_umi_count',
		file = paste0(dirpath, experiment, '.umi_distribution.RData'))
	pdf(paste0(dirpath, "plots/", experiment, '.umi_distribution', suff, '.pdf'),
		width = 10, height = 10)

	# Prepare X/Y-indexes
	xvals <- 1:max(unlist(lapply(cum_umi_count, length)))
	yvals <- as.numeric(unique(unlist(cum_umi_count)))

	# PLOT ---------------------------------------------------------------------
	p <- ggplot()
	p <- p + lapply(1:length(cum_umi_count), FUN=function(i) {
		data <- data.frame(
			x = 1:length(cum_umi_count[[i]]),
			y = cum_umi_count[[i]],
			col = conditions[i]
		)
		geom_line(data = data, aes(x = x, y = y, colour = col))
	})
	p <- p + scale_colour_manual("", breaks = conditions,
		values = 1:length(conditions))
	p <- p + scale_y_log10()
	p <- p + scale_x_log10()
	p <- p + labs(
		x = 'Number of unique UMIs (thr)',
		y = 'Number of cutsites with #UMIs >= thr'
		, title = experiment
	)
	print(p)

	graphics.off()
}

plotField = function(
	chr_list, conditions, field, xlab, ylab,
	default = NULL,
	leg_title = '',
	col = c('#74C957', '#0068CB', '#CC2343', '#009688', '#FF9800'),
	dt = NULL, dt_field = NULL, dt_levels = NULL, dt_title = NULL,
	ncores = 1, verbose = 0
) {
	# Plot the specified column (or formula using column names)
	# of the umi_table.
	# 
	# Args:
	# 	chr_list (list): list of chromosomes
	# 	conditions (list): list of condition-specific data-frames
	# 	field (string): which field(s) or field-base operation to plot
	# 	xlab (string): x-axis label
	# 	ylab (string):y-axis label
	# 	default (*): default data series to plot
	# 	leg_title (string): legend title
	# 	col (list): list of colors (hexadec)
	# 	dt (data.frame): data-track for background
	# 	dt_field (string): data-track field
	# 	dt_levels (list): data-track levels
	# 	dt_title (string): data-track title
	# 	ncores (int): number of threads for parallelization
	# 	verbose (int): verbosity level
	# 	

	# Work on one chromosome at a time
	l <- lapply(chr_list,
		FUN=function(chr) {
			if ( verbose > 0 ) cat(paste0(' >>> Working on chr', chr, ' ...\n'))

			# Select conditions with data on the current chromosome
			sel_conds <- unlist(mclapply(conditions,
				FUN=function(condition) {
					if ( chr %in% names(umi_tab[[condition]]) )
						return(condition)
				}
				, mc.cores = ncores
			))
			sel_conds <- sel_conds[!is.null(sel_conds)]
			data <- lapply(sel_conds, FUN=function(x) umi_tab[[x]][[chr]])
			names(data) <- sel_conds

			# Add data-track background
			if ( !is.null(dt) & !is.null(dt_field) ) {

				# Extract ylims
				ylim <- do.call(rbind, mclapply(1:length(data),
					FUN=function(i) {
						tmp <- data[[i]]
						attach(tmp, warn.conflicts=F)
						v <- eval(parse(text = field))
						return(c(min(v, na.rm = T), max(v, na.rm = T)))
						detach(tmp)
					}
					, mc.cores = ncores
				))

				if ( !is.null(default) )
					ylim <- rbind(
						ylim,
						c(
							min(data[[1]][, default]),
							max(data[[1]][, default])
						)
					)
				ylim <- c(min(ylim[,1]), max(ylim[,2]))

				# Add data-track
				if ( !is.null(dt_levels) )
					lev <- dt_levels
				else
					lev <- unique(dt[, dt_field])

				# Select data_track fields & chromosome
				dt_sel <- data.frame(
					dt$start, dt$end,
					ylim[1], ylim[2],
					dt[, dt_field]
				)
				dt_sel <- dt_sel[dt$chr == chr,]
				colnames(dt_sel) <- c('start', 'end', 'ymin', 'ymax', dt_field)

				# Set colors for background
				gray_col <- gray.colors(length(lev))
				names(gray_col) <- rev(lev)

				# Initialize plot and add background
				p <- ggplot() + theme_classic() + geom_rect(data = dt_sel,
					aes(xmin = start, xmax = end,
						ymin = ymin, ymax = ymax,
						fill = eval(parse(text = dt_field))
				)) + scale_fill_manual(
					breaks = lev,
					labels = lev,
					values = gray_col,
					guide = guide_legend(title = dt_title)
				)

			} else {
				# Initialize plot
				p <- ggplot()
			}

			# Add tracks
			p <- p + mclapply(1:length(data),
				FUN=function(i) {
					data[[i]]$cond <- sel_conds[i]
					geom_line(
						data=data[[i]],
						aes(
							x = mid,
							y = eval(parse(text = field)),
							colour = cond
						),
						na.rm = T
					)
				}
				, mc.cores = ncores
			)

			# Add default track if present
			if ( !is.null(default) )
				p <- p + geom_line(
					data=data[[1]],
					aes(
						x = mid,
						y = eval(parse(text = default)),
						colour = 'default'
					),
					linetype='dotted', na.rm = T
				)

			# Select tracks colors
			sel_col <- c('#000000', col[which(conditions %in% sel_conds)])
			names(sel_col) <- c('default', sel_conds)
			p <- p + scale_colour_manual("",
				breaks = c('default', sel_conds),
				labels = c('default', sel_conds),
				values = sel_col
			)

			# Add labels
			p <- p + guides(colour = guide_legend(leg_title))
			p <- p + labs(
				x = xlab, y = ylab,
				title = paste0('[', experiment, '.', chr, '] ',
					'bin.size = ', bin_size, '; bin.step = ', bin_step)
			)

			# Output
			print(p)
			return(p)
		}
	)
}

plotDiffField = function(
	chr_list, conditions, field, xlab, ylab,
	default = NULL,
	leg_title = '', y_hline = NULL,
	col = c('#74C957', '#0068CB', '#CC2343', '#009688', '#FF9800'),
	dt = NULL, dt_field = NULL, dt_levels = NULL, dt_title = NULL,
	neg_id = 0, rm_neg = F, only_extr = F,
	rm.X = F, rm.Y = F,
	ncores = 1, verbose = 0
) {
	# Plot the difference between the specified conditions
	# of the value in the specified column (or formula using column names)
	# from the umi_table.
	# 
	# Args:
	# 	chr_list (list): list of chromosomes
	# 	conditions (list): list of condition-specific data-frames
	# 	field (string): which field(s) or field-base operation to plot
	# 	xlab (string): x-axis label
	# 	ylab (string):y-axis label
	# 	default (*): default data series to plot
	# 	leg_title (string): legend title
	# 	y_hline (numeric): y-coordinate for hline (baseline)
	# 	col (list): list of colors (hexadec)
	# 	dt (data.frame): data-track for background
	# 	dt_field (string): data-track field
	# 	dt_levels (list): data-track levels
	# 	dt_title (string): data-track title
	# 	neg_id (int): negative condition id
	# 	rm_neg (bool): remove negative condition, requires neg_id
	# 	only_extr (bool): plot only lower and higher condition
	# 	ncores (int): number of threads for parallelization
	# 	verbose (int): verbosity level
	# 
	
	# Remove negative
	if ( rm_neg & 0 != neg_id & neg_id <= length(conditions) ) {
		conditions = conditions[-neg_id]
		col = col[-neg_id]
	}

	# Select extreme conditions
	if ( only_extr ) {
		conditions = c(conditions[1], conditions[length(conditions)])
	}

	# Work on one chromosome at a time
	l <- lapply(chr_list,
		FUN=function(chr) {
			if ( verbose > 0 ) cat(paste0(' >>> Working on ', chr, ' ...\n'))

			# Select conditions with data on the current chromosome
			sel_conds <- unlist(mclapply(conditions,
				FUN=function(condition) {
					if ( chr %in% names(umi_tab[[condition]]) )
						return(condition)
				}
				, mc.cores = ncores
			))
			sel_conds <- sel_conds[!is.null(sel_conds)]

			# Check that at least two conditions are selected
			if ( 2 > length(sel_conds) ) return(NULL)

			# Prepare data table
			data <- lapply(sel_conds, FUN = function(x) umi_tab[[x]][[chr]])
			names(data) <- sel_conds

			# Couple of conditions to calculate difference between
			couples <- unlist(lapply(1:(length(sel_conds) - 1),
				FUN=function(i) {
					paste0(sel_conds[i + 1], '-', sel_conds[i])
				}
			))
			
			# Add data-track background
			if ( !is.null(dt) & !is.null(dt_field) ) {

				# Extract ylims
				ylim <- do.call(rbind, mclapply(2:length(data),
					FUN=function(i) {
						# Values for higher condition
						tmp <- data[[i]]
						attach(tmp, warn.conflicts=F)
						v <- eval(parse(text = field))
						detach(tmp)

						# Values for lower condition
						tmp <- data[[i - 1]]
						attach(tmp, warn.conflicts=F)
						v <- v - eval(parse(text = field))
						detach(tmp)

						return(c(min(v, na.rm = T), max(v, na.rm = T)))
					}
					, mc.cores = ncores
				))
				if ( !is.null(default) )
					ylim <- rbind(
						ylim,
						c(
							min(data[[1]][, default]),
							max(data[[1]][, default])
						)
					)
				ylim <- c(min(ylim[,1]), max(ylim[,2]))

				# Add data-track
				if ( !is.null(dt_levels) )
					lev <- dt_levels
				else
					lev <- unique(dt[, dt_field])

				# Select data_track fields & chromosome
				dt_sel <- data.frame(
					dt$start, dt$end,
					ylim[1], ylim[2],
					dt[, dt_field]
				)
				dt_sel <- dt_sel[dt$chr == chr,]
				colnames(dt_sel) <- c('start', 'end', 'ymin', 'ymax', dt_field)

				# Set colors for background
				gray_col <- gray.colors(length(lev))
				names(gray_col) <- rev(lev)

				# Initialize plot and add background
				p <- ggplot() + theme_classic() + geom_rect(data=dt_sel,
					aes(xmin = start, xmax = end,
						ymin = ymin, ymax = ymax,
						fill = eval(parse(text = dt_field))
				)) + scale_fill_manual(
					breaks = lev,
					labels = lev,
					values = gray_col,
					guide = guide_legend(title = dt_title)
				)

			} else {
				# Initialize plot
				p <- ggplot()
			}
			
			# Add tracks
			p <- p + mclapply(1:(length(data) - 1),
				FUN=function(i) {
					data[[i]]$cond <- couples[i]

					attach(data[[i+1]], warn.conflicts = F)
					data[[i]]$diff <- eval(parse(text=field))
					detach(data[[i+1]])

					attach(data[[i]], warn.conflicts = F)
					data[[i]]$diff <- data[[i]]$diff - eval(parse(text=field))
					detach(data[[i]])

					geom_line(
						data=data[[i]],
						aes(x = mid, y = diff, colour = cond),
						na.rm = T
					)
				}
				, mc.cores = ncores
			)

			# Add default track if present
			if ( !is.null(default) )
				p <- p + geom_line(
					data=data[[1]],
					aes(
						x = mid,
						y = eval(parse(text = default)),
						colour = 'default'
					),
					linetype='dotted', na.rm = T
				)

			# Select tracks colors
			sel_col <- c('#000000', col[which(conditions %in% sel_conds)])
			names(sel_col) <- c('default', couples)
			p <- p + scale_colour_manual("",
				breaks = c('default', couples),
				labels = c('default', couples),
				values = sel_col
			)

			# Add labels
			p <- p + guides(colour = guide_legend(leg_title))
			p <- p + labs(
				x = xlab, y = ylab,
				title = paste0('[', experiment, '.', chr, '] ',
					'bin.size = ', bin_size, '; bin.step = ', bin_step)
			)

			# Add baseline if requested
			if ( !is.null(y_hline) )
				p <- p + geom_hline(yintercept = y_hline, linetype = 'dotted')

			# Output
			print(p)
			return(p)
		}
	)
}

plotChrwideVar = function(
	dirpath, conditions, experiment, bin_size, cutsites,
	var_type = 'FF', rm.X = F, rm.Y = F,
	col = c('#74C957', '#0068CB', '#CC2343', '#009688', '#FF9800'),
	neg_id = 0,	ncores = 1, suff = ''
) {
	# Plot chr-wide FF
	# 
	# Args:
	# 	dirpath (string): base directory path (outdata)
	# 	conditions (list): list of condition-specific data-frames
	# 	experiment (string): experiment ID
	# 	bin_size (int): bin size in nt
	# 	cutsites (string): path to file with cutsite list
	# 	var_type (string): type of variance measure (FF|CV)
	# 	rm.X (bool): remove X chromosome
	# 	rm.Y (bool): remove Y chromosome
	# 	col (list): color list (hexadec)
	# 	neg_id (int): negative condition id
	# 	ncores (int): number of threads for parallelization
	# 	suff (string): suffix used for input/output operations
	# 
	
	# Check specified variance type
	if ( !var_type %in% c('CV', 'FF') ) return()
	
	# Set computation based on specified variance type
	if ( 'CV' == var_type ) {
		mvar = function(x) { sqrt(var(x)) / mean(x) }
	} else {
		mvar = function(x) { var(x) / mean(x) }
	}

	# Retrieve cutsite list
	el <- read.delim(cutsites, as.is=T, header=F)
	colnames(el) <- c('chr', 'pos')

	# Retrieve chr-wide Variance -----------------------------------------------
	
	t <- do.call(cbind, mclapply(conditions,
		FUN=function(condition) {
			fname <- 'UMIpos.unique'
			if ( file.exists(cutsites) ) fname <- paste0(fname, '.atcs')
			fname = paste0(dirpath, condition, '/', fname, suff, '.txt')
			u <- read.delim(fname, as.is = T, header = F)
			colnames(u) <- c('chr', 'pos', 'seq')

			# Remove X/Y chromosomes if requested
			toRM = c()
			if ( rm.X ) toRM = which(u$chr == 23)
			if ( rm.Y ) toRM = c(toRM, which(u$chr == 24))
			if ( 0 != length(toRM) ) u = u[-unique(toRM),]

			# Calculate variance values
			varvs <- by(u, u$chr,
				FUN=function(st, ncores) {
					chr <- st$chr[1]
					if (chr == 23) chr <- 'X'
					if (chr == 24) chr <- 'Y'
					
					# Count UMIs
					n_umi <- unlist(lapply(strsplit(as.character(st$seq), ' '),
						FUN = length))

					# Add non-sensed cutsites
					nrep = length(which(el$chr == paste0('chr', chr)))
					nrep = nrep - length(n_umi)
					n_umi <- c(n_umi, rep(0, nrep))

					# Return variance measure
					return(mvar(n_umi))

				}, ncores
			)
			names(varvs) <- paste0('chr', sort(unique(u$chr)))

			# Add missing chr
			chrlist <- paste0('chr', 1:24)
			if ( rm.X ) chrlist = chrlist[-23]
			if ( rm.Y ) chrlist = chrlist[-length(chrlist)]
			varvs[chrlist[! chrlist %in% names(varvs)]] <- 0
			varvs <- varvs[order(as.numeric(unlist(lapply(names(varvs),
				FUN = function(x) { substr(x, 4, nchar(x)) }))))]
			names(varvs)[names(varvs) == 'chr23'] <- 'chrX'
			names(varvs)[names(varvs) == 'chr24'] <- 'chrY'

			return(varvs)
		}
		, mc.cores=ncores
	))
	colnames(t) <- conditions

	# Plot Var per condition barplot -------------------------------------------
	
	p <- ggplot(data = melt(t), aes(x = Var1, y = value, fill = Var2))
	p <- p + geom_bar(stat = 'identity', position = position_dodge())
	p <- p + guides(fill = guide_legend('Conditions'))
	p <- p + scale_fill_manual('',
		breaks = conditions,
		labels = conditions,
		values = col
	)
	p <- p + labs(
		x = 'Chromosome',
		y = paste0(var_type, ' of #UMI per cutsite'),
		title = paste0('[', experiment, '] Chromosome-wide ', var_type)
	)
	print(p)

	# Remove negative
	if ( 0 < neg_id & neg_id <= ncol(t) ) {
		t = t[,-neg_id]
		conditions = conditions[-neg_id]
		col = col[-neg_id]
	}

	# Plot deltaVar ------------------------------------------------------------
	
	t2 <- data.frame(
		chr = rownames(t),
		diff = t[,ncol(t)] - t[,1],
		stringsAsFactors = F)
	levels(t2$chr) <- t2$chr
	t2 <- t2[order(t2$diff),]
	t2$chr <- factor(t2$chr, levels=t2$chr)
	label <- paste0(conditions[ncol(t)], '-', conditions[2])

	p <- ggplot(data = t2, aes(x = chr, y = diff))
	p <- p + geom_bar(stat = 'identity', position = 'identity',
		fill = col[length(conditions)])
	p <- p + scale_fill_manual('',
		breaks = col[length(conditions)],
		labels = label,
		values = col[length(conditions)]
	)
	title = paste0('[', experiment, '] Chromosome-wide delta',
		var_type, ' (', label, ')')
	p <- p + labs(
		x = 'Chromosome',
		y = paste0('delta', var_type, ' of #UMI per cutsite'),
		title = title
	)
	p
}

plotGenome = function(
	dirpath, experiment, bin_size, bin_step,
	maskfile = NULL,
	rm.X = F, rm.Y = F,
	global = F
) {
	# Generate heatmap with fraction of binned chromosome at specified condition
	# 
	# Args:
	# 	condition (string): condition ID
	# 	dirpath (string): path to analysis main experiment
	# 	experiment (string): experiment ID
	# 	bin_size (int): bin size in nt
	# 	rm.X (bool): remove the X chromosome
	# 	rm.Y (bool): remove the Y chromosome
	# 	global (bool): perform condition- instead of chr-normalization
	# 
	
	# Input --------------------------------------------------------------------
	
	# Read umi_tab
	load(paste0(dirpath, '/', experiment, '.umi_table.bsize',
		bin_size, '.bstep', bin_step, '.RData'))

	# Read maskfile
	mf = maskfile

	# Prepare condition-specific UMI table -------------------------------------
	umis = lapply(names(umi_tab),
		FUN = function(condition) {

			# Merge chromosome-specific tables
			t = rbindlist(lapply(names(umi_tab[[condition]]),
				FUN = function(chr) {

					# Select chromosome table columns
					sel_cols = c('start', 'end', 'n.sum')
					chr_tab = umi_tab[[condition]][[chr]][,sel_cols]

					# Prepare chromosome number
					chr_id = substr(chr, 4, nchar(chr))
					chr_id = ifelse(chr_id == 'X', 23, chr_id)
					chr_id = ifelse(chr_id == 'Y', 24, chr_id)
					
					# Select chromosome masked regions
					mfchr = mf[mf$num == chr_id,]

					# Add chromosome info to chromosome table
					chr_tab = data.frame(
						chr = rep(chr, nrow(chr_tab)),
						chr_num = as.numeric(rep(chr_id, nrow(chr_tab))),
						start = chr_tab$start,
						end = chr_tab$end,
						n = chr_tab$n.sum,
						p = chr_tab$n.sum / sum(chr_tab$n.sum),
						stringsAsFactors = T
					)

					# Select non-overlapping bins
					breaks = seq(min(chr_tab$start, na.rm = T),
						max(chr_tab$end, na.rm = T), bin_size)
					breaks = breaks[-length(breaks)]
					sel_bins = paste0(breaks, '-', breaks + bin_size)
					chr_tab$sign = paste0(chr_tab$start, '-', chr_tab$end)
					chr_tab = chr_tab[chr_tab$sign %in% sel_bins,]
					chr_tab$sign = paste0(chr, ':', chr_tab$sign)

					# Identify masked bins
					toMask = unlist(lapply(1:nrow(chr_tab),
						FUN = function(i) {
							bin = c(chr_tab[i,]$start, chr_tab[i,]$end)
							cond = bin[1] >= mfchr$start
							cond = cond & bin[2] < mfchr$end
							cond = any(cond)
							return(cond)
						}
					))
					chr_tab[toMask, 'n'] = NA
					chr_tab[toMask, 'p'] = NA

					# Output
					return(chr_tab)
				}
			), use.names = T)

			# Remove chromosome X, if present
			if ( rm.X ) {
				toRM = which(t$chr_num == 23)
				if ( 0 != length(toRM) )
					t = t[-toRM,]
			}

			# Remove chromosome Y, if present
			if ( rm.Y ) {
				toRM = which(t$chr_num == 24)
				if ( 0 != length(toRM) )
					t = t[-toRM,]
			}

			# Normalize on the overall genome UMI count
			t$gp = t$n / sum(t$n, na.rm = T)

			# Output
			return(t)
		}
	)
	names(umis) = names(umi_tab)

	# Plot ---------------------------------------------------------------------
	ps = lapply(names(umis),
		FUN = function(condition) {
			title = paste0(experiment, '~', condition, ' : ')

			if ( global ) {
				title = paste0(title, 'normalized on the whole condition')
				p = ggplot(umis[[condition]],
					aes(x = chr, y = start, fill = gp))
			} else {
				title = paste0(title, 'normalized per chromosome')
				p = ggplot(umis[[condition]],
					aes(x = chr, y = start, fill = p))
			}
			p = p + geom_tile()
			p = p + scale_fill_gradient(low = "white", high = "red",
				na.value = 'black')
			p = p + labs(title = title,
				xlab = 'Chromosome', ylab = 'Coordinates')
			print(p)
		}
	)
}

# END --------------------------------------------------------------------------

################################################################################
