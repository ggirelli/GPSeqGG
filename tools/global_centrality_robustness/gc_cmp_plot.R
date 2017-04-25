#!/usr/bin/Rscript

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Description: plot shuffled global centrality comparison.
# 
# ------------------------------------------------------------------------------



# DEPENDENCIES =================================================================

suppressMessages(library(argparser))
suppressMessages(library(data.table))

# INPUT ========================================================================

# Create arguent parser
parser = arg_parser('Plot shuffled global centrality comparison.',
	name = 'gc_cmp_plot.R')

# Define mandatory arguments
parser = add_argument(parser, arg = 'percs', short = "-p", nargs = 1,
	help = 'Percentage of reshuffled reads, csv', type = class(''),
	default = '')
parser = add_argument(parser, arg = 'outFile', short = "-o", nargs = 1,
	help = 'Output file, pdf.', type = class(''),
	default = 'gc.cmp.pdf')

# Define optional arguments
parser = add_argument(parser, arg = '--cmpFiles', short = "-c", nargs = Inf,
	help = 'Kendall-tau distance matrices')
parser = add_argument(parser, arg = '--label', short = "-l", nargs = 1,
	help = 'Dataset label.', type = class(''),
	default = 'data')

# Parse arguments
p = parse_args(parser)

# Attach argument values to variables
attach(p['' != names(p)])

# Split percentages
percs = as.numeric(unlist(strsplit(percs, ',', fixed = T)))

# Check file list
if ( !all(unlist(lapply(cmpFiles, FUN = file.exists))) ) {
	stop("!!!ERROR!!! Cannot find specified files.\n")
}

# Color palette
colors = c('#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99', '#e31a1c',
	'#fdbf6f', '#ff7f00', '#cab2d6', '#6a3d9a', '#FF2299', '#b15928', '#000000',
	'#656565')

# RUN ==========================================================================

# Prepare table with metric value, metric label and shuffling percentage cols
t = rbindlist(lapply(seq(length(cmpFiles)),
	FUN = function(i) {
		t = read.delim(paste0(cmpFiles[i]), as.is = T)

		t2 = rbindlist(lapply(colnames(t),
			FUN = function(cn) {
				data.frame(
					value = t[,cn],
					label = rep(cn, nrow(t)),
					perc = percs[i]
				)
			}
		))

		return(as.data.frame(t2))
	}
))

# Identify number of iterations
nIter = unique(as.numeric(by(t, t$label,
	FUN = function(x) { unique(table(x$perc)) })))

# Identify max medians
max_median = max(unlist(by(t, t$label,
	FUN = function(x) { by(x, x$perc,
		FUN = function(y) { median(y$value) }) })))

# Prepare plot output
pdf(outFile, width = 20, height = 10)
par(mfrow = c(1,2), oma = c(0,0,2,0), mar = c(5,5,2.5,2))

# MEDIANS
# Plot per metric type
labels = as.character(unique(t$label))
l = lapply(seq(length(labels)),
	FUN = function(li) {
		# Select current metric type
		label = labels[li]
		subt = as.data.frame(t)[as.character(t$label) == label,]

		# Calculate medians
		medians = as.numeric(by(subt, subt$perc,
			FUN = function(x) {
				median(x$value)
			}
		))

		# Set plot window with first metric, then add to the plot
		if ( 1 == li ){
			plot(medians ~ percs, type = 'b', col = colors[li],
				ylim=c(0, max_median), xlab='Percentage of reshuffled reads',
				ylab='Median Kendall tau distance from original ranking')
		} else {
			lines(medians ~ percs, type = 'b', col = colors[li])
		}
	}
)
legend(10, max_median, col=colors, labels, lty=1, pch=1)

# MAXIMAS
# Plot per metric type
l = lapply(seq(length(labels)),
	FUN = function(li) {
		# Select current metric type
		label = labels[li]
		subt = as.data.frame(t)[as.character(t$label) == label,]

		# Calculate maximas
		maxs = as.numeric(by(subt, subt$perc,
			FUN = function(x) {
				max(x$value)
			}
		))

		# Set plot window with first metric, then add to the plot
		if ( 1 == li ){
			plot(maxs ~ percs, type = 'b', col = colors[li],
				ylim=c(0, max(t$value)), xlab='Percentage of reshuffled reads',
				ylab='Max Kendall tau distance from original ranking')
		} else {
			lines(maxs ~ percs, type = 'b', col = colors[li])
		}
	}
)
legend(10, max(t$value), col=colors, labels, lty=1, pch=1)

# Add more text
title("Study of global centrality metrics robustness",
	outer=TRUE)
mtext(paste0("(", label, "; full-chromosome windows; nIter = ", nIter, ")"),
	side = 3, line = -1, outer = T)
mtext("Shuffling performed on the whole genome at the single-reads level.",
	side = 3, line = -2, outer = T)

# Write output
graphics.off()

# END --------------------------------------------------------------------------

################################################################################
