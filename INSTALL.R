#!/usr/bin/env Rscript

# RUN ==========================================================================

# Required packages
required = c(
	"argparser",
	"data.table",
	"ggplot2",
	"grid",
	"gridExtra",
	"parallel",
	"readr",
	"reshape2"
)

# Checking if needed
to_install = c()
for ( package in required ) {
	test = require(package = package, character.only = T, quietly = T)
	if ( !test ) {
		to_install = c(to_install, package)
	}
}

# Install
if ( !is.null(to_install) ) {
	print(paste0("Installing the following packages: ",
		paste(to_install, collapse = ", ")))
	install.packages(to_install)
}

# END ==========================================================================

################################################################################
