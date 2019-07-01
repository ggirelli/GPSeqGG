# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).



## Unreleased
### ...
- ...



## [2.0.3] - 2019-07-01
### Changed
- Optimized UMI de-duplication step.

### Fixed
- Fraction of final output VS input with prefix.


## [2.0.2] - 2018-04-04
### Changed
- Updated submodules.

### Fixed
- Now properly summarizing preseq reports with multiqc.
- Fixed alignment fraction calculation.



## [2.0.1] - 2018-01-16
### Fixed
- Fixed secondary alignments removal during SAM filtering step.
- Fixed chromosome removal during SAM filtering step.

### Added
- `multiqc` report recap as a last step of the pipeline.
- Version tag to log.



## [2.0.0] - 2017-01-09
### Added
- `preseq` package as a git submodule.
- Library complexity estimation step, using `preseq` package.
- Version tracking to bed header and summary tables.

### Changed
- Clarified neatness level labeling.
- Promoted centrality estimate from `dev` to `tool`.
- Now using `samtools` and `awk` for SAM file filtering.
- Updated documentation.

### Removed
- Useless thread declaration when running FASTQC tool.
- Binning and plotting step, stored in `deprec` folder.

### Fixed
- Input option exceptions.



## [1.1.0] - 2017-07-26
### Added
- Neatness level for automatic clean-up of intermediate files.
- Multi-library analysis.
- Saving orphan reads in orphan.txt.

### Changed
- Cutsite sequence option, moved to single option from pattern file.
- pat_files to patterns.tsv (renamed).
- Log, now divided by library.

### Fixed
- Bug in SAM filter that crashed the script if no chromosome is removed.

### Removed
- Sub-step queries for filepaths.



## [1.0.1] - 2017-07-17
### Added
- Debug mode to global centrality script.
- Optional intermediate ranking to global centrality script.
- Counting UMI duplicates during strict deduplication step.

### Changed
- Read quality filter, optimized memory usage.
- Summary update method to avoid creationg of broken lines.

### Fixed
- Global centrality option.
- Bed file creation: multiple pipeline run would append new bed files to old ones.



## [1.0.0] - 2017-03-29.
### Added
- Input parameter support with getopt in main.
- Auxillary intermediate data `aux` subfolder.
- Counts of reads filtered from SAM in the summary.
- Command line saved in aux/CMD file.
- X/Y chromosomes filter during SAM filter step.
- MAPQ level in SAM filter output.
- SAM filter step log.
- UMI preparation step log.
- Header line to bed file.
- Global centrality script.
- Cumulative of ratio and ratio of cumulative approaches to global centrality script.
- Properly formatted bed files per condition.

### Changed
- File location, moved code to `bin` subfolder.
- Output structure, moved conditions into `conditions` subfolder.
- Output structure, moved output plots into `plots` subfolder.
- Read trimming script, optimized.
- Main script, now using modules.
- File structure, moved scripts to separate folder.
- Pos2cutsite script, optimized by calculating distance only from the 2 closest cutsite.
- Bed creation, moved to separate script.
- Barcode sequence option, moved to pat_files.

### Fixed
- Read trimming.
- rmX/rmY option report.



## [0.1.0] - 2017-03-13.



[Unreleased] https://github.com/ggirelli/gpseq-seq-gg
[2.0.3] https://github.com/ggirelli/gpseq-seq-gg/releases/tag/v2.0.3
[2.0.2] https://github.com/ggirelli/gpseq-seq-gg/releases/tag/v2.0.2
[2.0.1] https://github.com/ggirelli/gpseq-seq-gg/releases/tag/v2.0.1
[2.0.0] https://github.com/ggirelli/gpseq-seq-gg/releases/tag/v2.0.0
[1.1.0] https://github.com/ggirelli/gpseq-seq-gg/releases/tag/v1.1.0
[1.0.1] https://github.com/ggirelli/gpseq-seq-gg/releases/tag/v1.0.1
[1.0.0] https://github.com/ggirelli/gpseq-seq-gg/releases/tag/v1.0.0
[0.1.0] https://github.com/ggirelli/gpseq-seq-gg/releases/tag/v0.1.0
