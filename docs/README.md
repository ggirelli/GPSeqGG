Genomic loci Positioning by Sequencing (GPSeq) sequencing data analysis (v1.0.0)
===

## The `main` script

```
usage: ./main.sh [-h][-w][-t threads] -i inDir -o outDir -e expID
 [-n][-a aligner][-g refGenome][-d bwaIndex][-x][-y][-q mapqThr]
 [-p platform][-u umiLength][-r csRange][-j emax][-k eperc][-z binSize]
 [-b binStep][-l csList][-m maskFile][-s chrLengths]

 Description:
	Run a step-by-step interactive GPSeq sequencing data analysis.

 Required files:
	Requires R1 (and R2 if paired-end sequencing) and a pattern files in the
	input directory (inDir). The patterns file should have a condition per row
	with condition name, pattern (scan_for_mateches format), cutsite sequence and
	non-genomic portion length, separated by tabulations.

 Mandatory arguments:
	-i indir      Input directory.
	-o outdir     Output directory. Created if not found.
	-e expID      Experiment ID.

 Optional arguments:
	-h            Show this help page.
	-w            Perform every step of the pipeline without asking.
	-n            Negative present. Expected label: neg;
	-x            Remove X chromosome after alignment.
	-y            Remove Y chromosome after alignment.
	-t threads    Number of threads for parallelization.
	-a aligner    Aligner. Either 'bwa' (default) or 'bowtie2'.
	-g refGenome  Path to reference genome file. Default: 'hg19'.
	-d bwaIndex   Path to BWA index file. Required if BWA is the selected aligner.
	-q mapqThr    Mapping quality threshold. Default: 1.
	-p platform   Sequencing platform. Default: 'L'.
	-u umilength  UMI sequence length. Default: 8.
	-r csRange    Range around cutsite for UMI assignment. Default: 40.
	-j emax       Maximum error probability for read quality filtering. Default: 1e-3.
	-k eperc      Maximum % of bases with emax error probability. Default: 20.
	-z binSize    Bin size. Default: 1e6.
	-b binStep    Bin step. Default: 1e5.
	-l csList     File with cutsite list. Columns: chr|pos. No header.
	-m maskFile   File with masked regions. Columns: id|chr|start|end. No header.
	-s chrLengths File with chromosome lengths. chr|len. No header.
```

#### Pipeline steps

The steps performed by `main` are the following:

1. Run `./quality_control.sh`.
2. Run `./files_prepare.sh`
3. Perform pattern filtering:
		* Count total reads and prepare `summary` file
		* Run `./pattern_filtering.sh` per condition and update the `summary` file
4. Perform read alignment.
		* First, run `./reads_trim.sh` trim the pattern from the R1 reads (length specified in `pat_files`)
		* Run `./reads_align.sh` and update summary with fraction of mapped reads
		* Add the patterns to the SAM files (condition.linker.sam)
5. Filters the SAM file. (this step requires a specified `mapqthr` setting)
		* Run `./sam_filter`.
6. Prepare UMIs. If a `maskfile` was specified, UMIs are masked based on the information contained in such file.
		* Run `./umi_group` to group UMIs on the same position.
		* Run `./pos2cutsite` to group UMIs on a cutsite
		* Run `./umi_dedupl` to deduplicate UMIs.
7. Performs binning.
		* Run `./cs_bin` to bin the cutsites.
		* Run `./umi_bin` to bin the UMIs.
8. Analyzes UMIs.
		* Run `./umi_plot` to generate analysis results as figures.

#### Output

Description of different output files.

* [summary](summary/)

#### Input

At least two files are required: a sequencing platform output file (R1), and a `pat_files` file in the input directory, which contains the pattern to recognize the different conditions. If the sequencing is paired ended, also an R2 file is required.

The `patterns` file, which is expected to be located in the input folder, contains a row per condition per experiment. Every row contains the following tab-separated information:

```
experiment_ID	condition_label	scan_for_matches_pattern	cutsite_sequence	trimming_length
```

An example being the following:

```
TK20	neg	^ 8...8 CATCAGAA AAGCTT 1...1000 $	AAGCTT	22
TK20	1min	^ 8...8 CATCATCC AAGCTT 1...1000 $	AAGCTT	22
TK21	1min	^ 8...8 GTCGTTCC AAGCTT 1...1000 $	AAGCTT	22
TK21	2h	^ 8...8 TGATGTCC AAGCTT 1...1000 $	AAGCTT	22
```

## Other scripts

...