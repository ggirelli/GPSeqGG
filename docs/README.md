Genomic loci Positioning by Sequencing (GPSeq) sequencing data analysis
===

## The `main` script

```
usage: ./main.sh [-h][-w][-t threads] -i inDir -o outDir -e expID -c conditions
 [-n][-a aligner][-g refGenome][-d bwaIndex][-x][-y][-f cutsite][-q mapqThr]
 [-p platform][-u umiLength][-r csRange][-j emax][-k eperc][-z binSize]
 [-b binStep][-l csList][-m maskFile][-s chrLengths]

 Description:
  Run a step-by-step interactive GPSeq sequencing data analysis.

 Required files:
  Requires R1 (and R2 if paired-end sequencing) and a pattern files in the
  input directory (inDir). The patterns file should have a condition per row
  with condition name, pattern (scan_for_mateches format) and pattern length
  separated by tabulations.

 Mandatory arguments:
  -i indir      Input directory.
  -o outdir     Output directory. Created if not found.
  -e expID      Experiment ID.
  -c conditions Comma-separated conditions.

 Optional arguments:
  -h            Show this help page.
  -y            Perform every step of the pipeline without asking.
  -n            Negative present. Expected label: neg;
  -x            Remove X chromosome after alignment.
  -y            Remove Y chromosome after alignment.
  -t threads    Number of threads for parallelization.
  -a aligner    Aligner. Either 'bwa' (default) or 'bowtie2'.
  -g refGenome  Path to reference genome file. Default: 'hg19'.
  -d bwaIndex   Path to BWA index file. Required if BWA is the selected aligner.
  -f cutsite    Cutsite sequence. Default: 'AAGCTT' (HindIII).
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

#### Input

At least two files are required: a sequencing platform output file (R1), and a `pat_files` file in the input directory, which contains the pattern to recognize the different conditions. If the sequencing is paired ended, also an R2 file is required.

The `pat_files`, which is expected to be located in the input folder, contains a row per condition (as specified in the `-c conditions` argument). Every row contains the following tab-separated information: condition label, scan_for_matches pattern, pattern length. An example being the following:

    neg ^ 8...8 CATCAGAA AAGCTT 1...1000 $  22
    1min    ^ 8...8 CATCATCC AAGCTT 1...1000 $  22
    10min   ^ 8...8 GTCGTTCC AAGCTT 1...1000 $  22
    2h  ^ 8...8 TGATGTCC AAGCTT 1...1000 $  22

## Other scripts

### ./quality_control.sh

This script is run by `main` during step 1. It performs a quality control with the `fastqc` tool.

```
 usage: ./quality_control.sh [-h][-t threads] -o outDir -1 r1file [-2 r2file]

 Description:
  Run FASTQC tool.

 Mandatory arguments:
  -o outDir Output folder, created if not found.
  -1 r1file R1 fastq file.

 Optional arguments:
  -h        Show this help page.
  -t threads    Number of threads for parallelization
  -2 r2file R2 fastq file.
```

### ./files_prepare.sh

This script is run by `main` during step 2. It decompresses the sequencing data and converts the generated files in different 'formats'.

```
 usage: ./files_prepare.sh [-h][-t threads] -o outDir -1 r1file [-2 r2file]

 Description:
  Run FASTQC tool.

 Mandatory arguments:
  -o outDir Output folder, created if not found.
  -1 r1file R1 fastq file.

 Optional arguments:
  -h        Show this help page.
  -t threads    Number of threads for parallelization
  -2 r2file R2 fastq file.
```

The produced files are the following:

* `r1.fq` (and `r2.fq`): decompressed **fastq** file(s).
* `r2oneline.fq` (and `r2oneline.fq`): one read per line **fastq** file(s).
* `r1.fa` (and `r2.fa`): decompressed **fasta** file(s).
* `r1oneline.fa` (and `r2oneline.fa`): one read per line **fasta** file(s).

### ./pattern_filter

This script is run by `main` during step 3. It runs the `scan_for_matches` tool to select reads based on a pattern (e.g., containing barcodes), and divide them into different conditions.

```
 usage: ./pattern_filter.sh [-h][-t threads] -i inDir -o outDir -p patFile

 Description:
  Select reads based on the provided pattern by runnin scan_for_matches.

 Mandatory arguments:
  -i inDir  Input directory.
  -o outDir Output directory.
  -p patFile    Pattern file.

 Optional arguments:
  -h        Show this help page.
  -t threads    Number of threads for parallelization.
```

The script will create a `filtered.r1.fa` and a `filtered.r1.fq` (and respective r2 files in case of paired-end sequencing) in the specified `outDir`. These files will contain the reads that match the specified pattern, retrieved from `r1.fa` and `r1.fq` (and r2...).

Finally, the `summary` file is updated with the number and percentage of reads belonging to that condition.

### ./reads_trim.sh

This BASH script is run by `main` in step 4. It removes and stores the linkers from the R1 reads to a separate file.

```
 usage: ./reads_trim.sh [-h][-t threads] -o outDir -c cond -p patFile

 Description:
  Trim linker from reads.

 Mandatory arguments:
  -o outDir     Output directory. Created if not found.
  -c cond       Condition to analyze.
  -p patFile    Pattern file.

 Optional arguments:
  -h            Show this help page.
  -t threads    Number of threads for parallelization.
```

### ./reads_align.sh

This BASH script is run by `main` in step 4. It aligns the reads to a reference genome using either `bwa` or `bowtie2` and generates `bam` and `bai` files.

```
 usage: ./reads_align.sh [-h][-t threads] -o outDir -c cond
        [-p][-r ref][-a aln]

 Description:
  Align reads using either bowtie2 or bwa mem.

 Mandatory arguments:
  -o outDir     Output directory. Created if not found.
  -c cond       Condition to analyze.

 Optional arguments:
  -h            Show this help page.
  -t threads    Number of threads for parallelization.
  -p            Option for paired-end sequencing.
  -r ref        Reference ref. Default: hg19.
  -a aln        Aligner. Either 'bwa' or 'bowtie2'. Default: 'bwa'.
  -i index      Path to BWA index.
                Default: '$DATA'/BiCro-Resources/genomes/'$ref'bwt/'$ref'.fa
```

### ./sam_filter.R

This R script is run by `main` in step 5. It filters a SAM file based on a variety of settings.

```
usage:  sam_filter.R [--help][--cutsite CUTSITE][--mapq-thr MAPQ-THR]
        [--num-proc NUM-PROC][--suffix SUFFIX] dirpath experiment condition

Filter alignment output.

positional arguments:
  dirpath           Experiment ccondition directory, contains UMI counts.
  experiment        Experiment ID (e.g., TK26).
  condition         Condition folder name (e.g., 400U2h).

flags:
  -h, --help        show this help message and exit

optional arguments:
  -cs, --cutsite CUTSITE        Cutsite sequence (e.g., for HindIII AAGCTT).
                                [default: AAGCTT]
  -mt, --mapq-thr MAPQ-THR      Mapping quality threshold.
                                [default: 30]
  -c, --num-proc NUM-PROC       Number of cores for parallel computation.
                                [default: 1]
  -s, --suffix SUFFIX           Suffix to be added to output files.
                                [default: '']
```

Specifically, the filters do the following:

* Remove secondary alignments (based on SAM file flags).
* If paired-end sequencing, remove chimeras (where R1 and R2 are mapped on different chromosomes).
* Remove unmapped reads.
* If paired-end sequencing, remove R2 and keep only R1.
* Remove reads with a mapping quality lower than `mapqthr`.

Then, the generated output consists of two files:

* `condition.filtered.sam.RData` contains the raw filtered SAM file as an RData structure.
* `condition.filtered.umi.pos.txt` is a table with the following columns: chromosome, 3'position, linker sequence, linker quality, read name. The position does not correspond to the SAM file one, but is instead shifted accordingly to the CIGAR string in the SAM file or to the length of the cutsite.

### ./umi_group.py

This Python script is run by `main` in step 6. It groups reads mapped on the same genomic locus.

```
usage: umi_group.py [-h] [--mask-file MF] dirpath expID cond umiLen

Group UMIs based on their location.

positional arguments:
  dirpath         Experiment condition directory, contains UMI counts.
  expID           Experiment ID (e.g., TK26).
  cond            Condition folder name (e.g., 400U2h).
  umiLen          Length of the UMIs in nt (e.g., 8).

optional arguments:
  -h, --help      show this help message and exit
  --mask-file MF  File with regions to be masked.
```

The output consists of:

* `UMIpos.txt`, containing the following columns:
    - chromosome
    - position
    - space-separated UMIs
    - quality strings
* If a `maskfile` was specified, then an additional `UMIpos.all.txt` file will be generated, containing all the UMIs (also those withing masked regions).

### ./pos2cutsite.R

This R script is run by `main` in step 6. It groups reads mapped on the same cutsite. Every read is assigned to the closest cutsite. Then, a read is discarded and considered 'orphan', if the closest cutsite is more than bin_size/2 nt away from the aligned locus. In other words, a window of bin_size nt is built around every cutsite.

The `pos2cutsite` step is not mandatory and, if the list of custite in the current reference genome is not known, the analysis can be performed anyway.

```
usage:  pos2cutsite.R [--help][--bin-size BIN-SIZE][--num-proc NUM-PROC]
        [--suffix SUFFIX] dirpath experiment condition cutsites

Convert genomic coordinate to cutsite.

positional arguments:
  dirpath           Experiment ccondition directory, contains UMI counts.
  experiment        Experiment ID (e.g., TK26).
  condition         Condition folder name (e.g., 400U2h).
  cutsites          File containing the cutsite positions. (chr | pos)

flags:
  -h, --help        show this help message and exit

optional arguments:
  -i, --bin-size BIN-SIZE   Range size in bp around the cutsite.
                            [default: 40]
  -c, --num-proc NUM-PROC   Number of cores for parallel computation.
                            [default: 1]
  -s, --suffix SUFFIX       Suffix to be added to output files.
                            [default: '']
```

### ./umi_dedupl.R

This R script is run by `main` in step 6. It removes duplicated (through amplification) reads based on the UMIs.

```
usage:  umi_dedupl.R [--help][--cutsites CUTSITES][--platform PLATFORM]
        [--cutoff CUTOFF][--emax EMAX][--eperc EPERC][--num-proc NUM-PROC]
        [--suffix SUFFIX] dirpath experiment condition

Deduplicate UMIs.

positional arguments:
  dirpath           Experiment condition directory, contains UMI counts.
  experiment        Experiment ID (e.g., TK26).
  condition         Condition folder name (e.g., 400U2h).

flags:
  -h, --help        show this help message and exit

optional arguments:
  -cs, --cutsites CUTSITES  Binary flag for cutsite assignment.
                            1 for, and 0 for no, cutsites
                            [default: 0]
  -p, --platform PLATFORM   Sequencing platform identifier.
                            [default: L]
  -co, --cutoff CUTOFF      Probability cutoff,
                            compared to the automatic for filtering.
                            [default: 1]
  -em, --emax EMAX          Maximum error probability for filtering.
                            [default: 0.001]
  -ep, --eperc EPERC        Maximum percentage of bases
                            with emax error probability. [default: 20]
  -c, --num-proc NUM-PROC   Number of cores for parallel computation.
                            [default: 1]
  -s, --suffix SUFFIX       Suffix to be added to output files.
                            [default: '']
```

First, the length of the UMIs is checked for inconsistencies.

Then, if a manual automatic `cutoff` higher than 0 was specified, a self-match probability filter is applied. Also, a single-base quality filter is applied that discards any read containing more than `eperc`% bases with error probability higher than or equal to `emax`.

Finally, a strict unique operation is performed that removes duplciated UMIs merely based on their sequence.

A figure report is generated (`umi_dedup.report.png`), showing the intensity of the detected UMI duplication.

The output consists of a single file, either `UMIpos.unique.txt` or `UMIpos.unique.atcs.txt`. `atcs` stands for **AT CutSite**. In other words, the `pos2cutsite` step is not mandatory and, if the list of custite in the current reference genome is not known, the analysis can be performed anyway.

The output text files contains the following columns:

* chromosome
* position
* space-separated UMI sequences

### ./cs_bin.R

This R script is run by `main` in step 7, and performs the cutsite binning. The output is a `cs_per_bin.bsize.bstep.RData` file containing the binned cutsite.

```
usage:  cs_bin.R [--help][-i BIN-SIZE][-t BIN-STEP]
        [-c NUM-PROC][-s SUFFIX] dirpath cutsites chrlengths

Bin cutsites.

positional arguments:
  dirpath                     Directory to which the binned cutsites will be saved.
  cutsites                    File containing the cutsite positions. (chr | pos)
  chrlengths                  File containing the chromosome lengths. (chr | len)

flags:
  -h, --help                  show this help message and exit

optional arguments:
  -i, --bin-size BIN-SIZE     Bin size in bp.
                              [default: 1e+06]
  -t, --bin-step BIN-STEP     Distance between the starting point of consecutive bins in bp.
                              [default: 1e+05]
  -c, --num-proc NUM-PROC     Number of cores for parallel computation.
                              [default: 1]
  -s, --suffix SUFFIX         Suffix to be added to output files.
                              [default: ]
```

### ./umi_bin.R

This R script is run by `main` in step 7, and performs the UMI binning. The output is a `experiment.umi_table.bsize.bstep.RData` file containing the binned cutsite. Also, a `experiment.condition.umi_table.bsize.bstep.tsv` file is saved in every condition folder.

```
usage:  umi_bin.R [-h][-c][-i BIN-SIZE][-t BIN-STEP][-p NUM-PROC]
        [-s SUFFIX] dirpath experiment conditions chrlengths

Bin UMIs.

positional arguments:
  dirpath         Experiment main directory, contains binned cutsites if available.
  experiment      Experiment ID (e.g., TK26).
  conditions      Comma-separated conditions (e.g., 400U2h,200U1h).
  chrlengths      File containing the chromosome lengths. (chr | len)

flags:
  -h, --help           show this help message and exit
  -c, --cutsites       Whether cutsites where used.

optional arguments:
  -i, --bin-size BIN-SIZE     Bin size in bp.
                              [default: 1e+06]
  -t, --bin-step BIN-STEP     Distance between the starting point of consecutive bins in bp.
                              [default: 1e+05]
  -p, --num-proc NUM-PROC     Number of cores for parallel computation.
                              [default: 1]
  -s, --suffix SUFFIX         Suffix to be added to output files.
                              [default: ]

```

### ./umi_plot.R

This R script is executed by `main` in step 8. It generates the final plots of the UMI analysis. The actual plot code is stored in `umi_plot.functions`, which is imported and required by the script.

```
usage:  umi_plot.R [--] [--help] [--rmChrX] [--rmChrY] [--opts OPTS]
    [--neg NEG] [--bin-size BIN-SIZE] [--bin-step BIN-STEP]
    [--num-proc NUM-PROC] [--suffix SUFFIX]
    dirpath experiment conditions cutsites chrlengths mask

Generate UMI plots.

positional arguments:
dirpath       Experiment condition directory, contains UMI counts.
experiment    Experiment ID (e.g., TK26).
conditions    Comma-separated conditions (e.g., 400U2h,200U1h).
cutsites      File containing the cutsite positions. (chr|pos)
chrlengths    File containing the chromosome lengths. (chr|len)
mask          File containing the masked regions.
            (num|chr|start|end|type)

flags:
-h, --help        show this help message and exit
-r, --rmChrX      Remove ChrX from the analysis
--rmChrY          Remove ChrY from the analysis

optional arguments:
-x, --opts OPTS           RDS file containing argument values
-n, --neg NEG             Negative condition label.
                        [default: '']
-i, --bin-size BIN-SIZE   Bin size in bp.
                        [default: 1e+06]
-t, --bin-step BIN-STEP   Distance between the starting point
                        of consecutive bins in bp.
                        [default: 1e+05]
-c, --num-proc NUM-PROC   Number of cores for parallel computation.
                        [default: 1]
-s, --suffix SUFFIX       Suffix to be added to output files.
                        [default: '']
```

Generated plots:

* UMI distribution across condition.
    - Saved as `experiment.umi_distribution.suff.pdf`.
* Sensed cutsites profile.
    - Saved as `experiment.sensed_cutsites.bsize.bstep.pdf`.
* UMI count profiles.
    - Plots the absolute count (`umi_absolute`), mean (`umi_mean`), standard deviation (`umi_sd`), CV (`umi_cv`) and FF (`umi_FF`).
* UMI probability profiles.
    - Plots the absolute probability (`p_umi_absolute`), mean (`p_umi_mean`), standard deviation (`p_umi_sd`), CV (`p_umi_cv`) and FF (`p_umi_FF`).
    - If a negative condition is present, the `p_umi_mean` plot is generated also without the negative condition as `p_umi_mean.no_neg`
* Difference profiles.
    - Plots the FF and CV difference between consecutive conditions (`deltaFF` and `deltaCV`) and between the first and last conditions (`deltaFF.extr` and `deltaCV.extr`).
* Chromosome wide plots.
    - Plots the average FF and CV on chromosome-wide windows, as `chr_wide_ff` and `chr_wide_cv`.
* Genomic overview.
    - Produces a genome overview of the read location (binned), normalized either per chromosome (`genome_view.chromo`) or globally per condition (`genome_view.global`).