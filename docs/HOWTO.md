HOW-TO run gpseq-seq-gg
===

## Prepare your data

### 1. Input folder

Create an **input** folder, that will contain the pattern file and the sequencing output.

### 2. Sequencing data

Place the sequencing output in the **input** folder.

* If single-end sequencing, only one fastq.tar.gz file is expected (i.e., R1).
* If pair-end sequencing also a second fastq.tar.gz file is expected (i.e., R2).

The sequencing output file name must follow the [Illumina](https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/NamingConvention_FASTQ-files-swBS.htm) format `expID_S*_L***_R*_001.fastq.gz`, where `expID` is the *experiment ID*, `S*` is the sample number, `L***` is the lane number, and `R*` indicates which fragment is contained.

If multiple lane are available in the output, the file should be merged using the command:

```
cat expID_S1_L*_R1_001.fastq.gz > expID_S1_LALL_R1_001.fastq.gz
```

If the fastq files were merged, remove the original ones.

### 3. Pattern file

Create a pattern file called `pat_files` containing the condition name, pattern, cutsite sequence and linker length in nucleotides, per row and separated by tabulations.

For example:

```
 neg ^ 8...8 CATCAGAA AAGCTT 1...1000 $  AAGCTT  22
 1min    ^ 8...8 CATCATCC AAGCTT 1...1000 $  AAGCTT  22
 10min   ^ 8...8 GTCGTTCC AAGCTT 1...1000 $  AAGCTT  22
 2h  ^ 8...8 TGATGTCC AAGCTT 1...1000 $  AAGCTT  22
```

More information on the format of the pattern file is available in the README.

The pattern information must be compatible with `scan_for_matches`, more information in the [official README](http://iubio.bio.indiana.edu/soft/molbio/pattern/scan_for_matches.readme).

## Select settings

Before running the pipeline, a few settings should be prepared.

```
  -i indir      Input directory.
  -o outdir     Output directory. Created if not found.
  -e expID      Experiment ID.
  -h            Show the help page.
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

## Run

To run the pipeline, open a terminal and go to the `gpseq-seq-gg` folder. Type `./main.sh -h` to see the help page and instructions on how to run the pipeline.

The selected settings will be showed and a confirmation is required before the pipeline starts.