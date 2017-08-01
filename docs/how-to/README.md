HOW-TO run gpseq-seq-gg
===

## Prepare your data

### 1. Input folder

Create an **input** folder, that will contain the pattern file and the sequencing output of your experiment(s).

### 2. Sequencing data

Place the sequencing data in the **input** folder.

* If single-end sequencing, only one fastq.tar.gz file is expected (i.e., R1).
* If pair-end sequencing also a second fastq.tar.gz file is expected (i.e., R2).

The sequencing output file name (i.e., the fastq file(s)) **must** follow the [Illumina](https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/NamingConvention_FASTQ-files-swBS.htm) format `expID_S*_L***_R*_001.fastq.gz`, where `expID` is the *experiment ID* (this is provided to the pipeline through the pattern file), `S*` is the sample number, `L***` is the lane number, and `R*` indicates which fragment is contained.

If multiple lanes are available in the output, the file should be merged using the command:

```
cat expID_S1_L*_R1_001.fastq.gz > expID_S1_LALL_R1_001.fastq.gz
```

If the fastq files were merged, remove the original ones. Remember to double check that the merged fastq file still follows the Illumina name format.

Please, note that the pipeline will be able to recognize the fastq files **only and only if** they contain the `expID` **AND** either "R1" or "R2" in their name.

### 3. Pattern file

Create a pattern file called `patterns.tsv` containing experiment id, condition name, pattern, and linker length in nucleotides. One condition per row, with columns separated **by tabulations**.

```
experiment_ID condition_label scan_for_matches_pattern  trimming_length
```

An example being the following:

```
TK20  neg ^ 8...8 CATCAGAA AAGCTT 1...1000 $  22
TK20  1min  ^ 8...8 CATCATCC 1...1000 $ 18
TK21  1min  ^ 8...8 GTCGTTCC 1...1000 $ 16
TK21  2h  ^ 8...8 TGATGTCC AAGCTT 1...1000 $  22
```

The pattern field (3rd) must be compatible with `scan_for_matches`, more information in the [official README](http://iubio.bio.indiana.edu/soft/molbio/pattern/scan_for_matches.readme).

**The pattern file should NOT have a header line.**

## Select settings

Before running the pipeline, a few settings should be prepared.

```
  -w  Perform every step of the pipeline without asking.
  -x  Remove X chromosome after alignment.
  -y  Remove Y chromosome after alignment.
  -t threads  Number of threads for parallelization.
  -g refGenome  Path to reference genome file. Default: 'hg19'.
  -a aligner  Aligner. Either 'bwa' (default) or 'bowtie2'.
  -d bwaIndex Path to BWA index file. Required if BWA is the selected aligner.
  -q mapqThr  Mapping quality threshold. Default: 30.
  -p platform Sequencing platform. Default: 'L'.
  -u umilength  UMI sequence length. Default: 8.
  -r csRange  Range around cutsite for UMI assignment. Default: 40.
  -j emax Maximum error probability for read quality filtering. Default: 1e-3.
  -k eperc  Maximum % of bases with emax error probability. Default: 20.
  -z binSize  Bin size. Default: 1e6.
  -b binStep  Bin step. Default: 1e5.
  -c cutsite  Cutsite sequence, needed for read re-position. Default: AAGCTT.
  -l csList File with cutsite list. Columns: chr|pos. No header.
  -m maskFile File with masked regions. Columns: id|chr|start|end. No header.
  -s chrLengths File with chromosome lengths. chr|len. No header.
  -n neatness Neatness level: 0 (dirty), 1 (neat), 2 (neatest). Default: 1.
```

You don't need to specify a setting if the default value is compatible with your analysis. Also, the selected settings will be showed and a confirmation is required before the pipeline starts. Thus, we suggest to add the settings one by one, and start the program to check what is missing and if everything is correct.

## Run

To run the pipeline, open a terminal and go to the `gpseq-seq-gg` folder. Type `./main.sh -h` to see the help page and instructions on how to run the pipeline. The script guides the user into building a proper run, by specifying which options are missing (if any) and eventual mistakes in the pattern file. A description of possible errors is available in the [F.A.Q.s](../faq/).