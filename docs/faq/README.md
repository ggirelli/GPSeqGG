F.A.Q.
===

* What do the `summary` columns mean?
* What can I do if I get the following "!!! ERROR!" message?
	 * [Invalid -i option, folder not found.](#invalid--i-option-folder-not-found)
	 * [Invalid -a option. Available values: 'bwa', 'bowtie2'.](#invalid--a-option-available-values-bwa-bowtie2)
	 * [Missing mandatory -d option. OR Invalid -d option, file not found.](#missing-mandatory--d-option-or-invalid--d-option-file-not-found)
	 * [Invalid -l option, file not found.](#invalid--l-option-file-not-found)
	 * [Invalid -m option, file not found.](#invalid--m-option-file-not-found)
	 * [Invalid -s option, file not found.](#invalid--s-option-file-not-found)
	 * [Missing mandatory -i option.](#missing-mandatory--i-option)
	 * [Missing mandatory -o option.](#missing-mandatory--o-option)
	 * [Missing mandatory -e option.](#missing-mandatory--e-option)
	 * [Missing columns in pat_files, row X. Expected 4 columns, found Y.](#missing-columns-in-pat_files-row-x-expected-4-columns-found-y)
	 * [Missing pat_files.](#missing-pat_files)
	 * [No R1/R2 files found for $expID.](#no-r1r2-files-found-for-expid)

## What do the `summary` columns mean?

Visit the [summary description page](../summary/) for more details.

## What can I do if I get the following "!!! ERROR!" message?

### Invalid -i option, folder not found.

The specified **input folder** does not exist. Double-check if you provided the proper path.

### Invalid -a option. Available values: 'bwa', 'bowtie2'.

Currently, `gpseq-seq-gg` supports only `bwa mem` (with `bwa` option) and `bowtie2` (with `bowtie2` option) as aligners. `bwa` is the default.

### Missing mandatory -d option. OR Invalid -d option, file not found.

The specified **BWA index file** does not exist. Double-check if you provided the proper path. Remember that the `-d` option is mandatory if you use `-a bwa` (default).

### Invalid -l option, file not found.

The specified **cutsite list file** does not exist. Double-check if you provided the proper path. The cutsite list file is expected to have two tab-separated columns with the chromosome ('chrN' format) and the position, with no header.

```
chr1	214151
chr2	24125
chr2	65643
```

### Invalid -m option, file not found.

The specified **mask file** does not exist. Double-check if you provided the proper path. The mask file is expected to have one line per masked region and four tab-separated columns, with no header: chromosome ID, chromosome name ('chrN' format), start and end positions.

```
1	chr1	10000	120000
X	chrX	350	4500
```

### Invalid -s option, file not found.

The specified **chromosome length file** does not exist. Double-check if you provided the proper path. The chromosome length file is expected to have one line per chromosome and two tab-separated columns, with no header: chromsome ('chrN' format) and chromosome length.

```
chr1    249250621
chr2    243199373
chr3    198022430
```

### Missing mandatory -i option.

The -i option is mandatory, as the pipeline needs input data to analyze.

### Missing mandatory -o option.

The -o option is mandatory, as the pipeline needs to know where to save the output.

### Missing mandatory -e option.

The -e option is necessary, as it is used to identify the fastq files and to generate proper labeling of plots and in the reports. The sequencing output file name (i.e., the name of the fastq file(s)) must follow the [Illumina](https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/NamingConvention_FASTQ-files-swBS.htm) format `expID_S*_L***_R*_001.fastq.gz`, where `expID` is the *experiment ID* provided with the -e option, `S*` is the sample number, `L***` is the lane number, and `R*` indicates which fragment is contained.

### Missing columns in pat_files, row X. Expected 4 columns, found Y.

The pattern file (pat_files) format expects one line per condition and 4 tab-separated columns.For example:

```
 neg	^ 8...8 CATCAGAA AAGCTT 1...1000 $	AAGCTT	22
 1min	^ 8...8 CATCATCC AAGCTT 1...1000 $	AAGCTT	22
 10min	^ 8...8 GTCGTTCC AAGCTT 1...1000 $	AAGCTT	22
 2h	^ 8...8 TGATGTCC AAGCTT 1...1000 $	AAGCTT	22
```

If you get this error you can double-check that:

* columns should be tab-separated
* no empty row should be present at the end of the file

### Missing pat_files.

No `pat_files` is present in the input folder. Also, a mock pattern file with the expected structure should be printed after the error message.

### No R1/R2 files found for $expID.

Either no fastq file is present in the specified input folder, or they do not follow the [Illumina](https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/NamingConvention_FASTQ-files-swBS.htm) format `expID_S*_L***_R*_001.fastq.gz`, where `expID` is the *experiment ID* provided with the -e option, `S*` is the sample number, `L***` is the lane number, and `R*` indicates which fragment is contained.
