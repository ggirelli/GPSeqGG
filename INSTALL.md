# (1) Get source code

Clone the repository and git submodules to a local folder using:

```
git clone http://github.com/ggirelli/gpseq-seq-gg/
cd gpseq-seq-gg
git submodule init
git submodule update
```

# (2) Third party software

* `bwa` v0.7.12-r1039
* `bowtie2` v2.2.6
* `fastqc` v0.11.4
* `scan_for_matches`
* `samtools` v1.3.1

On Ubuntu, run: `sudo apt install bwa bowtie2 fastqc`

### scan_for_matches

Download scan_for_matches from [here](http://www.theseed.org/servers/downloads/scan_for_matches.tgz). Then, follow [the instructions](http://blog.theseed.org/servers/2010/07/scan-for-matches.html) to install it. Finally, add the path of the uncompressed folder to your `$PATH` environment variable with:

```
echo -e "\nexport PATH=$PATH:path_to_scan_for_matches" >> ~/.bash_profile
```

### fastqc

If, after installation, fastqc complains of missing files, download them from [here](https://github.com/ggirelli/configs/tree/master/fastqc/Configuration).

# (3) Required R libraries

* argparser
* data.table
* ggplot2
* grid
* gridExtra
* parallel
* readr
* reshape2

Install them from within R with the `install.packages()` command.
Or run the `./INSTALL.R` script.

# (4) Required Python libraries

* argparse
* sys

Should be pre-installed with Python. Install them with `pip install` otherwise.

# (5) FOLDERS

Need a ~/tmp folder for sorting.

```
mkdir ~/tmp
```
