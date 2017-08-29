Estimate centrality
===

The **estimate centrality** tool aims at estimating the 3D spatial nuclear centrality of genomic regions.

### v1.0.0

Two main families of metrics are implemented: probability-based and variance-based.

```
usage: ./estimate_centrality.sh [-h][-d][-s binSize][-p binStep][-g groupSize]
                                [-r prefix][-u suffix] -o outdir -c csBed
                                [BEDFILE]...

 Description:
  Estimate global centrality. The script performs the following steps:
   (1) Identify & sort chromosomes
   (2) Generate bins
   (3) Group cutsites (intersect)
   (4) Assign grouped reads to bins (intersect)
   (5) Calculate bin statistics
   (6) Combine condition into a single table
   (7) Estimate centrality
   (8) Rank bins
   (9) Write output
 
 Requirements:
  - bedtools for bin assignment
  - datamash for calculations
  - gawk for text manipulation.

 Notes:
  (A) Statistics (mean, variance) metrics take into account only cutsites
      sensed in that condition. The script ignores 'zero' loci (with no reads).
  (B) Depending on the sequencing resolution, it might not be feasible to go for
      single-cutsite resolution. Thus, cutsite can be grouped for the statistics
      calculation using the -g option.
  (C) In case of non chromosome-wide bins, the ranking is done in an ordered
      chromosome-wise manner.

 Mandatory arguments:
  -o outdir     Output folder.
  -c csBed      Cutsite bedfile.
  BEDFILE       At least two (2) GPSeq condition bedfiles, space-separated and
                in increasing order of restriction conditions intensity.
                Expected to be ordered per condition. As BEDFILE is a positional
                argument, it should be provided after any other argument.

 Optional arguments:
  -h            Show this help page.
  -d            Debugging mode: save intermediate results.
  -n            Use last condition for normalization.
  -s binSize    Bin size in bp. Default to chromosome-wide bins.
  -p binStep    Bin step in bp. Default to bin sizeinStep.
  -g groupSize  Group size in bp. Used to group bins for statistics calculation.
                binSize must be divisible by groupSize. Not used by default.
  -r prefix     Output name prefix.
  -u suffix     Output name suffix.
```
