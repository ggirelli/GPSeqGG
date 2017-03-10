
# -------- #
# SETTINGS #
# -------- #


# GENERAL
# ------------------------------

experiment='TK20'				# Experiment ID
conds='neg,1min,10min,2h'		# Comma-separated conditions list
neg='neg'						# Negative condition
numbproc=8						# Number of threads for parallelization


# ALIGNMENT
# ------------------------------

genome='hg19'					# Reference genome
aligner='bwa'					# Aligner: 'bwa' or 'bowtie2'
rmX=false						# Remove X chromosome after alignment
rmY=true						# Remove Y chromosome after alignment


# SAM filtering
# ------------------------------

cutsite='AAGCTT'				# Cutsite
mapqthr=30						# MAPQ threshold for SAM filtering
platform='L'					# Sequencing platform used
								# 	S - Sanger, Phred+33
								# 	X - Solexa, Solexa+64
								# 	I - Illumina 1.3+, Phred+64
								# 	J - Illumina 1.5+, Phred+64
								# 	L - Illumina 1.8+, Phred+33

# UMI analysis
# ------------------------------

umi_length=8					# UMI length in nt
binsize=1e6						# Window size for UMI analysis
binstep=1e5						# Distance from consecutive windows (smaller than binsize)
csrange=40						# Window around cutsites (if cutsitelist is provided)
pthr=1							# Threshold for UMI identity probability
								# Set pthr=1 to use the automatic cutoff
								# The cutoff used is the smaller between pthr and the automatic one.


# FILE PATHS
# ------------------------------
# P.s.:	If you don't want to provide a file just assign a non existent path (e.g., '-')
# 		If you assign an empty string (i.e., '') then a user input might be required

# Cutsite list file path, required columns: |chr|pos|
cutsitelist='/media/gire/Data2/BiCro-Resources/HindIII_hg19.txt'

# Maskfile path, required columns: |id|chr|start|end|
maskfile='/home/gire/Desktop/BiCro-Data/Sequencing/TK20/consensusBlacklist.centro.telo.chr5Peak.tsv'

# Chromosome lengths file path, required columns: |chr|len| (no header)
chrlengths='/media/gire/Data2/BiCro-Resources/chr_size.hg19.txt'