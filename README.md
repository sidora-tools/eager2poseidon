# eager2poseidon
A tool to fill in the janno file of a poseidon package by pulling archaeological information from Pandora, 
and genetic information from an nf-core/eager run (prepared with pandora2eager).

### eager2poseidon is designed to only fill in empty fields (`n/a`) in the input janno file. Any fields in the input janno that have been previously filled in will NOT be edited/updated.

# Installation instructions
To install the package, open an R session and run the following code.
```
if(!require('remotes')) install.packages('remotes')
remotes::install_github("sidora-tools/eager2poseidon")
```

# Quickstart
The Rscript `eager2poseidon.R` can be ran directly from the command line. The helptext provided with `-h` 
explains all the arguments required for the script. Default values to optional arguments are shown in square 
brackets `[]` following the argument explanation.
```
$ ./exec/eager2poseidon.R -h
Loading required package: optparse
Loading required package: eager2poseidon
Usage: ./exec/eager2poseidon.R [options]


Options:
	-h, --help
		Show this help message and exit

	-j INPUT_JANNO, --input_janno=INPUT_JANNO
		The input janno file.

	-e EAGER_TSV, --eager_tsv=EAGER_TSV
		Path to the TSV file used as input for the eager run of your package data.

	-g GENERAL_STATS_TABLE, --general_stats_table=GENERAL_STATS_TABLE
		Path to the MultiQC general stats table. Can be found in multiqc/multiqc_data/multiqc_general_stats.txt within the specified eager output directory (--outdir).

	-c CREDENTIALS, --credentials=CREDENTIALS
		Path to a credentials file containing four lines listing the database host, the port of the database server, user and password, respectively.

	-t, --trust_uncalibrated_dates
		Should any uncalibrated dates in pandora be trusted? If provided, then quickcalibrate() is used to calibrate these dates on the fly. [False]

	-k KEEP_ONLY, --keep_only=KEEP_ONLY
		Can be set to 'none', 'single', or 'double'. If set to 'single' or 'double', will keep only information for libraries with the specified strandedness. If 'none', all information is retained. ['none']

	-s SNP_CUTOFF, --snp_cutoff=SNP_CUTOFF
		The snp cutoff for nuclear contamination results. Nuclear contamination results with fewer than this number of SNPs will be ignored when calculating the values for 'Contamination_*' columns. [100]

	-p PLOIDY, --genotypePloidy=PLOIDY
		The genotype ploidy of the genotypes produced by eager. This value will be used to fill in all missing entries in the 'Genotype_Ploidy' in the output janno file.

	-o OUTPUT_JANNO, --output_janno=OUTPUT_JANNO
		By default, the input janno is overwritten. Providing a path to this option instead writes the new janno file to the specified location.

```
