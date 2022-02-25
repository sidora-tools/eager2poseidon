#!/usr/bin/env Rscript

require(optparse)
library(magrittr)
require(eager2poseidon)


## Parse arguments ----------------------------
parser <- OptionParser()
parser <- add_option(parser, c("-j", "--input_janno"),
  type = "character",
  action = "store", dest = "janno_fn",
  help = "The input janno file."
)
parser <- add_option(parser, c("-e", "--eager_tsv"),
  type = "character",
  action = "store", dest = "eager_tsv_fn",
  help = "Path to the TSV file used as input for the eager run of your package data."
)
parser <- add_option(parser, c("-g", "--general_stats_table"),
  type = "character",
  action = "store", dest = "general_stats_fn",
  help = "Path to the MultiQC general stats table. Can be found in multiqc/multiqc_data/multiqc_general_stats.txt within the specified eager output directory (--outdir)."
)
parser <- add_option(parser, c("-c", "--credentials"),
  type = "character",
  action = "store", dest = "credentials",
  help = "Path to a credentials file containing four lines listing the database host, the port of the database server, user and password, respectively."
)
parser <- add_option(parser, c("-t", "--trust_uncalibrated_dates"),
  type = "logical",
  action = "store_true", default = F, dest = "trust_uncalibrated_dates",
  help = "Should any uncalibrated dates in pandora be trusted? If provided, then quickcalibrate() is used to calibrate these dates on the fly. [False]"
)
parser <- add_option(parser, c("-k", "--keep_only"),
  type = "character",
  action = "store", default = "none", dest = "keep_only",
  help = "Can be set to 'none', 'single', or 'double'. If set to 'single' or 'double', will keep only information for libraries with the specified strandedness. If 'none', all information is retained. ['none']"
)
parser <- add_option(parser, c("-s", "--snp_cutoff"),
  type = "integer",
  action = "store", default = "100", dest = "snp_cutoff",
  help = "The snp cutoff for nuclear contamination results. Nuclear contamination results with fewer than this number of SNPs will be ignored when calculating the values for 'Contamination_*' columns. [100]"
)
parser <- add_option(parser, c("-p", "--genotypePloidy"),
  type = 'character',
  action = "store", dest = "genotype_ploidy",
  metavar="ploidy",
  help = "The genotype ploidy of the genotypes produced by eager. This value will be used to fill in all missing entries in the 'Genotype_Ploidy' in the output janno file."
)
parser <- add_option(parser, c("-o", "--output_janno"),
  type = "character",
  action = "store", dest = "output_fn", default = "",
  help = "By default, the input janno is overwritten. Providing a path to this option instead writes the new janno file to the specified location."
)


args <- parse_args(parser)

## DEBUG For debugging ease
# print(args, file=stderr())

## If no output is provided, output_fn is the input janno path.
if (args$output_fn == "") {
  output_fn <- args$janno_fn
} else {
  output_fn <- args$output_fn
}

input_janno_table <- standardise_janno(args$janno_fn)

external_results_table <- collate_external_results(
  sample_ids = dplyr::select(input_janno_table, Poseidon_ID),
  eager_tsv_fn = args$eager_tsv_fn,
  general_stats_fn = args$general_stats_fn,
  credentials = args$credentials,
  keep_only = args$keep_only,
  trust_uncalibrated_dates = args$trust_uncalibrated_dates,
  snp_cutoff = args$snp_cutoff
)

output_janno <- fill_in_janno(
  input_janno_table = input_janno_table,
  external_results_table = external_results_table,
  genotype_ploidy = args$genotype_ploidy
)

## Replace NAs with "n/a" (requires re-typing everything to character)
output_janno <- output_janno %>%
  dplyr::mutate(
    dplyr::across(
      .fn = as.character
    )
  ) %>%
  replace(is.na(.), "n/a") %>%
  poseidonR::as.janno()

## Write output file
poseidonR::write_janno(output_janno, output_fn)
