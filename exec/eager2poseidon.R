#!/usr/bin/env Rscript

require(optparse)
library(magrittr)

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
parser <- add_option(parser, c("-p", "--prefer"),
  type = "character",
  action = "store", default = "none", dest = "prefer",
  help = "Can be set to 'none', 'single', or 'double'. If set to 'single' or 'double', will keep only information for libraries with the specified strandedness. If 'none', all information is retained. [none]"
)
parser <- add_option(parser, c("-o", "--output_janno"),
  type = "character",
  action = "store", dest = "output_fn", default = "",
  help = "By default, the input janno is overwritten. Providing a path to this option instead writes the new janno file to the specified location."
)


args <- parse_args(parser)

## If no output is provided, output_fn is the input janno path.
if (args$output_fn == "") {
  output_fn <- args$janno_fn
}

input_janno_table <- standardise_janno(args$input)

external_results_table <- collate_external_results(
  sample_ids = dplyr::select(input_janno_table, Poseidon_ID),
  eager_tsv_fn = args$eager_tsv_fn,
  general_stats_fn = args$general_stats_fn,
  credentials = args$credentials,
  prefer = args$prefer,
  trust_uncalibrated_dates = args$trust_uncalibrated_dates
)

output_janno <- fill_in_janno(
  input_janno = input_janno_table,
  external_results_table = external_results_table
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
