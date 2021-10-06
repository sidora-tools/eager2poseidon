#!/usr/bin/env Rscript
if (!require('sidora.core')) {
  if(!require('remotes')) install.packages('remotes')
remotes::install_github('sidora-tools/sidora.core', quiet=T)
} else {library(sidora.core)}
# if (!require(tidyverse)) {install.packages('tidyverse')}
require('dplyr')
require('readr')
require('purrr')
require('optparse')

## Validate genotype source option input
validate_geno_source <- function(option, opt_str, value, parser) {
  valid_entries=c("single", "double")
  ifelse(value %in% valid_entries, return(value), stop(call.=F, "\nInvalid genotype source: '", value, "'\nAccepted values: ", paste(valid_entries,collapse=", "),"\n\n"))
}

## Infer genetic sex from x-/y-rates. Caveat, no error bar taken into account.
infer_genetic_sex <- function(x_rate, y_rate) {
    ifelse(`x_rate` > 0.70 & `y_rate` < 0.10, "F", 
        ifelse(`x_rate` < 0.6 & `x_rate` > 0.3 & `y_rate` > 0.3 & `y_rate` < 0.6, "M",
            "U")
    )
}

## Infer row type by ID length and TF/SG
# infer_data_type <- function()
  
## MAIN ##

## Parse arguments ----------------------------
parser <- OptionParser()
parser <- add_option(parser, c("-i", "--input"), type = 'character', 
                     action = "store", dest = "input", 
                     help = "A path to the eager output 
                    directory to turn into a poseidon package.")
parser <- add_option(parser, c("-n", "--outPackageName"), type = 'character',
                     action = "store", dest = "package_name", 
                     help = "The output package name")
parser <- add_option(parser, c("-o", "--outPackagePath"), type = 'character',
                     action = "store", dest = "package_path", 
                     help = "The output package directory path")
parser <- add_option(parser, c("-g", "--genotypeSource"), type = 'character',
                     action = "callback", dest = "genotype_source", 
                     callback = validate_geno_source,
                     help = "If both single and double stranded genotypes exist in the eager output, which set should be converted into a poseidon package?")
args <- parse_args(parser)


# eager_outputs <- read_tsv(args$input, col_names = "paths", col_types = 'c')
# general_stats_tables <- eager_outputs$paths %>% purrr::map_chr(~ paste0(., "/multiqc/multiqc_data/multiqc_general_stats.txt"))
# con <- get_pandora_connection(cred_file = args[2])

## Load eager output data
eager_data <- read_tsv(paste0(args$input, "/multiqc/multiqc_data/multiqc_general_stats.txt"), na = c("", "N/A", "NA")) 
eager_data <- eager_data %>% 
  ## Only Keep columns of interest
  select(
    `Sample`, 
    `SexDetErrmine_mqc-generalstats-sexdeterrmine-RateX`, 
    `SexDetErrmine_mqc-generalstats-sexdeterrmine-RateY`, 
    `endorSpy_mqc-generalstats-endorspy-endogenous_dna`, 
    `snp_coverage_mqc-generalstats-snp_coverage-Covered_Snps`, 
    `snp_coverage_mqc-generalstats-snp_coverage-Total_Snps`, 
    `DamageProfiler_mqc-generalstats-damageprofiler-5_Prime1`, 
    `nuclear_contamination_mqc-generalstats-nuclear_contamination-Method1_ML_estimate`, 
    `nuclear_contamination_mqc-generalstats-nuclear_contamination-Method1_ML_SE`
  ) %>% 
  ## give columns better names
  rename(
    rate_x=`SexDetErrmine_mqc-generalstats-sexdeterrmine-RateX`,
    rate_y=`SexDetErrmine_mqc-generalstats-sexdeterrmine-RateY`,
    endogenous_dna=`endorSpy_mqc-generalstats-endorspy-endogenous_dna`,
    convered_snps=`snp_coverage_mqc-generalstats-snp_coverage-Covered_Snps`,
    total_snps=`snp_coverage_mqc-generalstats-snp_coverage-Total_Snps`, 
    damage_5p1=`DamageProfiler_mqc-generalstats-damageprofiler-5_Prime1`,
    x_contamination=`nuclear_contamination_mqc-generalstats-nuclear_contamination-Method1_ML_estimate`,
    x_contamination_error=`nuclear_contamination_mqc-generalstats-nuclear_contamination-Method1_ML_SE`
  ) %>%
  ## Remove rows with all but the sample name being NA (sequencing lanes)
  filter(rowSums(is.na(.)) != ncol(.) - 1) %>%
  ##Set Xcont to 0 if negative
  mutate(
    x_contamination=case_when(
      x_contamination < 0 ~ 0,
      x_contamination < 1 ~ x_contamination
      ## Sometimes Angsd returns contamination over 1 when no SNPs pass the depth filter.
    )
  )

