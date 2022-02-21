if (getRversion() >= "2.15.1") utils::globalVariables(c(".")) ## Disables notes about '.' due to magrittr

#' Collect poseidon janno data from an eager TSV and its corresponding general stats table.
#'
#' @param eager_tsv character. Path to the TSV file used as input for the eager run.
#' @param general_stats_fn character. Path to the root eager output dir (`--outdir`) of the eager run.
#' @param prefer character. Can be set to "none", "single", or "double".
#' If set to 'single; or 'double', will keep only information for libraries with the specified strandedness.
#' If 'none', all information is retained.
#'
#' @return A tibble containing the poseidon janno fields of Nr_Libraries, Capture_Type, UDG, Library_Built, Damage, Nr_SNPs, Endogenous_DNA, Contamination,
#' Contamination_Err, Contamination_Meas, and Contamination_Note
#'
#' @export
#' @importFrom rlang .data
#' @importFrom magrittr "%>%"

import_eager_results <- function(eager_tsv, general_stats_fn, prefer) {
  ## Parse TSV
  tsv_data <- read_eager_tsv(eager_tsv, prefer)
  write("Parsing of eager input TSV completed.", file = stderr())

  ## Parse general stats table
  eager_data <- read_eager_stats_table(general_stats_fn, tsv_data)
  write("Parsing of eager general stats table completed.", file = stderr())

  ## Collate data
  result <- tsv_data %>%
    dplyr::group_by(.data$Sample_Name) %>%
    dplyr::summarise(
      # ## Keep the library IDs as they appear in the tsv, to filter for relevant lines from the general stats table
      # Libraries_Included=paste(Library_ID,collapse=";"),
      Nr_Libraries = dplyr::n(),
      Capture_Type = paste(.data$Seq_Type, collapse = ";"),
      ## UDG gets formatted here to account for 'mixed' UDG treatment when multiple libs of the same strandedness have different UDG treatment.
      UDG = paste(.data$UDG_Treatment, collapse = ";") %>% format_for_poseidon(., "UDG"),
      Library_Built = paste(.data$Strandedness, collapse = ";")
    ) %>%
    dplyr::full_join(eager_data, by = c("Sample_Name" = "Sample"))

  result
}

#' Read content of an eager TSV and convert information to poseidon janno columns
#'
#' @inheritParams import_eager_results
#'
#' @return A tibble with columns for information for janno columns Sample_Name, Nr_Libraries, Capture_Type, UDG,
#' and Library_Built, and an additional column with the IDs of the libraries included for each sample. The latter
#' is used to pull only relevant information from the eager results tables.
#'
#' @export
#' @importFrom magrittr "%>%"

read_eager_tsv <- function(eager_tsv, prefer = "none") {
  ## End execution if TSV file is not found
  if (!file.exists(eager_tsv)) {
    stop(paste0("File '", eager_tsv, "' not found."))
  }
  if (!prefer %in% c("none", "single", "double")) {
    stop(paste0("Invalid library strandedness preference: '", eager_tsv, "'"))
  }
  tsv_data <- readr::read_tsv(eager_tsv, col_types = "cciiccccccc") %>%
    dplyr::select(.data$Sample_Name, .data$Library_ID, .data$Strandedness, .data$UDG_Treatment, .data$R1, .data$R2, .data$BAM)
  ## Filter for preferred library strandedness
  if (prefer != "none") {
    tsv_data <- dplyr::filter(tsv_data, .data$Strandedness == prefer)
  }
  ## Separate the Library ID to pandora Library ID and Sequencing ID if the latter was used
  tsv_data <- tsv_data %>%
    dplyr::mutate(
      Library = substr(.data$Library_ID, 1, 12),
      Seq_Type = dplyr::case_when(
        ## First try to infer from Lib_ID
        substr(.data$Library_ID, 14, 15) != "" ~ substr(.data$Library_ID, 14, 15) %>% format_for_poseidon(., "Capture_Type"),
        ## Then check R1, if provided. Assumes R1 file name starts with pandora full sequencing ID
        !is.na(.data$R1) ~ basename(.data$R1) %>%
          substr(., 14, 15) %>%
          format_for_poseidon(., "Capture_Type"),
        ##  Else, check the BAM input (must be there if no R1 is provided)
        !is.na(.data$BAM) ~ basename(.data$BAM) %>%
          substr(., 14, 15) %>%
          format_for_poseidon(., "Capture_Type"),
        ## Finally, a catch-all, just in case we have unexpected behaviours
        TRUE ~ NA_character_
      ),
      Strandedness = format_for_poseidon(.data$Strandedness, "Library_Built")
    ) %>%
    dplyr::select(-c(.data$R1, .data$R2, .data$BAM)) %>%
    dplyr::distinct()
  tsv_data
}


#' Internal function to convert information from eager lingo to poseidon accepted strings
#'
#' @param x character. A sequencing ID.
#' @param field character. The field to format for. One of 'UDG', 'Library_Built', or 'Capture_Type'.
#'
#' @return character. A valid Poseidon Capture_Type, Library_Built, or UDG value.
#'
format_for_poseidon <- function(x, field) {
  if (field == "Capture_Type") {
    result <- dplyr::case_when(
      x == "SG" ~ "Shotgun",
      x == "TF" ~ "1240k",
      TRUE ~ "Other"
    )
  }

  if (field == "Library_Built") {
    result <- dplyr::case_when(
      x == "single" ~ "ss",
      x == "double" ~ "ds",
      TRUE ~ "other"
    )
  }

  if (field == "UDG") {
    result <- dplyr::case_when(
      x == "none" ~ "minus",
      x == "half" ~ "half",
      x == "full" ~ "plus",
      TRUE ~ "mixed"
    )
  }
  result
}

#' Extract poseidon information from the general stats table of a run
#'
#' @inheritParams import_eager_results
#'
#' @param tsv_data A tibble containing information from the eager TSV. Created with \link[eager2poseidon:read_eager_tsv]{read_eager_tsv}
#'
#' @return A tibble containing the poseidon janno fields of Damage, Nr_SNPs, Endogenous_DNA, Contamination,
#' Contamination_Err, Contamination_Meas, and Contamination_Note
#' @export
#' @importFrom magrittr "%>%"

read_eager_stats_table <- function(general_stats_fn, tsv_data) {
  ## Check the file exists
  if (!file.exists(general_stats_fn)) {
    stop(paste0("File '", general_stats_fn, "' not found."))
  }

  ## Create a character vector of the sample and library IDs to filter out information from the general stats table
  pull_samples <- c(unique(tsv_data$Sample_Name), tsv_data$Library_ID) %>% sort()

  ## Load stats table and keep relevant columns and lines.
  general_stats <- readr::read_tsv(general_stats_fn, na = c("", "N/A", "NA")) %>%
    dplyr::filter(.data$Sample %in% pull_samples) %>%
    dplyr::select(
      .data$Sample,
      endogenous_dna = .data$`endorSpy_mqc-generalstats-endorspy-endogenous_dna`,
      covered_snps = .data$`snp_coverage_mqc-generalstats-snp_coverage-Covered_Snps`,
      total_snps = .data$`snp_coverage_mqc-generalstats-snp_coverage-Total_Snps`,
      damage_5p1 = .data$`DamageProfiler_mqc-generalstats-damageprofiler-5_Prime1`,
      x_contamination_snps = .data$`nuclear_contamination_mqc-generalstats-nuclear_contamination-Num_SNPs`,
      x_contamination = .data$`nuclear_contamination_mqc-generalstats-nuclear_contamination-Method1_ML_estimate`,
      x_contamination_error = .data$`nuclear_contamination_mqc-generalstats-nuclear_contamination-Method1_ML_SE`,
      duplication_rate = .data$`Picard_mqc-generalstats-picard-PERCENT_DUPLICATION`,
      filtered_mapped_reads = .data$`Samtools Flagstat (post-samtools filter)_mqc-generalstats-samtools_flagstat_post_samtools_filter-mapped_passed`
    ) %>%
    dplyr::mutate(
      x_contamination = dplyr::case_when(
        .data$x_contamination_snps == 0 ~ "n/a",
        # is.na(x_contamination) ~ "n/a",
        TRUE ~ format(.data$x_contamination, nsmall = 3, digits = 0, trim = T)
      ),
      x_contamination_error = tidyr::replace_na("n/a")
    )

  ## For endogenous, only look at Shotgun libraries. If none exist return "n/a"s
  if (nrow(tsv_data %>% dplyr::filter(.data$Seq_Type == "Shotgun")) == 0) {
    samples <- unique(tsv_data$Sample_Name)
    endogenous <- rep("n/a", length(samples))
    endogenous_per_sample <- tibble::tibble(Sample = samples, Endogenous_DNA = endogenous)
  } else {
    endogenous_per_sample <- dplyr::right_join(general_stats, tsv_data %>% dplyr::filter(.data$Seq_Type == "Shotgun"), by = c("Sample" = "Library_ID")) %>%
      tidyr::separate(.data$Sample, sep = "\\.", into = c("Sample", "Library"), fill = "right", extra = "merge") %>%
      dplyr::group_by(.data$Sample) %>%
      dplyr::summarise(
        Endogenous_DNA = format(max(.data$endogenous_dna), nsmall = 2, digits = 0)
      )
  }

  ## Library info for library based fields
  contamination_per_sample <- general_stats %>%
    tidyr::separate(.data$Sample, sep = "\\.", into = c("Sample", "Library"), fill = "right", extra = "merge") %>%
    dplyr::filter(!is.na(.data$Library)) %>%
    ## Keep only Library Entries
    dplyr::group_by(.data$Sample) %>%
    dplyr::summarise(
      Contamination = paste(.data$x_contamination, collapse = ";"),
      Contamination_Err = paste(.data$x_contamination_error, collapse = ";"),
      Contamination_Meas = paste(rep("ANGSD", length(.data$x_contamination)), collapse = ";"),
      Contamination_Note = paste("Nr Snps:", paste(.data$x_contamination_snps, collapse = ";"))
    )

  ## For the rest use everything
  sample_stats <- general_stats %>%
    tidyr::separate(.data$Sample, sep = "\\.", into = c("Sample", "Library"), fill = "right", extra = "merge") %>%
    dplyr::group_by(.data$Sample) %>%
    dplyr::summarise(
      Damage = stats::weighted.mean(stats::na.omit(.data$damage_5p1), stats::na.omit(.data$filtered_mapped_reads)),
      Nr_SNPs = max(.data$covered_snps, na.rm = T)
    )

  ## Join it all together
  result <- dplyr::full_join(sample_stats, endogenous_per_sample, by = "Sample") %>%
    dplyr::full_join(contamination_per_sample, by = "Sample")
  result
}
