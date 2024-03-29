if (getRversion() >= "2.15.1") utils::globalVariables(c(".")) ## Disables notes about '.' due to magrittr

#' Compile eager results tables
#'
#' Compile together various tables of eager results read into poseidon-compatible tables.
#'
#' @param tsv_table A tibble with the eager TSV input information as returned by \link[eagerR]{infer_merged_bam_names}.
#' @param sexdet_table A tibble with the eager sex determination information as returned by \link[eagerR]{read_sexdet_json}.
#' @param snpcov_table A tibble with the eager eigenstrat snp coverage information as returned by \link[eagerR]{read_snp_coverage_json}.
#' @param dmg_table A tibble with the eager damageprofiler results, as returned by \link[eagerR]{read_damageprofiler_jsons_from_dir}.
#' @param endogenous_table A tibble with the eager endogenous DNA information as returned by \link[eagerR]{read_endorspy_jsons_from_dir}.
#' @param nuccont_table A tibble withthe eager nuclear contamination results as returned by \link[eagerR]{read_angsd_cont_json}.
#' @param contamination_method The contamination method to keep. Can be either "1" or "2".
#' @param contamination_algorithm The estimation algorithm to keep. Can be "ml" or "mom".
#' @param capture_type character. The Capture_Type. Assumes all libraries have the same one. Defaults to NA, which will disable inference of the Capture_Type and return a dummy column.
#' @inheritParams infer_genetic_sex
#'
#' @return A tibble with the compiled information from the provided eager output tables.
#' @export
#'
compile_eager_result_tables <- function(tsv_table=NULL, sexdet_table=NULL, snpcov_table=NULL, dmg_table=NULL, endogenous_table=NULL, nuccont_table=NULL, contamination_method="1", contamination_algorithm = "ml", XX_cutoffs = c(0.7, 1.2, 0.0, 0.1), XY_cutoffs = c(0.2, 0.6, 0.3, 0.6), capture_type = NA_character_) {

  ## Without a TSV info table, nothing can be compiled
  if (is.null(tsv_table)) {stop("[compile_eager_result_tables()]: tsv_table is required.")}

  ## Validate input params
  if (! contamination_method %in% c("1","2")) { stop(paste0("[compile_eager_result_tables()]: Invalid preferred nuclear contamination estimation method provided.\n\tValid entries are c('1','2').\n\tYou provided: ",contamination_method))}

  if (! contamination_algorithm %in% c("ml","mom")) { stop(paste0("[compile_eager_result_tables()]: Invalid preferred nuclear contamination estimation method provided.\n\tValid entries are c('ml','mom').\n\tYou provided: ", contamination_algorithm))}

  if (! capture_type %in% c("Shotgun", "1240K", "OtherCapture", "ReferenceGenome", NA_character_)) {stop(paste0("[compile_eager_result_tables()]: Invalid capture type provided.\n\tValid entries are c('Shotgun', '1240K', 'OtherCapture', 'ReferenceGenome', NA_character_).\n\tYou provided: ", capture_type))}

  ###################
  ## TSV Lib infos ##
  ###################

  ## Add number of libraries, capture type, overall UDG treatment and Strandedness columns
  tsv_table <- tsv_table %>%
    dplyr::select(.data$Sample_Name, .data$Library_ID, .data$Strandedness, .data$UDG_Treatment) %>%
    dplyr::group_by(.data$Sample_Name, .data$Strandedness) %>%
    dplyr::summarise(.groups='keep',
                     UDG=dplyr::case_when(
                       unique(.data$UDG_Treatment) %>% length(.) > 1 ~ 'mixed',
                       TRUE ~ format_for_poseidon(unique(.data$UDG_Treatment), "UDG")
                     ),
                     Nr_Libraries=vctrs::vec_unique_count(.data$Library_ID),
                     Capture_Type=dplyr::if_else(
                       is.na(capture_type),
                       ##TRUE
                       NA_character_,
                       ## FALSE
                       paste0(rep(capture_type, .data$Nr_Libraries), collapse=";")
                     ),
                     Library_Built=format_for_poseidon(.data$Strandedness, "Library_Built")
    ) %>%
    dplyr::left_join(tsv_table, ., by=c("Sample_Name", "Strandedness")) %>%
    ## Keep one line per library, otherwise samples with multiple libraries get duplicate entries per library
    dplyr::distinct()


  ############
  ## SEXDET ##
  ############
  if (is.null(sexdet_table)) {
    ## If no sexdet table was provided, create dummy columns so output table columns are consistent
    tsv_table <- tsv_table %>%
      tibble::add_column(
        "Genetic_Sex"=NA_character_,
        "Sex_Determination_Note"=NA_character_
      )
  } else {
    ## Use sexdet results to infer genetic sex based on provided cutoffs, and put together the associated Note field.
    sexdet <- sexdet_table %>%
      dplyr::group_by(.data$sexdet_input_bam) %>%
      dplyr::mutate(
      Genetic_Sex = infer_genetic_sex(.data$sexdet_ratex, .data$sexdet_ratey, XX_cutoffs, XY_cutoffs),
      Sex_Determination_Note = paste0(
          "x-rate: ",
          format(.data$sexdet_ratex, nsmall = 3, digits = 1, trim = T),
          " +- ",
          format(.data$sexdet_rateerrx, nsmall = 3, digits = 1, trim = T),
          ", y-rate: ",
          format(.data$sexdet_ratey, nsmall = 3, digits = 1, trim = T),
          " +- ",
          format(.data$sexdet_rateerry, nsmall = 3, digits = 1, trim = T)
        )
      ) %>%
      dplyr::select(.data$sexdet_input_bam, .data$Genetic_Sex, .data$Sex_Determination_Note) %>%
      dplyr::ungroup()
    tsv_table <- dplyr::left_join(tsv_table, sexdet, by=c("sexdet_bam"="sexdet_input_bam"))
  }

  #############
  ## SNP COV ##
  #############

  if (is.null(snpcov_table)) {
    ## If no snp coverage table was provided, create dummy columns so output table columns are consistent
    tsv_table <- tsv_table %>%
      tibble::add_column(
        "Nr_SNPs"=NA_integer_
      )
  } else {
    ## Add number of snps by sample/strandedness combo
    snpcov_table <- snpcov_table %>% dplyr::select(.data$snpcov_sample_name, "Nr_SNPs"=.data$snpcov_covered_snps, .data$snpcov_strandedness)
    tsv_table <- dplyr::left_join(tsv_table, snpcov_table, by=c("Sample_Name"="snpcov_sample_name", "Strandedness"="snpcov_strandedness"))
  }

  ################
  ## Endogenous ##
  ################

  if (is.null(endogenous_table)) {
    ## If no endogenous table was provided, create dummy columns so output table columns are consistent
    tsv_table <- tsv_table %>%
      tibble::add_column(
        "Endogenous"=NA_real_
      )
  } else {
    ## Use endogenous post if available, else endogenous pre.
    endogenous <- endogenous_table %>%
      dplyr::mutate(
        Endogenous=dplyr::case_when(
          "endorspy_endogenous_dna_post" %in% colnames(endogenous_table) ~ endorspy_endogenous_dna_post,
          TRUE ~ endorspy_endogenous_dna
        )
      ) %>%
      dplyr::select(.data$endorspy_library_id, .data$Endogenous)
    tsv_table <- dplyr::left_join(tsv_table, endogenous, by=c("Library_ID"="endorspy_library_id"))
  }

  #############
  ## DNA DMG ##
  #############
  ## This bit returns damage at the individual level, with a weighted sum of 5p1bp % across all libs.
  if (is.null(dmg_table)) {
    ## If no damage table was provided, create dummy columns so output table columns are consistent
    tsv_table <- tsv_table %>%
      tibble::add_column(
        "Damage"=NA_real_,
        "damage_num_reads"=NA_integer_
      )
  } else {
    ## Keep only 5p 1st bp damage and number of reads (for weighted sum if needed?)
    damage <- dmg_table %>%
      dplyr::select(.data$damage_library_id, .data$damage_dmg_5p_1, .data$damage_num_reads) %>%
      dplyr::left_join(
        tsv_table %>% dplyr::select(.data$Sample_Name, .data$Library_ID),
        .,
        by=c("Library_ID"="damage_library_id")
      )
    mean_damage <- damage %>%
      dplyr::group_by(.data$Sample_Name) %>%
      dplyr::summarise(
        Damage = stats::weighted.mean(stats::na.omit(.data$damage_dmg_5p_1), stats::na.omit(.data$damage_num_reads))
      )
    tsv_table <- dplyr::left_join(tsv_table, mean_damage, by=c("Sample_Name"="Sample_Name")) %>%
      dplyr::left_join(., damage %>% dplyr::select(.data$Library_ID, .data$damage_num_reads), by="Library_ID")
  }

  ##############
  ## NUC CONT ##
  ##############
  if (is.null(nuccont_table)) {
    ## If no contamination table was provided, create dummy columns so output table columns are consistent
    tsv_table <- tsv_table %>%
      tibble::add_column(
        "Contamination_NrSnps"=NA_integer_,
        "Contamination"=NA_real_,
        "Contamination_Err"=NA_real_,
        "Contamination_Meas"=NA_character_
      )
  } else {
    contamination <- nuccont_table %>%
      dplyr::select(.data$angsd_library_id, .data$angsd_num_snps, tidyselect::contains(contamination_algorithm) & tidyselect::contains(paste0("method",contamination_method))) %>%
      dplyr::rename(
        Contamination      = tidyselect::matches(paste0("^.*_method", contamination_method,"_",contamination_algorithm,"_estimate$")),
        Contamination_Err  = tidyselect::matches(paste0("^.*_method", contamination_method,"_",contamination_algorithm,"_se$")),
        Contamination_NrSnps = .data$angsd_num_snps
      ) %>%
      dplyr::mutate(
        Contamination_Meas = "ANGSD"
      )
    tsv_table <- dplyr::left_join(tsv_table, contamination, by=c("Library_ID"="angsd_library_id"))
  }
  tsv_table
}

#' Compile library level results to sample level
#'
#' Function to calculate sample_level damage and contamination using a weighted mean across libraries
#'
#' @param x A Tibble as given by \link[eager2poseidon]{compile_eager_result_tables}.
#' @param snp_cutoff integer. The minimum number of SNPs, below which nuclear contamination results should be ignored.
#'
#' @return A tibble with the following columns replaced by weighted sums on the sample level: c(Contamination, Contamination_Error, Contamination_Meas, Damage).
#' @export
compile_across_lib_results <- function(x, snp_cutoff=100) {
  result <- x %>%
    dplyr::mutate(
      Contamination = dplyr::case_when(
        is.na(.data$Contamination_NrSnps) ~ NA_real_, ## Only comes up if TSV results do not fully match final contamination (i.e. more libs added but not processed yet)
        .data$Contamination_NrSnps < snp_cutoff ~ NA_real_,
        TRUE ~ .data$Contamination
        # TRUE ~ format(.data$x_contamination, nsmall = 3, digits = 0, trim = T)
      ),
      Contamination_Err = dplyr::case_when(
        is.na(.data$Contamination_NrSnps) ~ NA_real_, ## Only comes up if TSV results do not fully match final contamination (i.e. more libs added but not processed yet)
        .data$Contamination_NrSnps < snp_cutoff ~ NA_real_,
        TRUE ~ .data$Contamination_Err
        # TRUE ~ format(.data$x_contamination_error, nsmall = 3, digits = 0, trim = T)
      )
    ) %>%
    dplyr::group_by(.data$Sample_Name) %>%
    dplyr::summarise(
      .groups = 'keep',
      ## Contamination columns are weighted mean from libraries.
      contamination_weighted_mean = library_unique_weighted_mean(.data$Library_ID, .data$Contamination, .data$Contamination_NrSnps, na.rm = T),
      contamination_weigthed_mean_err = library_unique_weighted_mean(.data$Library_ID, .data$Contamination_Err, .data$Contamination_NrSnps, na.rm = T),
      sample_Contamination = dplyr::if_else(
        is.nan(.data$contamination_weighted_mean),
        #TRUE
        NA_character_,
        #FALSE
        .data$contamination_weighted_mean %>% format(., nsmall = 3, digits = 1, trim = T) ## Change to type 'char' and format to three decimals
      ),
      sample_Contamination_Err = dplyr::if_else(
        is.nan(.data$contamination_weigthed_mean_err),
        #TRUE
        NA_character_,
        #FALSE
        .data$contamination_weigthed_mean_err %>% format(., nsmall = 3, digits = 1, trim = T) ## Change to type 'char' and format to three decimals
      ),
      sample_Contamination_Note = dplyr::if_else(
        is.na(.data$sample_Contamination),
        #TRUE
        paste0("No results found, or no libraries exceeded cutoff of ", snp_cutoff," SNPs after depth filtering."),
        #FALSE
        create_contamination_note(.data$Library_ID, .data$Contamination_NrSnps, snp_cutoff)
      ),
      sample_Contamination_Meas = dplyr::if_else(
        is.na(.data$sample_Contamination),
        #TRUE
        NA_character_,
        #FALSE
        "ANGSD"
      ),
      ## Damage is weighted mean of damage from each library
      mean_damage = library_unique_weighted_mean(.data$Library_ID, .data$Damage, .data$damage_num_reads, na.rm=T),
      sample_Damage = dplyr::if_else(
        is.nan(.data$mean_damage),
        ## TRUE
        NA_character_,
        ## FALSE
        .data$mean_damage %>% format(., nsmall = 3, digits = 1, trim = T) ## Change to type 'char' and format to three decimals
      ),
      ## Endogenous is max value across all libraries.
      sample_Endogenous = dplyr::if_else(
        ## If all values are NA, then return NA (max complains is only NAs are given). Otherwise return max after removing NAs.
        sum(is.na(.data$Endogenous)) == length(.data$Endogenous),
        ## TRUE
        NA_real_,
        ## FALSE
        ## For some reason, if all values are NA, this still prints out an error that it will return -Inf,
        ##    but the endogenous value is correctly inferred.
        max(.data$Endogenous, na.rm=T)
        ),
      ## Keep track of libraries in the results
      sample_Library_Names=vctrs::vec_unique(.data$Library_ID) %>% paste0(., collapse=";")
    ) %>%
    ## Keep sample level columns, and convert to their poseidon names
    dplyr::select( tidyselect::starts_with("sample_") ) %>%
    dplyr::rename_with( ~sub('^sample_', '', .)) %>%
    ## Keep only distinct rows to avoid row duplications
    dplyr::distinct() %>%
    dplyr::left_join(x, ., by="Sample_Name", suffix = c(".x", ".y")) %>%
    ## Remove duplicate columns from input data, and rename new columns to poseidon names again
    dplyr::select(-tidyselect::ends_with(".x")) %>%
    dplyr::rename_with(.fn=~sub('\\.y$', '', .), .cols=tidyselect::ends_with(".y"))
  result
}

#' Create contamination note from eager results for poseidon
#'
#' Create contamination note for the sample, ignoring duplicated library entries.
#'
#' @param lib_ids character. The names of the libraries the contamination estimates correspond to.
#' @param nr_snps character. The number of SNPs for the contamination estimate.
#' @inheritParams compile_across_lib_results
#'
#' @return character. An appropriate Contamination_Note
#' @export
create_contamination_note <- function(lib_ids, nr_snps, snp_cutoff=100) {
  x <- tibble::tibble(Library_ID=lib_ids, Contamination_NrSnps=nr_snps) %>%
    dplyr::distinct()
  paste0("Nr Snps (per library): ", paste(x$Contamination_NrSnps, collapse = ";"), ". Estimate and error are weighted means of values per library. Libraries with fewer than ", snp_cutoff, " were excluded.")
}

#' Calculate the weighted mean at the library level
#'
#' This function calculates the weighted mean of the measures provided, given the weights provided.
#' Should be provided the measures for one Sample_Name/Poseidon_ID at a time.
#' Unlike stats::weighted.mean, this function first filters for distinct rows, to avoid counting duplicated results
#' as independent observations.
#'
#' @param lib_ids character. The names of the libraries the measures correspond to.
#' @param measure character. The measures to get a weighted mean for
#' @param weight character. The weight of each measure.
#' @param na.rm logical. Should NAs be removed from the weighted mean. passed on to stats::weighted.mean(na.rm)
#'
#' @return character. The weighted mean of the unique library-level measures.
library_unique_weighted_mean <- function(lib_ids, measure, weight, na.rm = T) {
  x <- tibble::tibble(Library_ID=lib_ids, Measure=measure, Weight=weight) %>%
    dplyr::distinct()

  result <- stats::weighted.mean(x$Measure, x$Weight, na.rm = na.rm)

  result
}
