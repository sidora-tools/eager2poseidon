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
#' @inheritParams infer_genetic_sex
#'
#' @return A tibble with the compiled information from the provided eager output tables.
#' @export
#'
compile_eager_result_tables <- function(tsv_table=NULL, sexdet_table=NULL, snpcov_table=NULL, dmg_table=NULL, endogenous_table=NULL, nuccont_table=NULL, contamination_method="1", contamination_algorithm = "ml", XX_cutoffs = c(0.7, 1.2, 0.0, 0.1), XY_cutoffs = c(0.2, 0.6, 0.3, 0.6)) {

  ## Without a TSV info table, nothing can be compiled
  if (is.null(tsv_table)) {stop("[compile_eager_result_tables()]: tsv_table is required.")}

  ## Validate contamination input params
  if (! contamination_method %in% c("1","2")) { stop(paste0("[compile_eager_result_tables()]: Invalid preferred nuclear contamination estimation method provided.\n\tValid entries are c('1','2').\n\tYou provided: ",contamination_method))}

  if (! contamination_algorithm %in% c("ml","mom")) { stop(paste0("[compile_eager_result_tables()]: Invalid preferred nuclear contamination estimation method provided.\n\tValid entries are c('ml','mom').\n\tYou provided: ", contamination_algorithm))}

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
    sexdet <- sexdet_table %>% dplyr::mutate(
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
      dplyr::select(.data$sexdet_input_bam, .data$Genetic_Sex, .data$Sex_Determination_Note)
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
