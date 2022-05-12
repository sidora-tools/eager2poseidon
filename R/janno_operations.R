#' Fill in a janno table
#'
#' Fill in missing information in a janno table with relevant data from pandora,
#' an eager TSV and the MultiQC general stats table of an eager run using that TSV.
#' Any cells that are already filled in the janno table will be kept as is.
#'
#' @param input_janno_table A standardised janno tibble, as returned by \link[eager2poseidon:standardise_janno]{standardise_janno}.
#' @param external_results_table A tibble with the collected results form pandora and eager, as returned by \link[eager2poseidon:collate_external_results]{collate_external_results}.
#' @param genotype_ploidy character. The genotype ploidy of the genotypes in the dataset. Can be either 'diploid' or 'haploid'.
#'
#' @return A tibble with missing values filled in for the fields
#' Collection_ID, Country, Site, Location, Longitude,
#' Latitude, Date_C14_Labnr, Date_BC_AD_Start, Date_BC_AD_Stop, Date_BC_AD_Median,
#' Date_C14_Uncal_BP, Date_C14_Uncal_BP_Err, Date_Type, Nr_Libraries, Capture_Type,
#' UDG, Library_Built, Damage, Nr_SNPs, Endogenous, Contamination, Contamination_Err
#' @export
fill_in_janno <- function(input_janno_table, external_results_table, genotype_ploidy) {

  ## Validate genotype ploidy option input
  valid_ploidy_entries <- c("haploid", "diploid")
  if (!genotype_ploidy %in% valid_ploidy_entries) {
    stop(call. = F, "\nInvalid genotype ploidy: '", genotype_ploidy, "'\nAccepted values: ", paste(valid_ploidy_entries, collapse = ", "))
  }

  ## Set individual order so it matches the input janno
  ind_order <- dplyr::pull(input_janno_table, .data$Poseidon_ID)
  ## If no note field is present add one full of NAs
  if ( ! "Note" %in% names(input_janno_table)) {
    janno_table <- input_janno_table %>% dplyr::mutate(Note=NA_character_)
  } else {
    janno_table <- input_janno_table
  }

  ## Separate genetic sex from the rest of the external data, since it needs special handling.
  sexdet_table <- external_results_table %>% dplyr::select(.data$Poseidon_ID, .data$Genetic_Sex)
  external_results_table <- external_results_table %>% dplyr::select(-.data$Genetic_Sex)

  output_janno <- dplyr::full_join(janno_table, external_results_table, by = "Poseidon_ID") %>%
    dplyr::mutate(
      Collection_ID = dplyr::coalesce(.data$Collection_ID.x, .data$Collection_ID.y),
      Country = dplyr::coalesce(.data$Country.x, .data$Country.y),
      Site = dplyr::coalesce(.data$Site.x, .data$Site.y),
      Location = dplyr::coalesce(.data$Location.x, .data$Location.y),
      Longitude = dplyr::coalesce(.data$Longitude.x, .data$Longitude.y),
      Latitude = dplyr::coalesce(.data$Latitude.x, .data$Latitude.y),
      Date_C14_Labnr = dplyr::coalesce(.data$Date_C14_Labnr.x, .data$Date_C14_Labnr.y),
      Date_BC_AD_Start = dplyr::coalesce(.data$Date_BC_AD_Start.x, .data$Date_BC_AD_Start.y),
      Date_BC_AD_Stop = dplyr::coalesce(.data$Date_BC_AD_Stop.x, .data$Date_BC_AD_Stop.y),
      Date_BC_AD_Median = dplyr::coalesce(.data$Date_BC_AD_Median.x, .data$Date_BC_AD_Median.y),
      Date_C14_Uncal_BP = dplyr::coalesce(.data$Date_C14_Uncal_BP.x, .data$Date_C14_Uncal_BP.y),
      Date_C14_Uncal_BP_Err = dplyr::coalesce(.data$Date_C14_Uncal_BP_Err.x, .data$Date_C14_Uncal_BP_Err.y),
      Date_Type = dplyr::coalesce(.data$Date_Type.x, .data$Date_Type.y),
      Nr_Libraries = dplyr::coalesce(.data$Nr_Libraries.x, .data$Nr_Libraries.y),
      Capture_Type = dplyr::coalesce(.data$Capture_Type.x, .data$Capture_Type.y),
      UDG = dplyr::coalesce(.data$UDG.x, .data$UDG.y),
      Library_Built = dplyr::coalesce(.data$Library_Built.x, .data$Library_Built.y),
      Damage = dplyr::coalesce(.data$Damage.x, .data$Damage.y),
      Nr_SNPs = dplyr::coalesce(.data$Nr_SNPs.x, .data$Nr_SNPs.y),
      Endogenous = dplyr::coalesce(.data$Endogenous.x, .data$Endogenous.y),
      Contamination = dplyr::coalesce(.data$Contamination.x, .data$Contamination.y),
      Contamination_Err = dplyr::coalesce(.data$Contamination_Err.x, .data$Contamination_Err.y),
      Contamination_Meas = dplyr::coalesce(.data$Contamination_Meas.x, .data$Contamination_Meas.y),
      Contamination_Note = dplyr::coalesce(.data$Contamination_Note.x, .data$Contamination_Note.y),
      Genotype_Ploidy = dplyr::coalesce(.data$Genotype_Ploidy, genotype_ploidy),
      Note = dplyr::coalesce(.data$Note.x, .data$Note.y)
    ) %>%
    dplyr::select(
      .data$Poseidon_ID,
      .data$Collection_ID,
      .data$Country,
      .data$Site,
      .data$Location,
      .data$Longitude,
      .data$Latitude,
      .data$Date_C14_Labnr,
      .data$Date_BC_AD_Start,
      .data$Date_BC_AD_Stop,
      .data$Date_BC_AD_Median,
      .data$Date_C14_Uncal_BP,
      .data$Date_C14_Uncal_BP_Err,
      .data$Date_Type,
      .data$Nr_Libraries,
      .data$Capture_Type,
      .data$UDG,
      .data$Library_Built,
      .data$Damage,
      .data$Nr_SNPs,
      .data$Endogenous,
      .data$Contamination,
      .data$Contamination_Err,
      .data$Contamination_Meas,
      .data$Contamination_Note,
      tidyselect::any_of(
        c(
          "Alternative_IDs",
          "Group_Name",
          tidyselect::starts_with("Relation_*"),
          "Genetic_Sex",
          "MT_Haplogroup",
          "Y_Haplogroup",
          "Source_Tissue",
          "Genotype_Ploidy",
          "Coverage_on_Target_SNPs",
          "Genetic_Source_Accession_IDs",
          "Primary_Contact",
          "Publication",
          "Note"
        )
      )
    ) %>%
    fill_genetic_sex(., sexdet_table) %>%
    ## Set the order so it matches the input janno
    dplyr::mutate(
      Poseidon_ID = factor(.data$Poseidon_ID, levels = ind_order),
      Data_Preparation_Pipeline_URL = "https://nf-co.re/eager",
    ) %>%
    dplyr::arrange(.data$Poseidon_ID)

  return(output_janno)
}

#' Standardise janno file
#'
#' Read a janno file into a tibble with standardised data types, and update the column names to poseidon v2.5.0 if necessary
#'
#' @param janno_fn A tibble Path to the source .janno file. Passed to \link[poseidonR:janno]{read_janno}
#'
#' @return A list containing the standardised janno tibble and a tibble with the sample Ids of the janno table
#' @export
standardise_janno <- function(janno_fn) {

  ## Read janno and normalise column names for Individual/Poseidon IDs
  input_janno <- poseidonR::read_janno(janno_fn, to_janno = F, validate = F) %>%
    dplyr::mutate(
      dplyr::across(
        .cols = tidyselect::any_of(c(
          "Date_C14_Uncal_BP",
          "Date_C14_Uncal_BP_Err",
          "Date_BC_AD_Median",
          "Date_BC_AD_Start",
          "Date_BC_AD_Stop",
          "No_of_Libraries",
          "Nr_autosomal_SNPs",
          "Nr_Libraries",
          "Nr_SNPs"
        )),
        .fn = as.integer
      ),
      dplyr::across(
        .cols = tidyselect::any_of(c(
          "Latitude",
          "Longitude",
          "Coverage_1240K",
          "Endogenous_DNA",
          "Damage",
          "Endogenous",
          "Coverage_on_Target_SNPs"
        )),
        .fn = as.numeric
      )
    )

  sample_ids <- input_janno %>%
    dplyr::select(tidyselect::any_of(c("Individual_ID", "Poseidon_ID")))
  ## Rename to Individual_ID to Poseidon_ID if the janno is from an older Poseidon version. Throw warning.
  if (names(sample_ids) == "Individual_ID") {
    warning("Your janno file appears to be from an older poseidon version. It will be updated to poseidon version '2.5.0'.\n\t\tEnsure you UPDATE POSEIDON VERSION on your POSEIDON.yaml as well!")
    ## Change column names in imported janno and in sample_ids to poseidon 2.5.0 versions
    input_janno <- input_janno %>% dplyr::rename(
      "Poseidon_ID" = "Individual_ID",
      "Nr_Libraries" = "No_of_Libraries",
      "Capture_Type" = "Data_Type",
      "Nr_SNPs" = "Nr_autosomal_SNPs",
      # "Endogenous" = "Endogenous_DNA",
      "Contamination" = "Xcontam",
      "Contamination_Err" = "Xcontam_stderr",
      "Coverage_on_Target_SNPs" = "Coverage_1240K"
    ) %>%
      ## Add missing columns that we then coalesce onto.
      dplyr::mutate(
        Contamination_Meas = NA_character_,
        Contamination_Note = NA_character_,
      )
  }

  return(input_janno)
}

#' Collect Pandora and Eager results
#'
#' @inheritParams import_pandora_data
#' @inheritParams import_eager_results
#'
#' @return A tibble with the collected results from Pandora and eager.
#' @export
#'
collate_external_results <- function(sample_ids, eager_tsv_fn, general_stats_fn, credentials, keep_only = "none", trust_uncalibrated_dates = F, snp_cutoff) {
  pandora_table <- import_pandora_data(sample_ids, credentials, trust_uncalibrated_dates)
  eager_table <- import_eager_results(eager_tsv_fn, general_stats_fn, keep_only, snp_cutoff)

  external_results <- dplyr::full_join(pandora_table, eager_table, by = c("Poseidon_ID" = "Sample_Name"))

  return(external_results)
}

#' Add genetic sex information to the input janno, throwing a warning if any changes happen, so users know to update their .fam/.ind files
#'
#' @param input_janno_table tibble. A standardised janno tibble, as returned by \link[eager2poseidon:standardise_janno]{standardise_janno}.
#' @param genetic_sex_table tibble. A tibble containing the Poseidon_Id and Sex of inferred genetic sex of each sample
#'
#' @return tibble. A tibble with filled in Genetic_Sex information.
#' A warning is printed if any entries have been updated, requesting manual changes to be made to the .fam/.ind files of the poseidon package.
#'
fill_genetic_sex <- function(input_janno_table, genetic_sex_table) {

  input <- input_janno_table %>%
    dplyr::mutate(
      old_gsex = .data$Genetic_Sex,
      Genetic_Sex = dplyr::na_if(.data$Genetic_Sex, "U")
      )

  output_janno <- dplyr::full_join(input, genetic_sex_table, by = "Poseidon_ID") %>%
    dplyr::mutate(
      Genetic_Sex = dplyr::coalesce(.data$Genetic_Sex.x, .data$Genetic_Sex.y)
    ) %>%
    dplyr::select(-.data$Genetic_Sex.x, -.data$Genetic_Sex.y, -.data$old_gsex)

  if ( ! identical(input$old_gsex, output_janno$Genetic_Sex)) {
    warning("
  The Genetic_Sex field has been updated for a number of individuals, based on sexdeterrmine results found in the eager output tables.
            Please update the .fam/.ind file of the package to reflect this change!", call=F)

  }
  return(output_janno)
}
