#' Collect poseidon janno data from pandora
#'
#' @param janno A tibble Path to the source .janno file. Passed to \link[poseidonR:janno]{read_janno}
#' @param credentials character. Path to a credentials file containing four lines listing the database host,
#' the port of the database server, user and password, respectively.
#' Passed to \link[sidora.core:get_pandora_connection]{get_pandora_connection}
#' @param trust_uncalibrated_dates logical. Should any uncalibrated dates in pandora be trusted?
#' If set to TRUE, then \link[poseidonR:quickcalibrate]{quickcalibrate()} is used to calibrate these dates on the fly.
#'
#' @return A tibble containing the poseidon janno fields of
#' @export
#'
#' @importFrom magrittr "%>%"

import_pandora_data <- function(janno, credentials, trust_uncalibrated_dates) {
  ## Read janno
  input_janno <- poseidonR::read_janno(janno, to_janno = F, validate = TRUE)
  sample_ids <- input_janno %>%
    dplyr::select(tidyselect::any_of(c('Individual_ID', 'Poseidon_ID')))
  ## Rename to Individual_ID to Poseidon_ID if the janno is from an older Poseidon version. Throw warning.
  if (names(sample_ids) == "Individual_ID" ) {
    warning("Your janno file appears to be from an older poseidon version. It will be updated to poseidon version '2.5.0'.")
    ## Change column name in imported janno and sample names
    input_janno <- input_janno %>% dplyr::rename("Poseidon_ID"="Individual_ID")
    names(sample_ids) <- c("Poseidon_ID")
  }

  write("Querying Pandora for individual data.", file=stderr())
  pandora_table <- get_individual_pandora_data(sample_ids, credentials)
  pandora_table <- pandora_table %>%
    dplyr::mutate(
      Poseidon_ID=.data$individual.Full_Individual_Id,
      Collection_ID=.data$individual.Archaeological_ID,
      Country=.data$site.Country,
      Location=dplyr::if_else(.data$site.Locality == "", "n/a", .data$site.Locality),
      Longitude=.data$site.Longitude,
      Latitude=.data$site.Latitude,
      Date_C14_Labnr=dplyr::na_if(.data$individual.C14_Id,""),
      Date_BC_AD_Start_pandora=dplyr::case_when(
        ## If no C14 ID is given in pandora, don't trust the Calibrated date field, else take it as is.
        .data$Date_C14_Labnr %in% c("", NA) ~ NA_integer_,
        TRUE ~ .data$individual.C14_Calibrated_From
      ),
      Date_BC_AD_Stop_pandora=dplyr::case_when(
        ## If no C14 ID is given in pandora, don't trust the Calibrated date field, else take it as is.
        .data$Date_C14_Labnr %in% c("", NA) ~ NA_integer_,
        TRUE ~ .data$individual.C14_Calibrated_To
      ),
      Date_Type_pandora=dplyr::case_when(
        ## If a C14 Lab number is provided, set date type to C14.
        .data$Date_C14_Labnr %in% c("", NA) ~ NA_character_,
        TRUE ~ "C14"
      ),
      Date_Note_pandora=dplyr::case_when(
        ## If a C14 Lab number is provided, propagate the date info field.
        .data$Date_C14_Labnr %in% c("", NA) ~ NA_character_,
        TRUE ~ .data$individual.C14_Info
      ),
      ## Initialise These as empty for calibration to potentially fill in.
      ##    Need to exist for dplyr::coalesce() to not complain about columns not existing.
      Date_BC_AD_Start = NA_integer_,
      Date_BC_AD_Stop = NA_integer_,
      Date_BC_AD_Median = NA_integer_,
      Date_C14_Uncal_BP = NA_integer_,
      Date_C14_Uncal_BP_Err = NA_integer_,
      Date_Type_quickcal=NA_character_,
      Date_Note_quickcal=NA_character_
    )
  ## If user trusts uncalibrated dates (and C14 IDs are set up correctly) for their samples,
  ##    then use that data to quickcalibrate (calibrated dates are retained if existing)
  if (trust_uncalibrated_dates) {
    pandora_table <- pandora_table %>%
      dplyr::mutate(
        ## If no Lab Nr is given in pandora, the C14 data is not considered a valid entry.
        Date_C14_Uncal_BP=dplyr::if_else(.data$Date_C14_Labnr %in% c("", NA), NA_integer_, .data$individual.C14_Uncalibrated),
        Date_C14_Uncal_BP_Err=dplyr::if_else(.data$Date_C14_Labnr %in% c("", NA), NA_integer_, .data$individual.C14_Uncalibrated_Variation),
        ## Only use calibration values if calibrated values are not in pandora
        poseidonR::quickcalibrate(
          .data$Date_C14_Uncal_BP,
          .data$Date_C14_Uncal_BP_Err,
          select_calibration_curve(.data$site.Latitude),
          allowOutside=T
        ),
        Date_Type_quickcal = dplyr::if_else(.data$Date_C14_Labnr %in% c("", NA), NA_character_, "C14"),
        Date_Note_quickcal = dplyr::if_else(.data$Date_C14_Labnr %in% c("", NA), NA_character_, "quickcalibration from Date_C14_Uncal_BP and Date_C14_Uncal_BP_Err")
      )
  }
  pandora_table <- pandora_table %>%
    dplyr::mutate(
      ## Prioritise pandora's calibrated dates when available.
      Date_BC_AD_Start=dplyr::coalesce(.data$Date_BC_AD_Start_pandora,.data$Date_BC_AD_Start),
      Date_BC_AD_Stop=dplyr::coalesce(.data$Date_BC_AD_Stop_pandora,.data$Date_BC_AD_Stop),
      Date_Type=dplyr::coalesce(.data$Date_Type_pandora,.data$Date_Type_quickcal),
      Date_Note=dplyr::coalesce(.data$Date_Note_pandora,.data$Date_Note_quickcal)
      ) %>%
    dplyr::select(
      .data$Poseidon_ID,
      .data$Collection_ID,
      .data$Country,
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
      .data$Date_Note
      )

  write("Parsing of Pandora table completed.", file=stderr())

  pandora_table
}

#' Pull pandora information for specific individual IDs
#'
#' Query pandora for information across tabs Site and Sequencing, for a provided set of individual IDs.
#'
#' @inheritParams import_pandora_data
#' @param sample_ids character. A vector of individual IDs to pull pandora data for.
#'
#' @return A tibble containing the pandora information for individuals present in the janno file, from Site to Sequencing.
#' @export
#'
#' @importFrom magrittr "%>%"
get_individual_pandora_data <- function(sample_ids, credentials) {

  ## Pull pandora data
  con <- sidora.core::get_pandora_connection(credentials)

  ## Get complete pandora table, then filter for samples of interest
  pandora_table <- sidora.core::join_pandora_tables(
    sidora.core::get_df_list(
      c(sidora.core::make_complete_table_list(
        c("TAB_Site", "TAB_Sequencing")
      )), con = con
    )
  ) %>% sidora.core::convert_all_ids_to_values(con = con)
  DBI::dbDisconnect(con)
  pandora_table <- pandora_table %>% dplyr::filter(.data$individual.Full_Individual_Id %in% sample_ids$Poseidon_ID)

  write("Information successfully pulled from Pandora.", file=stderr())

  return(pandora_table)
}

#' Infer the correct calibration curve for quickcalibrate, based on the Pandora Latitude
#'
#' Pick a calibration curve for
#'
#' @param Latitude numeric. Latitude in decimal degrees.
#'
#' @return A string with the preferred calibration curve (intcal20 for northern hemisphere and shcal20 for the southern hemisphere).
#' @export
#'
#' @examples
#' select_calibration_curve(-50)  ## Northern hemisphere
#' select_calibration_curve(50)   ## Southern Hemisphere
#' select_calibration_curve(0)    ## 0 defaults to Northern hemisphere
#'
#' \dontrun{
#' select_calibration_curve("banana")
#' }
select_calibration_curve <- function(Latitude) {
  ## Bchron (the backend for quickcalibrate) does not recognise NAs in the calibration curves.
  ##   Therefor, if a site has no latitude set, or it is wrongly set, stop and complain.
  if (!is.numeric(Latitude) || is.na(Latitude)) { stop("Latitude value should be of class 'numeric', and NOT NA.") }
  dplyr::case_when(
    Latitude >= 0 ~ "intcal20",
    Latitude < 0 ~ "shcal20",
    TRUE ~ "none"
  )
}
