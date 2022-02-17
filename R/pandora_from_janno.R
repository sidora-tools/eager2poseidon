#' pandora_from_janno
#'
#' Query pandora for information for all individuals present in a janno file.
#'
#' @param janno character. Path to the source .janno file. Passed to \link[poseidonR:janno]{read_janno}
#' @param credentials character. Path to a credentials file containing four lines listing the database host, the port of the database server, user and password, respectively. Passed to \link[sidora.core:get_pandora_connection]{get_pandora_connection}
#'
#' @return A tibble containing the pandora information for individuals present in the janno file, from Site to Sequencing.
#' @export
#'
#' @importFrom magrittr "%>%"
pandora_from_janno <- function(janno, credentials) {
  ## Read janno
  input_janno <- poseidonR::read_janno(janno, to_janno = TRUE, validate = TRUE)
  sample_ids <- poseidonR::read_janno(janno, to_janno = TRUE, validate = TRUE) %>%
                  dplyr::select(tidyselect::any_of(c('Individual_ID', 'Poseidon_ID')))
  if (names(sample_ids) == "Individual_ID" ) {
    warning("Your janno file appears to be from an older poseidon version. It will be updated to poseidon version '2.5.0'.")
    }
  names(sample_ids) <- c("Poseidon_ID") ## Rename to Individual_ID to Poseidon_ID if the janno is from an older Poseidon version

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
  pandora_table <- pandora_table %>% dplyr::filter(individual.Full_Individual_Id %in% sample_ids$Poseidon_ID)

  return(pandora_table)
}
#
