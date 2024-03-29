% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/janno_operations.R
\name{standardise_janno}
\alias{standardise_janno}
\title{Standardise janno file}
\usage{
standardise_janno(janno_fn)
}
\arguments{
\item{janno_fn}{A tibble Path to the source .janno file. Passed to \link[poseidonR:janno]{read_janno}}
}
\value{
A list containing the standardised janno tibble and a tibble with the sample Ids of the janno table
}
\description{
Read a janno file into a tibble with standardised data types, and update the column names to poseidon v2.5.0 if necessary
}
