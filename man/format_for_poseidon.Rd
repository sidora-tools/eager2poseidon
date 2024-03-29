% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_eager_results.R
\name{format_for_poseidon}
\alias{format_for_poseidon}
\title{Internal function to convert information from eager lingo to poseidon accepted strings}
\usage{
format_for_poseidon(x, field)
}
\arguments{
\item{x}{character. A string to format from eager language to poseidon language.}

\item{field}{character. The field to format for. One of 'UDG', 'Library_Built', or 'Capture_Type'.}
}
\value{
character. A valid Poseidon Capture_Type, Library_Built, or UDG value.
}
\description{
Internal function to convert information from eager lingo to poseidon accepted strings
}
