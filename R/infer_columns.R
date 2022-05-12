#' Infer Poseidon encoded genetic sex from relative X/Y coverage
#'
#' @param x_rate double. The relative coverage on the X chromosome.
#' @param y_rate double. The relative coverage on the Y chromosome.
#' @param XX_cutoffs numeric. A vector containing the minimum x_rate, maximum x_rate, minimum y_rate and maximum y_rate for XX calls.
#' @param XY_cutoffs numeric. A vector containing the minimum x_rate, maximum x_rate, minimum y_rate and maximum y_rate for XY calls.
#'
#' @return character. The inferred poseidon-encoded genetic sex.
#' @export
#'
infer_genetic_sex <- function(x_rate, y_rate, XX_cutoffs = c(0.7, 1.2, 0.0, 0.1), XY_cutoffs = c(0.2, 0.6, 0.3, 0.6)) {

  if (is.na(x_rate) || is.na(y_rate)) {
    gsex <- "U"
    return (gsex)
  }

  if ( x_rate >= XX_cutoffs[1] && x_rate <= XX_cutoffs[2] && y_rate >= XX_cutoffs[3] && y_rate <= XX_cutoffs[4]) {
    gsex <- "F"
  } else if ( x_rate >= XY_cutoffs[1] && x_rate <= XY_cutoffs[2] && y_rate >= XY_cutoffs[3] && y_rate <= XY_cutoffs[4]) {
    gsex <- "M"
  } else {
    gsex <- "U"
  }
  return (gsex)
}
