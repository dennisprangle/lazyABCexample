#' Calculate ESS
#'
#' Calculates effective sample size from a vector of weights
#'
#' @param w Vector of weights
#'
#' @return Effective sample size
#' @export
ess <- function(w) {
    length(w) * mean(w)^2 / mean(w^2)
}
