#' Relative Standard Deviation
#'
#' @param values a `numeric` vector of values.
#' @noRd
rsd <- function(values) {
  ifelse(mean(values) != 0, (sd(values) / mean(values)), NA_real_)
}

permutations <- function(x) {
  permutations <- vector("list", length(x))
  for(i in seq_len(length(x))) {
    permutations[[i]] <- x[-i]
  }
  permutations
}
