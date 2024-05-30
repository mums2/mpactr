rsd <- function(values) {
  ifelse(mean(values) != 0, (sd(values) / mean(values)), NA_real_)
}