#' @title Limit Cores
#'
#' @description
#' Limits the amount of cores used to two. Mostly used in test
#' and examples so we can pass cran checks.
#' @return does not return anything.
#' @export
#' @examples
#' limit_cores()
#'
limit_cores <- function() {
  Sys.setenv("OMP_THREAD_LIMIT" = 1)
  setDTthreads(threads = 1)
}
