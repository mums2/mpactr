filter_pactr <- R6Class("filter_pactr", public = list(
  mpactr_data = NA,
  logger = NA, # Might have to add getter/setters
  initialize = function(mpactr) {
    self$mpactr_data <- mpactr
    self$logger <- new.env(hash = TRUE)
    self$logger[["list_of_summaries"]] <- list()
  }
  )
)

#' @export
print.filter_pactr <- function(mpactr_object)
{
  print(mpactr_object$mpactr_data$get_peak_table())
}
