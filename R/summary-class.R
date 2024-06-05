# filter
# failed_ions
# passing_ions

summary <- R6Class("summary", public = list(
  filter  = NA,
  failed_ions = NA,
  passed_ions = NA,
  initialize = function(filter, failed_ions, passed_ions)
  {
    stopifnot(any(class(filter) == "character"))
    stopifnot(any(class(failed_ions) == "numeric"))
    stopifnot(any(class(passed_ions) == "numeric"))
    self$filter = filter
    self$failed_ions = failed_ions
    self$passed_ions = passed_ions

  }
))