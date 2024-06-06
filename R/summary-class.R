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
    stopifnot(class(failed_ions) == "numeric")
    stopifnot(class(passed_ions) == "numeric")
    self$filter = filter
    self$failed_ions = failed_ions
    self$passed_ions = passed_ions
  },
  summarize = function(x)
  {
    print(paste0("For ", self$filter, ": ", length(self$failed_ions), " ions failed insource ion filter. ",
                              length(self$passed_ions), " ions remain."))
  }
))

