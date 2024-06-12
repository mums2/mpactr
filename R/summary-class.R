summary <- R6::R6Class("summary", public = list(
  initialize = function(filter, failed_ions, passed_ions)
  {
    stopifnot(any(class(filter) == "character"))
    stopifnot(class(failed_ions) == "numeric")
    stopifnot(class(passed_ions) == "numeric")
    private$filter = filter
    private$failed_ions = failed_ions
    private$passed_ions = passed_ions
  },
  summarize = function(x)
  {
    print(paste0("For ", private$filter, ": ", length(private$failed_ions), " ions failed insource ion filter. ",
                 length(private$passed_ions), " ions remain."))
  },
  get_failed_ions = function()
  {
    return(private$failed_ions)
  },
  get_passed_ions = function()
  {
    return(private$passed_ions)
  },
  get_filter = function()
  {
    return(private$filter)
  }
),
 private = list(
   filter = NA,
   failed_ions = NA,
   passed_ions = NA
 )
)

