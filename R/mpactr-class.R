mpactr <- R6Class("mpactr", public = list(
  # Properties

  # Constructor
  initialize = function(peak_table_path, meta_data_path) {
    private$peak_table = data.table(readr::read_csv(peak_table_path, skip = 2, show_col_types = FALSE))
    private$meta_data = data.table(readr::read_csv(meta_data_path, show_col_types = FALSE))
    stopifnot(any(class(private$peak_table) == "data.table"))
    stopifnot(any(class(private$meta_data) == "data.table"))
  }
),
  private = list(
    peak_table = NA,
    meta_data = NA
))
