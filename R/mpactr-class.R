mpactr <- R6Class("mpactr", public = list(
  # Properties

  # Constructor
  initialize = function(peak_table_path, meta_data_path) {
    if((any(class(peak_table_path) %in% c("data.table", "data.frame"))) &&
    (any(class(meta_data_path) %in% c("data.table", "data.frame"))))
    {
      private$peak_table = data.table(peak_table_path)
      private$meta_data = data.table(meta_data_path)
    }
    else{
      private$peak_table = data.table(readr::read_csv(peak_table_path, skip = 2, show_col_types = FALSE))
      private$meta_data = data.table(readr::read_csv(meta_data_path, show_col_types = FALSE))
    }
    stopifnot(any(class(private$peak_table) == "data.table"))
    stopifnot(any(class(private$meta_data) == "data.table"))
  },
  isMultipleTechReps = function() {
    
    any(private$meta_data[ , .N, by = Sample_Code][["N"]])
  }
),
  private = list(
    peak_table = NA,
    meta_data = NA,
    raw_peak_table = NA
))
