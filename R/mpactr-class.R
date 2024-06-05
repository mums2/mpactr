mpactr <- R6Class("mpactr", public = list(
  # Properties
  peak_table = NA,
  meta_data = NA,
  # sample_list_csv = NA,
  # Constructor
  initialize = function(peak_table_path, meta_data_path) {
    self$peak_table = data.table(readr::read_csv(peak_table_path, skip = 2))
    self$meta_data = data.table(readr::read_csv(meta_data_path))
    stopifnot(any(class(self$peak_table) == "data.table"))
    stopifnot(any(class(self$meta_data) == "data.table"))
    #self$sample_list_csv <- readr::read_csv(here::here(sample_list_path))
  }
))
