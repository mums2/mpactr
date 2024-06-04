mpactr <- R6Class("mpactr", public = list(
  # Properties
  peak_table_csv = NA,
  meta_data_csv = NA,
  # sample_list_csv = NA,
  # Constructor
  initialize = function(peak_table_path, meta_data_path){
    self$peak_table_csv = data.table(readr::read_csv(peak_table_path, skip =2))
    self$meta_data_csv = data.table(readr::read_csv(meta_data_path))
    stopifnot(any(class(self$peak_table_csv) == "data.table"))
    stopifnot(any(class(self$meta_data_csv) == "data.table"))
    #self$sample_list_csv <- readr::read_csv(here::here(sample_list_path))
  },
  setup = function() {
    self$initialize_data() # TODO Redo functions to deal with references
  },
  filter_data = function(filter_function, params = ...) {
    summary <- filter_function(self$peak_table_csv, params)
    self$logger(summary)
  },
  logger = function(summary){
    print(summary)
    print("Some other random statistics")
  }
),
private = list(
   set_kmd <- function() {
    self$peak_table_csv[ , kmd := mz - floor(mz)]
  },
  initialize_data <- function() {
    self$peak_table_csv[, .(mz = `m/z`, rt = `Retention time (min)`)]
    # data_frame$Compound <- as.numeric(data_frame$Compound)
    self$peak_table_csv <- self$peak_table_csv[which( rowSums(self$peak_table_csv[, !(colnames(self$peak_table_csv) %in% c("Compound", "mz", "rt")) ] ) > 0), ] # need to fix to work with data.table
    # calc kmd
    self$set_kmd()
  }
)
)
