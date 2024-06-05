mpactr <- R6Class("mpactr", public = list(
  # Properties
  peak_table = NA,
  meta_data = NA,
  # sample_list_csv = NA,
  # Constructor
  initialize = function(peak_table_path, meta_data_path) {
    self$peak_table = data.table(readr::read_csv(peak_table_path, skip =2))
    self$meta_data = data.table(readr::read_csv(meta_data_path))
    stopifnot(any(class(self$peak_table) == "data.table"))
    stopifnot(any(class(self$meta_data) == "data.table"))
    #self$sample_list_csv <- readr::read_csv(here::here(sample_list_path))
  },
  setup = function() {
    private$initialize_data() # TODO Redo functions to deal with references
  },
  logger = function(summary){
    print(summary)
    print("Some other random statistics")
  }
),
  private = list(
    set_kmd = function() {
    self$peak_table[, kmd := mz - floor(mz)]
  },
  initialize_data = function() {
    setnames(self$peak_table,
             c("m/z", "Retention time (min)"), c("mz", "rt"))
    self$peak_table <- self$peak_table[which(rowSums(
                                                             self$peak_table[, .SD, .SDcols = self$meta_data$Injection]) > 0), ]
    private$set_kmd()
  }
  ))
