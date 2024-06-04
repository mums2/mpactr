mpactr <- R6Class("mpactr", public = list(
  # Properties
  peak_table_csv = NA,
  meta_data_csv = NA,
  sample_list_csv = NA,
  # Constructor
  initialize = function(peak_table_path, meta_data_path, sample_list_path){
    self$peak_table_csv <- readr::read_csv(here::here(peak_table_path), skip =2)
    self$meta_data_csv <- readr::read_csv(here::here(meta_data_path))
    self$sample_list_csv <- readr::read_csv(here::here(sample_list_path))
  },
  setup = function() {
    initialize_data(self$peak_table_csv) # TODO Redo functions to deal with references
  },
  filter_data = function(filter_function, params = ...) {
    summary <- filter_function(self$peak_table_csv, params)
    self$logger(summary)
  },
  logger = function(summary){
    print(summary)
    print("Some other random statistics")
  },
  transform_data = function()
  {
    for(i in 1:nrow(self$peak_table_csv))
    {
      self$peak_table_csv$`m/z`[i] <- self$peak_table_csv$`m/z`[i] * 2
    }
  }
))

# For this to work, we need to ensure each function modifies the references of data within this class.
# This will allow us to make the code even faster and manipulate the data the user will use all inside of a
# Simple class function
# From here there are many other things we can do


graph_pactr <- R6Class("graph_pactr", public = list(
  mpactr_data = NA,
  initialize = function(mpactr){
    self$mpactr_data <- mpactr
  },
  bar_graph = function()
  {
    # Graph a bar graph
    # self$mpactr_data graph
  },
  volcano_plot = function()
  {
    # Graph a volcano_plot
    # self$mpactr_data graph
  }
))

# I could also give this class the mpact object and use this to call all of the filters on the mpact data.
# That might be more efficent

filter_pactr <- R6Class("filter_pactr", public = list(
  initialize = function(val){
    print(val)
  },
  cv_filter = function(x)
  {
    # etc
  },
  merge_mismtached_peaks = function(mpactr_object, ringwin, isowin, trwin, max_iso_shift, merge_peaks)
  {
    check_mismatched_peaks(mpactr_object, ringwin, isowin, trwin, max_iso_shift, merge_peaks)
    # We would pass the mpact object over to the class (as a reference) and work with the reference of it. Therefore,
    # there will be no copy and it will be passing by a reference.
  }
))



mpact_class <- mpactr$new("tests/exttestdata/102623 peaktable coculture simple.csv",
                          "tests/exttestdata/102623 metadata simple.csv",
                          "tests/exttestdata/102623 samplelist.csv")

by_val(mpact_class)
peak_table <- readr::read_csv(here::here("tests/exttestdata/102623 peaktable coculture simple.csv"), skip = 2)
by_val(peak_table)
mpact_class$transform_data()

by_val <- function(df)
{
  # rows <- nrow(df)
    for(i in 1:nrow(df$peak_table_csv))
    {
      df$peak_table_csv$`m/z`[i] <- df$peak_table_csv$`m/z`[i] * 2
    }
}
library(microbenchmark)
microbenchmark::microbenchmark(by_val(mpact_class$clone()),mpact_class$transform_data(), times = 10)

