graph_qc_pactr <- R6Class("graph_qc_pactr", public = list(
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

# env <- new.env(HASH = TRUE)
# env["statistic_one"] <- 
# env["statistic_two"] <- 
# env["statistic_three"] <-
# env["statistics_four"] <- list_of_compounds_effected

# x <- function(x)
# {
#   data1 <- env["statistic_one"] 
# }
# I could also give this class the mpact object and use this to call all of the filters on the mpact data.
# That might be more efficent