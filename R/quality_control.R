qc_summary <- function(mpactr_object)
{
  graph_pactr_object <- graph_qc_pactr$new(mpactr_object)
  graph_pactr_object$generate_QC_Summary()
  return(graph_pactr_object$get_summarized_dt())
}

plot_qc_tree <- function(mpactr_object)
{
  graph_pactr_object <- graph_qc_pactr$new(mpactr_object)
  graph_pactr_object$generate_QC_Summary()
  return(graph_pactr_object$plot_QC_Tree())
}