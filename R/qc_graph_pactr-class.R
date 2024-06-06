graph_qc_pactr <- R6Class("graph_qc_pactr", public = list(
  filter_pactr_data = NA,
  initialize = function(filter_pactr){
    self$filter_pactr_data = filter_pactr
  },
  generate_QC_Summary = function()
  {
    # logger[["list_of_summaries"]] < c("cv_filter_summary", "in_source_ion_summary")
    # cv_filtered failed x ions
    # blank filter failed x ions?
  },
  create_QC_Tree_Plot = function()
  {

  }
))


 # be more efficent
# overall QC summary
  # ions failing X filter
# graph tree plot
  # compound, status: mismatch, group, cv, decon, pass
# 1 mismatch
# 2 pass
# 3 pass


