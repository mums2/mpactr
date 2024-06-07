graph_qc_pactr <- R6Class("graph_qc_pactr", public = list(
  initialize = function(filter_pactr) {
    private$filter_pactr_data = filter_pactr
  },
  generate_QC_Summary = function()
  {
    list_of_failed_ions <- lapply(private$filter_pactr_data$logger$list_of_summaries, function(x) {
      x$get_failed_ions()
    })
    status <- gsub('[[:digit:]]+', '', names(unlist(list_of_failed_ions)))
    compounds <- unlist(unname(list_of_failed_ions))

    failed_ions_dt <- data.table(status, compounds)
    
    if (nrow(failed_ions_dt) != length(unique(failed_ions_dt$compounds))) {
      cli::cli_abort("At least one ion fails multiple filters. Please check.")
    }
    
    passed_ions_dt <- data.table("status" = rep("Passed", length
     (private$filter_pactr_data$mpactr_data$get_peak_table()$Compound)),
                    "compounds" = private$filter_pactr_data$mpactr_data$get_peak_table()$Compound)
                    
    private$filter_summarized <- rbind(passed_ions_dt, failed_ions_dt)
  },
  plot_QC_Tree = function()
  {
    ion_counts <- private$filter_summarized[ , .(count = .N), by = status][
       , percent := (count / sum(count) * 100)]

    return(ggplot(ion_counts) +
      aes(area = percent, fill = status) +
      geom_treemap() +
      geom_treemap_text(aes(label = paste(status, paste0(round(percent, 2), "%"), sep = "\n")), colour = "darkorchid1",
      fontface = c("bold")) +
      theme(legend.position = "none") +
      scale_fill_viridis(option = "G", discrete = TRUE))
  },
  get_summarized_dt = function()
  {
    return(private$filter_summarized)
  }
), private = list(
    filter_summarized = NA,
    filter_pactr_data = NA
))
