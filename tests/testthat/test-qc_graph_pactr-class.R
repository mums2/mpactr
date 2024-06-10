test_that("plot_QC_Tree generates the correct plot", {
    
  mpactr_class <- mpactr$new(here::here("tests/exttestdata/102623 peaktable coculture simple.csv"),
                            here::here("tests/exttestdata/102623_metadata_correct.csv"))
  mpactr_class$setup()
  filter_class <- filter_pactr$new(mpactr_class)
  filter_class$check_mismatched_peaks(ringwin = 0.5, isowin = 0.01, trwin = 0.005, max_iso_shift = 3, merge_peaks =
    TRUE)
  filter_class$filter_blank()
  filter_class$parse_ions_by_group(group_threshold = 0.01)
  filter_class$apply_group_filter("Blanks", remove_ions = TRUE)
  filter_class$filter_insource_ions(cluster_threshold = 0.95)

  graph_qc_pactr_class <- graph_qc_pactr$new(filter_class)
  graph_qc_pactr_class$generate_QC_Summary()
  plot <- graph_qc_pactr_class$plot_QC_Tree()
  
  expect_equal(class(plot), c("gg", "ggplot"))

  # ggsave(here::here("tests/exttestdata/actual_treemap.png"), plot)
  # expect_snapshot_file(ggsave(here::here("tests/exttestdata/actual_treemap.png"), plot), 
  #                      here::here("tests/exttestdata/expected_treemap.png"))
})

# ion_status <- graph_qc_pactr_class$get_summarized_dt()
# i <- ion_status[ , .(count = .N), by = status][ , percent := (count / sum(count) * 100)]

# library(treemapify)
# library(ggplot2)
# plot <- ggplot(i) +
#   aes(area = percent, fill = status) +
#   geom_treemap() +
#   geom_treemap_text(aes(label = paste(status, paste0(round(percent, 2), "%"), sep = "\n")), colour = "darkorchid1",
#   fontface = c("bold")) +
#   theme(legend.position = "none") +
#   scale_fill_viridis(option = "G", discrete = TRUE) 
  
# ggsave(here::here("tests/exttestdata/expected_treemap.png"), plot)
# saveRDS(plot, here::here("tests/exttestdata/expected_treemap.RDS"))

# expected_plot <- readRDS(here::here("tests/exttestdata/expected_treemap.RDS"))
# identical(plot, expected_plot)

# all(plot == expected_plot)
# plot$