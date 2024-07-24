test_that("mpactr class setup works properly", {
  meta <- data.table(read_csv(test_path("exttestdata", "102623_metadata_correct.csv"), show_col_types = FALSE))
  pt_list <- progenesis_formatter(test_path("exttestdata", "102623_peaktable_coculture_simple.csv"))

  mpactr_class <- mpactr$new(
    pt_list,
    meta
  )
  mpactr_class$setup()

  expect_true(all(c("kmd", "mz", "rt") %in% colnames(mpactr_class$get_peak_table())))
  expect_equal(nrow(mpactr_class$get_peak_table()), 1303)
})
