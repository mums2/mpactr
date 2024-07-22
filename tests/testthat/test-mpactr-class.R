test_that("mpactr class initialize works correctly", {
  meta <- data.table(read_csv(test_path("exttestdata", "102623_metadata_correct.csv"), show_col_types = FALSE))
  pt_list <- progenesis_formatter(test_path("exttestdata","102623_peaktable_coculture_simple.csv"))
  
  mpactr_class <- mpactr$new(pt_list,
                             meta)
  
  expect_true(all(class(mpactr_class) == c("mpactr", "R6")))
  expect_equal(length(mpactr_class$get_meta_data()), 3)
  expect_equal(length(mpactr_class$get_peak_table()), 21)
  expect_error(mpactr$new(peak_table_path = 2, meta_data_path = 5), NULL)
})

test_that("mpactr isMultipleTechReps correctly dtermines if there are technical replicates", {
    meta <- data.table(read_csv(test_path("exttestdata", "102623_metadata_correct.csv"), show_col_types = FALSE))
  pt_list <- progenesis_formatter(test_path("exttestdata","102623_peaktable_coculture_simple.csv"))
  
  mpactr_class <- mpactr$new(pt_list,
                             meta)
  expect_true(mpactr_class$isMultipleTechReps())
  # t <- c(1, 3, 3, 3, 3)
})


