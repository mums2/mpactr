test_that("import_data creates a proper mpactr and filter-pactr object", {

  data <- import_data(test_path("exttestdata","102623_peaktable_coculture_simple.csv"),
                      test_path("exttestdata", "102623_metadata_correct.csv"))

  expect_true(all(class(data) == c("filter_pactr", "R6")))
  expect_true(exists("list_of_summaries", data$logger))
  expect_true(all(sapply(data$logger[["list_of_summaries"]], is.null) == TRUE))
  expect_true(nrow(data$mpactr_data$get_peak_table()) > 1)
  expect_true(nrow(data$mpactr_data$get_meta_data()) > 1)
})
