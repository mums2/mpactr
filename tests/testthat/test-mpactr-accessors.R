test_that("get_peak_table correctly returns the table", {
  data <- import_data(test_path("exttestdata","102623_peaktable_coculture_simple.csv"),
                      test_path("exttestdata", "102623_metadata_correct.csv"))
  
  peak_table <- get_peak_table(data)
  expect_equal(class(peak_table), c("data.table", "data.frame"))
  expect_equal(nrow(peak_table), 1303)
})

test_that("get_meta_data correctly returns the table", {
  data <- import_data(test_path("exttestdata","102623_peaktable_coculture_simple.csv"),
                      test_path("exttestdata", "102623_metadata_correct.csv"))
  
  metadata <- get_meta_data(data)
  expect_equal(class(metadata), c("data.table", "data.frame"))
  expect_equal(nrow(metadata), 18)
})