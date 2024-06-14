
test_that("filter_summary returns an error when an incorrect fitler argument is provided", {
  data <- import_data(test_path("exttestdata","102623_peaktable_coculture_simple.csv"),
                      test_path("exttestdata", "102623_metadata_correct.csv"))

  expect_error(filter_summary(data, filter = "cv"), "`filter` must be one of mpactR's")
  
})

test_that("filter_summary returns an error when the fitler argument provided has not yet been run (e.g., not in the log)", {
  data <- import_data(test_path("exttestdata","102623_peaktable_coculture_simple.csv"),
                      test_path("exttestdata", "102623_metadata_correct.csv"))
  
  data_mpactr <- filter_mispicked_ions(data, ringwin = 0.5, isowin = 0.01, trwin = 0.005, max_iso_shift = 3, merge_peaks = TRUE)
  
  expect_error(filter_summary(data_mpactr, filter = "replicability"), "`filter` replicability has not yet been applied to the data")
})

test_that(" returns the correct fitler summary list", {
  data <- import_data(test_path("exttestdata","102623_peaktable_coculture_simple.csv"),
                      test_path("exttestdata", "102623_metadata_correct.csv"))
  
  data_mpactr <- filter_mispicked_ions(data, ringwin = 0.5, isowin = 0.01, trwin = 0.005, max_iso_shift = 3, merge_peaks = TRUE)
  
  mispicked_summary <- filter_summary(data_mpactr, filter = "mispicked")
  
  expect_type(mispicked_summary, "list")
  expect_equal(length(mispicked_summary), 2)
  expect_equal(length(mispicked_summary$passed_ions), 1233)
})
