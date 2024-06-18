test_that("filter mismatch ions wrapper works as expected when merge_peaks is TRUE", {
  data <- import_data(test_path("exttestdata","102623_peaktable_coculture_simple.csv"),
                      test_path("exttestdata", "102623_metadata_correct.csv"))

  data_mpactr_copy <- filter_mispicked_ions(data, ringwin = 0.5, isowin = 0.01, trwin = 0.005, max_iso_shift = 3,
                                      merge_peaks = TRUE, copy_object = TRUE)

  data_mpactr <- filter_mispicked_ions(data, ringwin = 0.5, isowin = 0.01, trwin = 0.005, max_iso_shift = 3,
                                      merge_peaks = TRUE, copy_object = FALSE)

  expected_cut_ions <- read_csv(test_path("exttestdata", "cut_ions.csv"), col_names = c("V1"), show_col_types = FALSE)
  expected_cut_ions <- as.integer(expected_cut_ions$V1)

  expect_equal(nrow(data_mpactr$mpactr_data$get_peak_table()), nrow((data_mpactr$mpactr_data$get_peak_table())))

  expect_equal(data_mpactr$logger[["check_mismatched_peaks"]][["cut_ions"]], expected_cut_ions)
  expect_equal(data_mpactr_copy$logger[["check_mismatched_peaks"]][["cut_ions"]], expected_cut_ions)

  expect_equal(nrow(data_mpactr$mpactr_data$get_peak_table()), 1233)
})

test_that("filter mismatch ions wrapper works as expected when merge_peaks is FALSE", {
  data <- import_data(test_path("exttestdata","102623_peaktable_coculture_simple.csv"),
                      test_path("exttestdata", "102623_metadata_correct.csv"))

   data_mpactr_copy <- filter_mispicked_ions(data, ringwin = 0.5, isowin = 0.01, trwin = 0.005, max_iso_shift = 3,
                                        merge_peaks = FALSE, copy_object = TRUE)

   data_mpactr <- filter_mispicked_ions(data, ringwin = 0.5, isowin = 0.01, trwin = 0.005, max_iso_shift = 3,
                                        merge_peaks = FALSE)

  expect_equal(nrow(data_mpactr$mpactr_data$get_peak_table()), 1303)
  expect_equal(nrow(data_mpactr_copy$mpactr_data$get_peak_table()), 1303)
})

test_that("group filter wrapper works as expected", {
  data <- import_data(test_path("exttestdata","102623_peaktable_coculture_simple.csv"),
                      test_path("exttestdata", "102623_metadata_correct.csv"))
  data_mpactr <- filter_mispicked_ions(data, ringwin = 0.5, isowin = 0.01, trwin = 0.005, max_iso_shift = 3, merge_peaks =
    TRUE)
  data_mpactr <- filter_group(data_mpactr, 0.01, "Blanks", FALSE)
  expect_equal(nrow(data_mpactr$mpactr_data$get_peak_table()), 1233)

  data_mpactr_copy <- filter_group(data_mpactr, 0.01, "Blanks", TRUE,
                                   copy_object = TRUE)

  data_mpactr <- filter_group(data_mpactr, 0.01, "Blanks", TRUE)
  expect_equal(nrow(data_mpactr$mpactr_data$get_peak_table()), 484)

  expect_equal(nrow(data_mpactr$mpactr_data$get_peak_table()),
               nrow(data_mpactr_copy$mpactr_data$get_peak_table()))

  expect_true(all(!(data_mpactr$logger[["group_filter-failing_list"]]$Blanks %in%
    data_mpactr$mpactr_data$get_peak_table()$Compound)))

   expect_true(all(!(data_mpactr_copy$logger[["group_filter-failing_list"]]$Blanks %in%
    data_mpactr_copy$mpactr_data$get_peak_table()$Compound)))
})

test_that("filter cv filter wrapper works as expected with cv_params mean", {
  data <- import_data(test_path("exttestdata","102623_peaktable_coculture_simple.csv"),
                      test_path("exttestdata", "102623_metadata_correct.csv"))

  data_mpactr <- filter_mispicked_ions(data, ringwin = 0.5, isowin = 0.01, trwin = 0.005, max_iso_shift = 3, merge_peaks =
    TRUE)
  data_mpactr <- filter_group(data_mpactr, 0.01, "Blanks", TRUE)
  data_mpactr_copy <- filter_cv(data_mpactr, 0.2 , "mean", copy_object = TRUE)
  data_mpactr <- filter_cv(data_mpactr, 0.2 , "mean")


  expect_equal(length(data_mpactr$logger[["list_of_summaries"]]$replicability$get_failed_ions()), 86)
  expect_equal(length(data_mpactr_copy$logger[["list_of_summaries"]]$replicability$get_failed_ions()), 86)
})

test_that("filter cv filter wrapper works as expected with cv_params median", {
  data <- import_data(test_path("exttestdata","102623_peaktable_coculture_simple.csv"),
                      test_path("exttestdata", "102623_metadata_correct.csv"))

  data_mpactr <- filter_mispicked_ions(data, ringwin = 0.5, isowin = 0.01, trwin = 0.005, max_iso_shift = 3, merge_peaks =
    TRUE)
  data_mpactr <- filter_group(data_mpactr, 0.01, "Blanks", TRUE)

  data_mpactr_copy <- filter_cv(data_mpactr, 0.2 , "median", copy_object = TRUE)
  data_mpactr <- filter_cv(data_mpactr, 0.2 , "median")

  expect_equal(length(data_mpactr$logger[["list_of_summaries"]]$replicability$get_failed_ions()), 61)

  expect_equal(length(data_mpactr$logger[["list_of_summaries"]]$replicability$get_failed_ions()),
               length(data_mpactr_copy$logger[["list_of_summaries"]]$replicability$get_failed_ions()))
})

test_that("filter insource ions wrapper works as expected", {
  data <- import_data(test_path("exttestdata","102623_peaktable_coculture_simple.csv"),
                      test_path("exttestdata", "102623_metadata_correct.csv"))

  data_mpactr <- filter_mispicked_ions(data, ringwin = 0.5, isowin = 0.01, trwin = 0.005, max_iso_shift = 3, merge_peaks =
    TRUE)
  data_mpactr <- filter_group(data_mpactr, 0.01, "Blanks", TRUE)

  data_mpactr_copy <- filter_insource_ions(data_mpactr, cluster_threshold = 0.95, copy_object = TRUE)

  data_mpactr <- filter_insource_ions(data_mpactr, cluster_threshold = 0.95)

  insource_ion_expected_list <- c(38, 204, 214, 993, 270, 1003, 271, 294, 331, 349, 382,
   447, 498, 1233, 644, 1307, 677, 675, 689,
   690, 688, 758, 985, 982, 981, 1297, 1311)

  expect_true(length(data_mpactr$logger[["list_of_summaries"]]$insource$get_failed_ions()) == 27)
  expect_true(length(data_mpactr_copy$logger[["list_of_summaries"]]$insource$get_failed_ions()) == 27)

  expect_true(all(!(insource_ion_expected_list %in% data_mpactr$mpactr_data$get_peak_table()$Compound)))
  expect_true(all(!(insource_ion_expected_list %in% data_mpactr_copy$mpactr_data$get_peak_table()$Compound)))
})