####  filter 1: mismatched peaks    ###
test_that("test that check_mistmatched_peaks works properly with filter_pactr-class data", {
  meta <- data.table(read_csv(test_path("exttestdata", "102623_metadata_correct.csv"), show_col_types = FALSE))
  pt_list <- progenesis_formatter(test_path("exttestdata", "102623_peaktable_coculture_simple.csv"))

  mpactr_class <- mpactr$new(
    pt_list,
    meta
  )
  mpactr_class$setup()
  filter_class <- filter_pactr$new(mpactr_class)
  filter_class$check_mismatched_peaks(
    ringwin = 0.5, isowin = 0.01, trwin = 0.005, max_iso_shift = 3, merge_peaks =
      TRUE, merge_method = "sum"
  )

  expected_cut_ions <- read_csv(test_path("exttestdata", "cut_ions.csv"), col_names = c("V1"), show_col_types = FALSE)
  expected_cut_ions <- as.integer(expected_cut_ions$V1)

  expect_equal(filter_class$logger[["check_mismatched_peaks"]][["cut_ions"]], expected_cut_ions)
  expect_equal(filter_class$mpactr_data$get_peak_table()[
    Compound == "153",
    "102623_UM1850B_ANGDT_71_1_5007"
  ][[1]], 2158.4)
  expect_equal(nrow(filter_class$mpactr_data$get_peak_table()), 1233)
  expect_equal(address(mpactr_class), address(filter_class$mpactr_data))
  expect_equal(mpactr_class$get_peak_table(), filter_class$mpactr_data$get_peak_table())
  expect_false(is.null(filter_class$logger$list_of_summaries$mispicked))
  expect_equal(class(filter_class$logger$list_of_summaries$mispicked), c("summary", "R6"))
})

test_that("test that check_mistmatched_peaks returns an error when no merge method is supplied", {
  meta <- data.table(read_csv(test_path("exttestdata", "102623_metadata_correct.csv"), show_col_types = FALSE))
  pt_list <- progenesis_formatter(test_path("exttestdata", "102623_peaktable_coculture_simple.csv"))

  mpactr_class <- mpactr$new(
    pt_list,
    meta
  )
  mpactr_class$setup()
  filter_class <- filter_pactr$new(mpactr_class)
  expect_error(filter_class$check_mismatched_peaks(ringwin = 0.5, isowin = 0.01, trwin = 0.005, max_iso_shift = 3, merge_peaks = TRUE))
})

####  filter 2: group filter    ###
test_that("blank filter works correctly", {
  meta <- data.table(read_csv(test_path("exttestdata", "102623_metadata_correct.csv"), show_col_types = FALSE))
  pt_list <- progenesis_formatter(test_path("exttestdata", "102623_peaktable_coculture_simple.csv"))

  mpactr_class <- mpactr$new(
    pt_list,
    meta
  )
  mpactr_class$setup()
  filter_class <- filter_pactr$new(mpactr_class)
  filter_class$check_mismatched_peaks(
    ringwin = 0.5, isowin = 0.01, trwin = 0.005, max_iso_shift = 3, merge_peaks =
      TRUE, merge_method = "sum"
  )
  filter_class$filter_blank()

  test_path("exttestdata", "102623_peaktable_coculture_simple_groupaverages.csv")
  error_prop <- read_csv(test_path("exttestdata", "102623_peaktable_coculture_simple_groupaverages.csv"),
    show_col_types = FALSE, skip = 1,
    col_names = c("Compound", "mz", "rt", "biologicalGroup", "average")
  )

  expect_true(all(filter_class$logger[["group_filter-group_stats"]]$Biological_Group %in%
    error_prop$biologicalGroup))

  # group_avgs <- group_avgs[order(group_avgs$Compound, group_avgs$Biological_Group), ]
  expect_true(all(round(filter_class$logger[["group_filter-group_stats"]]$average,
    digits = 5
  ) == round(error_prop$average, digits = 5)))
})

test_that("parse_ions_by_group flags the correct ions", {
  meta <- data.table(read_csv(test_path("exttestdata", "102623_metadata_correct.csv"), show_col_types = FALSE))
  pt_list <- progenesis_formatter(test_path("exttestdata", "102623_peaktable_coculture_simple.csv"))

  mpactr_class <- mpactr$new(
    pt_list,
    meta
  )
  mpactr_class$setup()
  filter_class <- filter_pactr$new(mpactr_class)
  filter_class$check_mismatched_peaks(
    ringwin = 0.5, isowin = 0.01, trwin = 0.005, max_iso_shift = 3, merge_peaks =
      TRUE, merge_method = "sum"
  )
  filter_class$filter_blank()
  filter_class$parse_ions_by_group(group_threshold = 0.01)

  ang_18 <- read_csv(test_path("exttestdata", "output_ANG18_monoculture.csv"),
    col_names = c("V1"),
    show_col_types = FALSE
  )
  angdt <- read_csv(test_path("exttestdata", "output_ANGDT_monoculture"),
    col_names = c("V1"),
    show_col_types = FALSE
  )
  blanks <- read_csv(test_path("exttestdata", "output_Blanks"), col_names = c("V1"), show_col_types = FALSE)
  coculture <- read_csv(test_path("exttestdata", "output_Coculture"), col_names = c("V1"), show_col_types = FALSE)
  jc1 <- read_csv(test_path("exttestdata", "output_JC1_monoculture"), col_names = c("V1"), show_col_types = FALSE)
  jc28 <- read_csv(test_path("exttestdata", "output_JC28_monoculture"), col_names = c("V1"), show_col_types = FALSE)
  group_filter_list <- filter_class$logger[["group_filter-failing_list"]]

  expect_false(all(sapply(group_filter_list, is.null)))
  expect_true(all(group_filter_list$`ANG18 monoculture` == as.character(ang_18$V1)))
  expect_true(all(group_filter_list$`ANGDT monoculture` == as.character(angdt$V1)))
  expect_true(all(group_filter_list$`Blanks` == as.character(blanks$V1)))
  expect_true(all(group_filter_list$`Coculture` == as.character(coculture$V1)))
  expect_true(all(group_filter_list$`JC28 monoculture` == as.character(jc28$V1)))
  expect_true(all(group_filter_list$`JC28 monoculture` == as.character(jc28$V1)))
})

test_that("apply_group_filter removes the correct ions", {
  meta <- data.table(read_csv(test_path("exttestdata", "102623_metadata_correct.csv"), show_col_types = FALSE))
  pt_list <- progenesis_formatter(test_path("exttestdata", "102623_peaktable_coculture_simple.csv"))

  mpactr_class <- mpactr$new(
    pt_list,
    meta
  )
  mpactr_class$setup()
  filter_class <- filter_pactr$new(mpactr_class)
  filter_class$check_mismatched_peaks(
    ringwin = 0.5, isowin = 0.01, trwin = 0.005, max_iso_shift = 3, merge_peaks =
      TRUE, merge_method = "sum"
  )
  filter_class$filter_blank()
  filter_class$parse_ions_by_group(group_threshold = 0.01)

  filter_class$apply_group_filter("Blanks", remove_ions = FALSE)
  expect_equal(nrow(filter_class$mpactr_data$get_peak_table()), 1233)

  filter_class$apply_group_filter("Blanks", remove_ions = TRUE)
  expect_true(all(!(filter_class$logger[["group_filter-failing_list"]]$Blanks %in%
    filter_class$mpactr_data$get_peak_table()$Compound)))

  expect_false(is.null(filter_class$logger$list_of_summaries[["group-Blanks"]]))
  expect_equal(class(filter_class$logger$list_of_summaries[["group-Blanks"]]), c("summary", "R6"))
})

####  filter 3: cv filter    ###
test_that("cv_filter filters out data properly", {
  meta <- data.table(read_csv(test_path("exttestdata", "102623_metadata_correct.csv"), show_col_types = FALSE))
  pt_list <- progenesis_formatter(test_path("exttestdata", "102623_peaktable_coculture_simple.csv"))

  mpactr_class <- mpactr$new(
    pt_list,
    meta
  )
  mpactr_class$setup()
  filter_class <- filter_pactr$new(mpactr_class)
  filter_class$check_mismatched_peaks(
    ringwin = 0.5, isowin = 0.01, trwin = 0.005, max_iso_shift = 3, merge_peaks =
      TRUE, merge_method = "sum"
  )
  filter_class$filter_blank()
  filter_class$parse_ions_by_group(group_threshold = 0.01)
  filter_class$apply_group_filter("Blanks", remove_ions = TRUE)
  filter_class_median <- filter_class$clone(deep = TRUE)
  filter_class$cv_filter(cv_threshold = 0.2, cv_params = c("mean"))
  cv_filter_passed_ions <- filter_class$logger[["list_of_summaries"]]$replicability$get_passed_ions()
  expect_equal(length(filter_class$logger[["list_of_summaries"]]$replicability$get_failed_ions()), 86)
  filter_class_median$cv_filter(cv_threshold = 0.2, cv_params = c("median"))
  cv_filter_passed_ions_median <- filter_class_median$logger[["list_of_summaries"]]$replicability$get_passed_ions()
  expect_equal(length(filter_class_median$logger[["list_of_summaries"]]$replicability$get_failed_ions()), 61)
  expect_false(length(cv_filter_passed_ions) == length(cv_filter_passed_ions_median))

  expect_false(is.null(filter_class$logger$list_of_summaries$replicability))
  expect_equal(class(filter_class$logger$list_of_summaries$replicability), c("summary", "R6"))
})
test_that("cv_filter errors without threshold", {
  meta <- data.table(read_csv(test_path("exttestdata", "102623_metadata_correct.csv"), show_col_types = FALSE))
  pt_list <- progenesis_formatter(test_path("exttestdata", "102623_peaktable_coculture_simple.csv"))

  mpactr_class <- mpactr$new(
    pt_list,
    meta
  )
  mpactr_class$setup()
  filter_class <- filter_pactr$new(mpactr_class)
  filter_class$check_mismatched_peaks(
    ringwin = 0.5, isowin = 0.01, trwin = 0.005, max_iso_shift = 3, merge_peaks =
      TRUE, merge_method = "sum"
  )
  filter_class$filter_blank()
  filter_class$parse_ions_by_group(group_threshold = 0.01)
  filter_class$apply_group_filter("Blanks", remove_ions = TRUE)

  expect_error(filter_class$cv_filter(cv_params = c("mean")))
})

test_that("cv_filter errors with incorrect paramter", {
  meta <- data.table(read_csv(test_path("exttestdata", "102623_metadata_correct.csv"), show_col_types = FALSE))
  pt_list <- progenesis_formatter(test_path("exttestdata", "102623_peaktable_coculture_simple.csv"))

  mpactr_class <- mpactr$new(
    pt_list,
    meta
  )
  mpactr_class$setup()
  filter_class <- filter_pactr$new(mpactr_class)
  filter_class$check_mismatched_peaks(
    ringwin = 0.5, isowin = 0.01, trwin = 0.005, max_iso_shift = 3, merge_peaks =
      TRUE, merge_method = "sum"
  )
  filter_class$filter_blank()
  filter_class$parse_ions_by_group(group_threshold = 0.01)
  filter_class$apply_group_filter("Blanks", remove_ions = TRUE)

  expect_error(filter_class$cv_filter(cv_threshold = 0.2, cv_params = ""))
})

test_that("cv_filter errors when there are no technical replicates", { # hmm
  meta <- data.table(read_csv(test_path("exttestdata", "102623_metadata_correct.csv"), show_col_types = FALSE))
  meta_sub <- meta[, head(.SD, 1), by = Sample_Code]

  pt_list <- progenesis_formatter(test_path("exttestdata", "102623_peaktable_coculture_simple.csv"))
  sub <- c("Compound", "mz", "rt", meta_sub$Injection)
  pt_list$peak_table <- pt_list$peak_table[, .SD, .SDcols = sub]
  pt_list$raw_table <- pt_list$raw_table[, .SD, .SDcols = sub]


  mpactr_class <- mpactr$new(
    pt_list,
    meta_sub
  )
  mpactr_class$setup()

  filter_class <- filter_pactr$new(mpactr_class)
  filter_class$check_mismatched_peaks(
    ringwin = 0.5, isowin = 0.01, trwin = 0.005, max_iso_shift = 3, merge_peaks =
      TRUE, merge_method = "sum"
  )
  filter_class$filter_blank()
  filter_class$parse_ions_by_group(group_threshold = 0.01)
  filter_class$apply_group_filter("Blanks", remove_ions = TRUE)


  expect_error(filter_class$cv_filter(cv_threshold = 0.2, cv_params = "mean"))
})




####  filter 4: insource ions    ###
test_that("filter_inscource_ions filters out data properly", {
  meta <- data.table(read_csv(test_path("exttestdata", "102623_metadata_correct.csv"), show_col_types = FALSE))
  pt_list <- progenesis_formatter(test_path("exttestdata", "102623_peaktable_coculture_simple.csv"))

  mpactr_class <- mpactr$new(
    pt_list,
    meta
  )
  mpactr_class$setup()
  filter_class <- filter_pactr$new(mpactr_class)
  filter_class$check_mismatched_peaks(
    ringwin = 0.5, isowin = 0.01, trwin = 0.005, max_iso_shift = 3, merge_peaks =
      TRUE, merge_method = "sum"
  )
  filter_class$filter_blank()
  filter_class$parse_ions_by_group(group_threshold = 0.01)
  filter_class$apply_group_filter("Blanks", remove_ions = TRUE)
  filter_class$filter_insource_ions(cluster_threshold = 0.95)

  insource_ion_expected_list <- c(
    38, 204, 214, 993, 270, 1003, 271, 294, 331, 349, 382,
    447, 498, 1233, 644, 1307, 677, 675, 689,
    690, 688, 758, 985, 982, 981, 1297, 1311
  )
  expect_true(length(filter_class$logger[["list_of_summaries"]]$insource$get_failed_ions()) == 27)
  expect_true(all(!(insource_ion_expected_list %in% filter_class$mpactr_data$get_peak_table()$Compound)))

  expect_false(is.null(filter_class$logger$list_of_summaries$insource))
  expect_true(is.null(filter_class$logger$list_of_summaries$replicability))
  expect_equal(class(filter_class$logger$list_of_summaries$insource), c("summary", "R6"))
})
