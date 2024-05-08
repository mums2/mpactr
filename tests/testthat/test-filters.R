test_that("Test that solvent blank filter removes correct columns", {
  test_ex_df <- data.frame("Compound" = c("1", "2", "3", "4", "5"),
                           "mz" = c(43.8, 46.2, 89.4, 101.8, 169.0),
                           "rt" = c(0.7748, 0.7847, 0.8, 0.8132, 0.8422),
                           "sample1_1" = c(0, 344, 592, 4, 1),
                           "sample1_2" = c(1, 350, 600, 8, 3),
                           "blank_1" = c(0, 0, 0, 3, 0),
                           "blank_2" = c(0, 0, 0, 0, 1))
  test_ex_df <- solvent_blank_filter(test_ex_df, c("blank_1", "blank_2"))
  expected_result <- data.frame("Compound" = c("1", "2", "3"),
                                "mz" = c(43.8, 46.2, 89.4),
                                "rt" = c(0.7748, 0.7847, 0.8),
                                "sample1_1" = c(0, 344, 592),
                                "sample1_2" = c(1, 350, 600))
  expect_equal(test_ex_df, expected_result)
})


test_that("mismatched peaks filter works", {
  peak_df <- readr::read_csv(here::here("tests/exttestdata/102623 peaktable coculture simple.csv"), skip = 2)
  colnames(peak_df)[which(colnames(peak_df) %in% c("m/z"))] <- "mz"
  colnames(peak_df)[which(colnames(peak_df) %in% c("Retention time (min)"))] <- "rt"

  peak_df <- initialize_data(peak_df)
  
  sample_df <- readr::read_csv(here::here("tests/exttestdata/102623 samplelist.csv"))
  meta <- readr::read_csv(here::here("tests/exttestdata/102623 metadata simple.csv"))
  
  relfil_ion_list <- check_mismatched_peaks(peak_df, ringwin = 0.5, isowin = 0.01, trwin = 0.005, max_iso_shift = 3, merge_peaks = FALSE)

  expected_cut_ions <- read.csv(here::here("tests/exttestdata/cut_ions.csv"), header = FALSE)
  expected_cut_ions <- as.integer(expected_cut_ions$V1)

  expect_equal(relfil_ion_list$cut_ions, expected_cut_ions)
})
