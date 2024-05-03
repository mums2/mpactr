test_that("Test that solvent blank filter removes correct columns", {
  test_ex_df <- data.frame("Compund" = c("1", "2", "3", "4", "5"),
                           "mz" = c(43.8, 46.2, 89.4, 101.8, 169.0),
                           "rt" = c(0.7748, 0.7847, 0.8, 0.8132, 0.8422),
                           "sample1_1" = c(0, 344, 592, 4, 1),
                           "sample1_2" = c(1, 350, 600, 8, 3),
                           "blank_1" = c(0, 0, 0, 3, 0),
                           "blank_2" = c(0, 0, 0, 0, 1))
  test_ex_df <- solvent_blank_filter(test_ex_df, c("blank_1", "blank_2"))
  expected_result <- data.frame("Compund" = c("1", "2", "3"),
                                "mz" = c(43.8, 46.2, 89.4),
                                "rt" = c(0.7748, 0.7847, 0.8),
                                "sample1_1" = c(0, 344, 592),
                                "sample1_2" = c(1, 350, 600))
  expect_equal(test_ex_df, expected_result)
})
