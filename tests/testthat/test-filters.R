test_that("Test that solvent blank filter removes correct columns", {
  test_ex_df <- data.frame("sample1_1" = c(0, 344, 592, 4, 1),
                           "sample1_2" = c(1, 350, 600, 8, 3),
                           "blank_1" = c(0, 0, 0, 3, 0),
                           "blank_2" = c(0, 0, 0, 0, 1))
  test_ex_df <- solvent_blank_filter(test_ex_df, c("blank_1", "blank_2"))
  expected_result <- data.frame("sample1_1" = c(0, 344, 592),
                                "sample1_2" = c(1, 350, 600))
  expect_equal(test_ex_df, expected_result)
})
