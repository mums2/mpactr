test_that("Set kmd properly creates a new kmd column", {
  test_ex_df <- data.frame("Compund" = c("1", "2", "3", "4", "5"),
                           "mz" = c(43.837890, 46.23978, 89.40925,
                                    101.87653, 169.02986),
                           "rt" = c(0.7748, 0.7847, 0.8, 0.8132, 0.8422),
                           "sample1_1" = c(0, 344, 592, 4, 1),
                           "sample1_2" = c(1, 350, 600, 8, 3),
                           "blank_1" = c(0, 0, 0, 3, 0),
                           "blank_2" = c(0, 0, 0, 0, 1))

  test_ex_df <- set_kmd(test_ex_df)
  expected_df <- data.frame("Compund" = c("1", "2", "3", "4", "5"),
                            "mz" = c(43.837890, 46.23978, 89.40925,
                                     101.87653, 169.02986),
                            "rt" = c(0.7748, 0.7847, 0.8, 0.8132, 0.8422),
                            "sample1_1" = c(0, 344, 592, 4, 1),
                            "sample1_2" = c(1, 350, 600, 8, 3),
                            "blank_1" = c(0, 0, 0, 3, 0),
                            "blank_2" = c(0, 0, 0, 0, 1),
                            "kmd" = c(0.83789, 0.23978, 0.40925, 0.87653,
                                      0.02986))

  expect_equal(test_ex_df, expected_df)
})
