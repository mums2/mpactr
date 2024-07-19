test_that("progensis_formatter reads data correctly", {
  pt <- progensis_formatter(test_path("exttestdata","102623_peaktable_coculture_simple.csv"))
  
  expect_equal(class(pt$peak_table), c("data.table", "data.frame"))
  expect_true(all(c("data.table", "data.frame") %in% class(pt$peak_table)))
  expect_true(nrow(pt$peak_table) > 0)
  
})
test_that("metaboscape reads data correctly", {

  samples <- c("UM1850B_ANGDT_0.25_mL_36_1_4792" , "UM1850B_ANGDT_0.25_mL_36_1_4802"  ,"UM1850B_ANGDT_0.25_mL_36_1_4815"  ,"MixedMonoculture_1mg_mL total_31_1_4799",
  "MixedMonoculture_1mg_mL total_31_1_4822","MixedMonoculture_31_1_4787","UM1852B_Coculture_42_1_4685", "UM1852B_Coculture_42_1_4695", "UM1852B_Coculture_42_1_4709", "UM1852B_Coculture_42_1_4785")        

  pt <- metaboscape_formatter(test_path("../../inst/extdata", "MJB_MonoVsCoculture_metaboscape_ft.csv"), sample_names = samples)
  
  expect_equal(class(pt$peak_table), c("data.table", "data.frame"))
  expect_true(all(c("data.table", "data.frame") %in% class(pt$peak_table)))
  expect_true(nrow(pt$peak_table) > 0)
  expect_equal(colnames(pt$peak_table), c("Compound", "mz", "rt", samples))
  
})
