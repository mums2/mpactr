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

############################################
####     filter 1: mismatched peaks     ####
############################################
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

test_that("merge_ions correctly updates intensity values and columns", {
  peak_df <- readr::read_csv(here::here("tests/exttestdata/102623 peaktable coculture simple.csv"), skip = 2)
  colnames(peak_df)[which(colnames(peak_df) %in% c("m/z"))] <- "mz"
  colnames(peak_df)[which(colnames(peak_df) %in% c("Retention time (min)"))] <- "rt"

  peak_df <- initialize_data(peak_df)
  
  sample_df <- readr::read_csv(here::here("tests/exttestdata/102623 samplelist.csv"))
  meta <- readr::read_csv(here::here("tests/exttestdata/102623 metadata simple.csv"))
  
  relfil_ion_list <- check_mismatched_peaks(peak_df, ringwin = 0.5, isowin = 0.01, trwin = 0.005, max_iso_shift = 3, merge_peaks = FALSE)
  
  peak_df_merged <- merge_ions(peak_df, relfil_ion_list)
  expect_true(all(!(as.character(relfil_ion_list$cut_ions) %in% colnames(peak_df_merged))))

  expect_equal(peak_df_merged[peak_df_merged$Compound == "153", "102623_UM1850B_ANGDT_71_1_5007"], 2158.4)
  
    
  expected_merged <- read_csv(here::here("tests/exttestdata/102623 peaktable coculture simple_merged.csv"), skip = 2) %>%
    as.data.frame() %>%
    select(-`...2`, -`...3`)
  expected_merged <- expected_merged[expected_merged$Compound %in% names(relfil_ion_list$merge_groups), ]
  
  merged_results_rows <- peak_df_merged[peak_df_merged$Compound %in% names(relfil_ion_list$merge_groups), ] %>%
    select(-mz, -rt, -kmd)
    
  # all(merged_results_rows == expected_merged)
  
  # dim <- c()
  # for(i in 1:nrow(merged_results_rows)) {
  #   dim <- c(dim, all(expected_merged[i, ] == merged_results_rows[i, ]))
  #   print(all(expected_merged[i, ] == merged_results_rows[i, ]))
  # }
  # table(dim)
  # merged_results_rows$Compound[which(dim == FALSE)]
  
  # merged_results_rows[2, "102623_UM1849B_ANG18_70_1_5021"] == expected_merged[2, "102623_UM1849B_ANG18_70_1_5021"]
  # print(merged_results_rows[2, "102623_UM1849B_ANG18_70_1_5021"], digits = 10)
  # print(expected_merged[2, "102623_UM1849B_ANG18_70_1_5021"], digits = 10)
  
  # group_avgs %>%
  #   filter(Compound == "1059") %>%
  #   select(Compound, Biological_Group, average)
})

############################################
####       filter 2: group/blank        ####
############################################

test_that("group/blank filtering works correctly", {
  # the blank filter is passed merged peaks from filter 1 (mistmatched peaks)
  peak_df <- readr::read_csv(here::here("tests/exttestdata/102623 peaktable coculture simple.csv"), skip = 2)
  colnames(peak_df)[which(colnames(peak_df) %in% c("m/z"))] <- "mz"
  colnames(peak_df)[which(colnames(peak_df) %in% c("Retention time (min)"))] <- "rt"

  peak_df <- initialize_data(peak_df)
  
  sample_df <- readr::read_csv(here::here("tests/exttestdata/102623 samplelist.csv"))
  meta <- readr::read_csv(here::here("tests/exttestdata/102623 metadata simple.csv"))
  
  full_meta <- sample_df %>% left_join(meta, by = "Sample_Code") %>%
    filter(Biological_Group != "NA") %>%
    select(Injection, Sample_Code, Biological_Group) 
  
  peak_df_relfil <- check_mismatched_peaks(peak_df, ringwin = 0.5, isowin = 0.01, trwin = 0.005, max_iso_shift = 3, merge_peaks = TRUE)
  
  group_avgs <- filter_blank(peak_df_relfil, full_meta)
  error_prop <- read.csv(here::here("tests/exttestdata/102623 peaktable coculture simple_groupaverages.csv"), header = TRUE)
  colnames(error_prop) <- c("Compound", "mz", "rt", "biologicalGroup", "average")
  error_prop$Compound <- as.character(error_prop$Compound)
  error_prop <- error_prop[order(error_prop$Compound), ]
  expect_true(all(group_avgs$Biological_Group == error_prop$biologicalGroup))
  expect_true(all(round(group_avgs$average, digits = 5) == round(error_prop$average, digits = 5)))
})


test_that("group/blank filtering3 works correctly", {
  # the blank filter is passed merged peaks from filter 1 (mistmatched peaks)
  peak_df <- readr::read_csv(here::here("tests/exttestdata/102623 peaktable coculture simple.csv"), skip = 2)
  colnames(peak_df)[which(colnames(peak_df) %in% c("m/z"))] <- "mz"
  colnames(peak_df)[which(colnames(peak_df) %in% c("Retention time (min)"))] <- "rt"

  peak_df <- initialize_data(peak_df)
  
  sample_df <- readr::read_csv(here::here("tests/exttestdata/102623 samplelist.csv"))
  meta <- readr::read_csv(here::here("tests/exttestdata/102623 metadata simple.csv"))
  
  full_meta <- sample_df %>% left_join(meta, by = "Sample_Code") %>%
    filter(Biological_Group != "NA") %>%
    select(Injection, Sample_Code, Biological_Group) 
  
  peak_df_relfil <- check_mismatched_peaks(peak_df, ringwin = 0.5, isowin = 0.01, trwin = 0.005, max_iso_shift = 3, merge_peaks = TRUE)
  
  group_avgs <- filter_blank_3(peak_df_relfil, full_meta)
  error_prop <- read.csv(here::here("tests/exttestdata/102623 peaktable coculture simple_groupaverages.csv"), header = TRUE)
  colnames(error_prop) <- c("Compound", "mz", "rt", "biologicalGroup", "average")
  error_prop$Compound <- as.character(error_prop$Compound)
  error_prop <- error_prop[order(error_prop$Compound, error_prop$biologicalGroup), ]
  
  expect_true(all(group_avgs$Biological_Group %in% error_prop$biologicalGroup))
  
  group_avgs <- group_avgs[order(group_avgs$Compound, group_avgs$Biological_Group), ]
  expect_true(all(round(group_avgs$average, digits = 5) == round(error_prop$average, digits = 5)))
})

test_that("blank_filter_4 works correctly", {
  # the blank filter is passed merged peaks from filter 1 (mistmatched peaks)
  peak_df <- readr::read_csv(here::here("tests/exttestdata/102623 peaktable coculture simple.csv"), skip = 2)
  colnames(peak_df)[which(colnames(peak_df) %in% c("m/z"))] <- "mz"
  colnames(peak_df)[which(colnames(peak_df) %in% c("Retention time (min)"))] <- "rt"

  peak_df <- initialize_data(peak_df)
  
  sample_df <- readr::read_csv(here::here("tests/exttestdata/102623 samplelist.csv"))
  meta <- readr::read_csv(here::here("tests/exttestdata/102623 metadata simple.csv"))
  
  full_meta <- sample_df %>% left_join(meta, by = "Sample_Code") %>%
    filter(Biological_Group != "NA") %>%
    select(Injection, Sample_Code, Biological_Group) 
  
  peak_df_relfil <- check_mismatched_peaks(peak_df, ringwin = 0.5, isowin = 0.01, trwin = 0.005, max_iso_shift = 3, merge_peaks = TRUE)
  
  group_avgs <- filter_blank_4(data.table(peak_df_relfil), full_meta)
  error_prop <- read.csv(here::here("tests/exttestdata/102623 peaktable coculture simple_groupaverages.csv"), header = TRUE)
  colnames(error_prop) <- c("Compound", "mz", "rt", "biologicalGroup", "average")
  error_prop$Compound <- as.character(error_prop$Compound)
  error_prop <- error_prop[order(error_prop$Compound, error_prop$biologicalGroup), ]
  
  expect_true(all(group_avgs$Biological_Group %in% error_prop$biologicalGroup))
  
  # group_avgs <- group_avgs[order(group_avgs$Compound, group_avgs$Biological_Group), ]
  expect_true(all(round(group_avgs$average, digits = 5) == round(error_prop$average, digits = 5)))
})

test_that("parse_ion_list return the correct compounds for each group", {
  peak_df <- readr::read_csv(here::here("tests/exttestdata/102623 peaktable coculture simple.csv"), skip = 2)
  colnames(peak_df)[which(colnames(peak_df) %in% c("m/z"))] <- "mz"
  colnames(peak_df)[which(colnames(peak_df) %in% c("Retention time (min)"))] <- "rt"

  peak_df <- initialize_data(peak_df)
  
  sample_df <- readr::read_csv(here::here("tests/exttestdata/102623 samplelist.csv"))
  meta <- readr::read_csv(here::here("tests/exttestdata/102623 metadata simple.csv"))
  
  full_meta <- sample_df %>% left_join(meta, by = "Sample_Code") %>%
    filter(Biological_Group != "NA") %>%
    select(Injection, Sample_Code, Biological_Group) 
  
  peak_df_relfil <- check_mismatched_peaks(peak_df, ringwin = 0.5, isowin = 0.01, trwin = 0.005, max_iso_shift = 3, merge_peaks = TRUE)
  
  group_avgs <- filter_blank(peak_df_relfil, full_meta)
  group_filter_list <- parse_ions_by_group(group_avgs, group_threshold = 0.01)

  ang_18 <- read.csv(here::here("tests/exttestdata/output_ANG18 monoculture.csv"), header = FALSE)
  angdt <- read.csv(here::here("tests/exttestdata/output_ANGDT monoculture"), header = FALSE)
  blanks <- read.csv(here::here("tests/exttestdata/output_Blanks"), header = FALSE)  
  coculture <- read.csv(here::here("tests/exttestdata/output_Coculture"), header = FALSE)
  jc1 <- read.csv(here::here("tests/exttestdata/output_JC1 monoculture"), header = FALSE)
  jc28 <- read.csv(here::here("tests/exttestdata/output_JC28 monoculture"), header = FALSE)


  expect_true(all(group_filter_list$`ANG18 monoculture` == as.character(ang_18$V1)))
  expect_true(all(group_filter_list$`ANGDT monoculture` == as.character(angdt$V1)))
  expect_true(all(group_filter_list$`Blanks` == as.character(blanks$V1)))
  expect_true(all(group_filter_list$`Coculture` == as.character(coculture$V1)))
  expect_true(all(group_filter_list$`JC28 monoculture` == as.character(jc28$V1)))
  expect_true(all(group_filter_list$`JC28 monoculture` == as.character(jc28$V1)))
})

test_that("Group filter filters out groups correctly",
{
  peak_df <- readr::read_csv(here::here("tests/exttestdata/102623 peaktable coculture simple.csv"), skip = 2)
  colnames(peak_df)[which(colnames(peak_df) %in% c("m/z"))] <- "mz"
  colnames(peak_df)[which(colnames(peak_df) %in% c("Retention time (min)"))] <- "rt"

  peak_df <- initialize_data(peak_df)
  
  sample_df <- readr::read_csv(here::here("tests/exttestdata/102623 samplelist.csv"))
  meta <- readr::read_csv(here::here("tests/exttestdata/102623 metadata simple.csv"))
  
  full_meta <- sample_df %>% left_join(meta, by = "Sample_Code") %>%
    filter(Biological_Group != "NA") %>%
    select(Injection, Sample_Code, Biological_Group) 
  
  peak_df_relfil <- check_mismatched_peaks(peak_df, ringwin = 0.5, isowin = 0.01, trwin = 0.005, max_iso_shift = 3, merge_peaks = TRUE)
  
  group_avgs <- filter_blank(peak_df_relfil, full_meta)
  group_filter_list <- parse_ions_by_group(group_avgs, group_threshold = 0.01)
  peak_df_filtered <- apply_group_filter(peak_df_relfil, group_filter_list, "Blanks", remove_ions = TRUE)

  expect_true(all(!(group_filter_list$Blanks %in% peak_df_filtered$Compound)))
  

  peak_df_filtered <- apply_group_filter(peak_df_relfil, group_filter_list, "Blanks", remove_ions = FALSE)
  expect_true(all(peak_df_filtered == peak_df_relfil))
})

############################################
####        filter 3: CV filter         ####
############################################

test_that("cv_filter return the correct number of ions failing cv filter", {
  peak_df <- readr::read_csv(here::here("tests/exttestdata/102623 peaktable coculture simple.csv"), skip = 2)
  colnames(peak_df)[which(colnames(peak_df) %in% c("m/z"))] <- "mz"
  colnames(peak_df)[which(colnames(peak_df) %in% c("Retention time (min)"))] <- "rt"

  peak_df <- initialize_data(peak_df)
  
  sample_df <- readr::read_csv(here::here("tests/exttestdata/102623 samplelist.csv"))
  meta <- readr::read_csv(here::here("tests/exttestdata/102623 metadata simple.csv"))
  
  full_meta <- sample_df %>% left_join(meta, by = "Sample_Code") %>%
    filter(Biological_Group != "NA") %>%
    select(Injection, Sample_Code, Biological_Group) 
  
  peak_df_relfil <- check_mismatched_peaks(peak_df, ringwin = 0.5, isowin = 0.01, trwin = 0.005, max_iso_shift = 3, merge_peaks = TRUE)
  
  group_avgs <- filter_blank(peak_df_relfil, full_meta)
  group_filter_list <- parse_ions_by_group(group_avgs, group_threshold = 0.01)
  peak_df_filtered <- apply_group_filter(peak_df_relfil, group_filter_list, "Blanks", remove_ions = TRUE)
  ions_failing_cv_mean <- cv_filter(peak_df_filtered, full_meta, cv_threshold = 0.2, cv_param = c("mean"))
  ions_failing_cv_median <- cv_filter(peak_df_filtered, full_meta, cv_threshold = 0.2, cv_param = c("median"))

  expect_equal(length(ions_failing_cv_median), 61)
  expect_equal(length(ions_failing_cv_mean), 86)
  
})