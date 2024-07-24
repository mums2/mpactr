# Code to generate metadata.csv
# metadata.csv is a file of sample metadata in mpactR format, for the data provided in cultures_peak_table.csv.

library(tidyverse)

ft <- read_csv(here::here("tests/testthat/exttestdata/102623_peaktable_MonoVsCoculture_expandedgradient.csv"), skip = 2)

data.frame(Injection = colnames(ft)) %>%
  filter(startsWith(Injection, "10")) %>%
  mutate(
    Sample_Code = str_split_i(Injection, "_", 2),
    Biological_Group = case_when(
      str_detect(Injection, "ANGDT") ~ "ANGDT",
      str_detect(Injection, "ANG18") ~ "ANG18",
      str_detect(Injection, "JC1") ~ "JC1",
      str_detect(Injection, "JC28") ~ "JC28",
      str_detect(Injection, "Coculture") ~ "Coculture",
      str_detect(Injection, "Media") ~ "Media",
      str_detect(Injection, "Mixed") ~ "Mixed_Monoculture",
      str_detect(Injection, "Blank") ~ "Solvent_Blank",
      TRUE ~ NA_character_
    ),
    dilution = if_else(str_detect(Injection, "0.25"), "0.25", "1")
  ) %>%
  filter(dilution == 1) %>%
  write_csv(here::here("inst/extdata/cultures_metadata.csv"))


# Filter the feature table by diltion for a simple example dataset
colnames(ft)
samples <- read_csv(here::here("inst/extdata/cultures_metadata.csv")) %>% pull(sample_id)

pt <- ft %>%
  select(Compound, `m/z`, `Retention time (min)`, all_of(samples)) %>%
  rbind(rep(NA, ncol(.)), rep(NA, ncol(.)), colnames(.), .)

write_csv(pt, here::here("inst/extdata/cultures_peak_table.csv"), col_names = FALSE)
