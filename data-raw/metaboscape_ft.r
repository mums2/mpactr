data_metaboscape <- read_csv(here::here("tests/testthat/exttestdata/MJB_MonoVSCoculture_metaboscape_ft.csv"))

meta <- data.frame(Injection = colnames(data_metaboscape)) %>%
    filter(startsWith(Injection, "UM") ) %>%
    filter(Injection == str_detect(Injection, "Media") | str_detect(Injection, "Coculture") | str_detect(Injection, "ANGDT")) %>%
    mutate(Sample_Code = str_split_i(Injection, "_", 2),
           Biological_Group = case_when(str_detect(Injection, "ANGDT") ~ "ANGDT",
                             str_detect(Injection, "Coculture") ~ "Coculture",
                             str_detect(Injection, "Media") ~ "Media",
                             TRUE ~ NA_character_)
        )
    
write_csv(meta, here::here("inst/extdata/cultures_metaboscape_metadata.csv"))

data_metaboscape %>%
    select(FEATURE_ID, RT, PEPMASS, CCS, SIGMA_SCORE, ADDUCT, KEGG, CAS, MaxIntensity, all_of(meta$Injection)) %>%
    sample_n(size = 1000) %>% # randomly subset the rows to reduce the size of the example df
    write_csv(here::here("inst/extdata/cultures_metaboscape_peaktable.csv"))
    

