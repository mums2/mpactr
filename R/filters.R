# TODO Generate a class-like structure for the tree-plot



#' @params data_frame with first three columns: Compound, mz, rt
solvent_blank_filter <- function(data_frame, sample_to_filter)
{
  result <- apply(data_frame[ , sample_to_filter, drop = FALSE], 1, function(x) {sum(x) <= 0})
  return(data_frame[result, !names(data_frame) %in% sample_to_filter])
}



# add params: df
check_mismatched_peaks <- function(data_frame, ringwin, isowin, trwin, max_iso_shift, merge_peaks) {

ion_filter_list <- list()
current_rt <- 0
current_mass <- 0
mass_diff <- 0
kmd_diff <- 0
cut_ions <- c() # list
merge_groups <- list() # dictonary

  data_frame <- data_frame[order(data_frame$mz, decreasing = FALSE), ]
  number_of_rows <- nrow(data_frame)
  for(i in seq_along(1:number_of_rows)) {
    current_feature <- data_frame[i, c("Compound", "mz", "rt")]
    current_rt <- current_feature$rt
    current_mass <- current_feature$mz
    current_ion <- current_feature$Compound
  if(current_ion == 1248)
        {
          print("stop")
        }
    if (!(current_ion %in% cut_ions))
    {
      for(j in (i + 1):number_of_rows) {
        if(j > number_of_rows)
        {
          break
        }
        if(data_frame$Compound[j] %in% cut_ions) {
          next
        }

        if(current_ion == 1248)
        {
          print("stop")
           if( data_frame$Compound[j] == 1250)
          {print("stop")}
        }
        mass_diff <- data_frame$mz[j] - current_mass
        kmd_diff <- mass_diff - floor(mass_diff)

        if (abs(mass_diff) > max_iso_shift - 0.4) { # BL  - why 0.4??
          break
        }

        rt_diff <- data_frame$rt[j] - current_rt
        ring_band <- floor(abs(mass_diff) * (1/ringwin)) %% (1/ringwin)
        double_band <- kmd_diff - .5004 -
        (floor(mass_diff) * .007) # BL - why 0.500 and 0.007?

        if((abs(rt_diff) <= trwin) && 
           (mass_diff <= max_iso_shift - 0.4) &&
           (ring_band == 0 ||
            double_band < isowin)) {
            cut_ions <- c(cut_ions, data_frame$Compound[j])
            merge_groups[[as.character(current_ion)]] <- c(merge_groups[[as.character(current_ion)]], data_frame$Compound[j])
        }
      }
    }
  }
  # TODO Look into removing merge groups
  ion_filter_list[["cut_ions"]] <- cut_ions
  ion_filter_list[["merge_groups"]] <- merge_groups
  
  if (isTRUE(merge_peaks)) {
    return(data_frame[which(!(data_frame$Compound %in% cut_ions)), ])
  }
  return(ion_filter_list)

}

# blank filter
# the current expected metadata format as input
#full_meta <- sample_df %>% left_join(meta, by = "Sample_Code") %>%
#  filter(Biological_Group != "NA") %>%
#  select(Injection, Sample_Code, Biological_Group) 

filter_blank <- function(data_frame, metadata) {
  
  # extract feature sample (rows are features, columns are samples)
  ft <- as.data.frame(peak_df[ , !(colnames(peak_df) %in% c("mz", "rt", "kmd"))])
  rownames(ft) <- as.character(ft$Compound)
  ft <- ft[ , !(colnames(ft) %in% c("Compound"))]
  
  # transpose so rows are samples and columns are features
  ft_t <- t(ft)
  
  # make sure metadata rows in the same order as ft_t rownames
  full_meta <- full_meta[order(match(full_meta$Injection, rownames(ft_t))), ]
  full_meta$Injection <- as.factor(full_meta$Injection)
  full_meta$Sample_Code <- as.factor(full_meta$Sample_Code)
  full_meta$Biological_Group <- as.factor(full_meta$Biological_Group)
  
  # calculate RSD for biological and technical groups
  biol_stats <- ft_t %>% 
    pivot_longer(cols = c(as.character(peak_df$Compound[1]):as.character(tail(peak_df$Compound, 1))), names_to = "Compound") %>%
    group_by(Compound, Biological_Group) %>%
    summarise(average = mean(value),
              biolRSD = if_else(average != 0, sd(value) / mean(value), 0),
              bioln = n()) %>%
    ungroup()
  
  tech_stats <- ft_t %>%
    pivot_longer(cols = c(as.character(peak_df$Compound[1]):as.character(tail(peak_df$Compound, 1))), names_to = "Compound") %>%
    group_by(Compound, Sample_Code, Biological_Group) %>%
    summarise(sd = if_else(mean(value) != 0, sd(value) / mean(value), 0),
              n = n()) %>%
    ungroup() %>%
    group_by(Compound, Biological_Group) %>%
    summarise(techRSD = mean(sd),
              techn = mean(n))
    
  group_stats <- biol_stats %>%
  left_join(tech_stats, by = c("Compound", "Biological_Group"))
  
  return(group_stats)
}

# tidy approach:
#biol_stats <- ft_t %>% 
#      pivot_longer(cols = c(as.character(peak_df$Compound[1]):as.character(tail(peak_df$Compound, 1))), names_to = "Compound") %>%
#      group_by(Compound, Biological_Group) %>%
#      summarise(average = mean(value),
#                biolRSD = if_else(average != 0, sd(value) / mean(value), 0),
                bioln = n()) %>%
#      ungroup()
#      
#tech_stats <- ft_t %>%
#  pivot_longer(cols = c(as.character(peak_df$Compound[1]):as.character(tail(peak_df$Compound, 1))), names_to = "Compound") %>%
#  group_by(Compound, Sample_Code, Biological_Group) %>%
#  summarise(sd = if_else(mean(value) != 0, sd(value) / mean(value), 0),
#            n = n()) %>%
#  ungroup() %>%
#  group_by(Compound, Biological_Group) %>%
#  summarise(techRSD = mean(sd),
#            techn = mean(n))
#  
#  group_stats <- biol_stats %>%
#  left_join(tech_stats, by = c("Compound", "Biological_Group"))


# Reshape approach - unstabe *in progress*
# stats::reshape(data = ft_t,
#                direction = "long",
#                varying = which(colnames(ft_t) == peak_df$Compound[1]):which(colnames(ft_t) == tail(peak_df$Compound, 1)),
#                v.names = colnames(ft_t)[which(colnames(ft_t) == peak_df$Compound[1]):which(colnames(ft_t) == tail(peak_df$Compound, 1))],
#                idvar = c("Injection", "Sample_Code", "Biological_Group")
# )

# stats::reshape(data = ft_t,
#                direction = "long",
#                varying = which(colnames(ft_t) == peak_df$Compound[1]):which(colnames(ft_t) == tail(peak_df$Compound, 1)),
#                v.names = colnames(ft_t)[which(colnames(ft_t) == peak_df$Compound[1]):which(colnames(ft_t) == tail(peak_df$Compound, 1))],
#                timevar = c("Injection", "Sample_Code", "Biological_Group"))

# ft_t_sub <- ft_t %>% select(`1`, `2`, `3`, Biological_Group, Sample_Code)

if analysis_params.grpave:
        stats.groupave(analysis_params)
        print('Parsing ion lists')
        groupionlists = filter.parsionlists(analysis_params)
if analysis_params.blnkfltr:
        msdata = pd.read_csv(analysis_params.outputdir / (analysis_params.filename.stem + '_formatted.csv'), sep=',', header=[0, 1, 2], index_col=None)
        msdata = filter.listfilter(msdata, groupionlists[analysis_params.blnkgrp], False)
        msdata = msdata.drop(analysis_params.blnkgrp, axis=1, level=0)
        msdata.to_csv(analysis_params.outputdir / (analysis_params.filename.stem + '_formatted.csv'), header=True, index=False)
        iondict = pd.read_csv(analysis_params.outputdir / 'iondict.csv', sep=',', header=[0], index_col=0)
        iondict['pass_blnkfil'] = ~iondict.index.isin(groupionlists[analysis_params.blnkgrp])
        iondict.to_csv(analysis_params.outputdir / 'iondict.csv', header=True, index=True)