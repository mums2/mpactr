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
    #data_frame_merged <- data_frame[which(!(data_frame$Compound %in% cut_ions)), ]
    #return(as.data.frame(data_frame_merged[order(data_frame_merged$Compound, decreasing = FALSE), ]))
    return(merge_ions(data_frame, ion_filter_list))
  }
  return(ion_filter_list)

}

merge_ions <- function(data_frame, ion_filter_list) {
  ion_dat <- data_frame %>%
    select(Compound, mz, rt, kmd) %>%
    mutate(Compound = as.character(Compound))
  
  dat <- data_frame %>%
    select(-mz, -rt, -kmd) %>%
    tibble::column_to_rownames(var = "Compound") %>%
    t() %>%
    as.data.frame()

  for (ion in names(ion_filter_list$merge_groups)) {
    dat[, ion] <- apply(dat[, c(ion,
                                ion_filter_list$merge_groups[[ion]])], 1, sum)
    dat[, as.character(ion_filter_list$merge_groups[[ion]])] <- NULL
  }
 
  data_frame <- dat %>%
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "Compound") %>%
    left_join(ion_dat, by = "Compound")

  return(data_frame)
  }

# peak_table <- peak_df %>%
#   select(-mz, -rt, -kmd) %>%
#   column_to_rownames(var = "Compound") %>%
#   t() %>%
#   as.data.frame()
  
# p2 <- peak_table

# for (ion in names(ion_filter_list$merge_groups)) {
#   p2[ , ion]  <- apply(p2[ , c(ion, relfil_ion_list$merge_groups[["1188"]])], 1, sum)
#   p2[, as.character(relfil_ion_list$merge_groups[["1188"]])] <- NULL
# }




# blank filter
# the current expected metadata format as input
#full_meta <- sample_df %>% left_join(meta, by = "Sample_Code") %>%
#  filter(Biological_Group != "NA") %>%
#  select(Injection, Sample_Code, Biological_Group) 

filter_blank <- function(data_frame, metadata) {
  
  # extract feature sample (rows are features, columns are samples)
  ft <- as.data.frame(data_frame[ , !(colnames(data_frame) %in% c("mz", "rt", "kmd"))])
  rownames(ft) <- as.character(ft$Compound)
  ft <- ft[ , !(colnames(ft) %in% c("Compound"))]
  
  # transpose so rows are samples and columns are features
  ft_t <- as.data.frame(t(ft))
  
  # make sure metadata rows in the same order as ft_t rownames
  metadata <- metadata[order(match(metadata$Injection, rownames(ft_t))), ]
  metadata$Injection <- as.factor(metadata$Injection)
  metadata$Sample_Code <- as.factor(metadata$Sample_Code)
  metadata$Biological_Group <- as.factor(metadata$Biological_Group)
  
  ft_t <- cbind(ft_t, metadata[, c("Sample_Code", "Biological_Group")])
  
  # calculate RSD for biological and technical groups
  biol_stats <- ft_t %>% 
    pivot_longer(cols = c(as.character(data_frame$Compound[1]):as.character(tail(data_frame$Compound, 1))), names_to = "Compound") %>%
    group_by(Compound, Biological_Group) %>%
    summarise(average = mean(value),
              biolRSD = if_else(average != 0, sd(value) / mean(value), 0),
              bioln = n()) %>%
    ungroup()
  
  tech_stats <- ft_t %>%
    pivot_longer(cols = c(as.character(data_frame$Compound[1]):as.character(tail(data_frame$Compound, 1))), names_to = "Compound") %>%
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
#                bioln = n()) %>%
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

# parse ions by group
parse_ions_by_group <- function(group_stats, group_threshold = 0.01) {
#  groups <- unique(group_stats$Biological_Group)
  group_stats <- group_avgs
  
  group_stats_tbl <- group_stats %>%
    mutate(Compound = as.numeric(Compound)) %>%
    arrange(Compound, desc = TRUE) %>%
    select(Compound, Biological_Group, average) %>%
    pivot_wider(names_from = Biological_Group, values_from = average) %>%
    column_to_rownames(var = "Compound") %>% # extract feature table where samples are columns 
    t() %>%
    as.data.frame()
    
  compound_max <- apply(group_stats_tbl, 2, max)
  
  # normalize group averages for each compound 
  group_stats_norm <- sapply(seq_along(1:ncol(group_stats_tbl)), 
         \(x){
           group_stats_tbl[ , x] <- group_stats_tbl[ , x] / compound_max[x]
         })
  group_stats_norm <- as.data.frame(group_stats_norm)
  # reset row and column names
  colnames(group_stats_norm) <- colnames(group_stats_tbl)
  rownames(group_stats_norm) <- rownames(group_stats_tbl)
  
  # need to pull filtered compounds by group
  biol_groups <- split(group_stats_norm, 1:nrow(group_stats_norm))
  names(biol_groups) <- rownames(group_stats_norm)
  
  group_filter_list <- lapply(biol_groups, \(x){names(x)[which(x > group_threshold)]})
    
}

# group_filter_list$`ANG18 monoculture` %in% group_stats$Compound

#lapply(biol_groups, \(x){which(x > group_threshold)})

#group_filter_list$`ANG18 monoculture`
#unique(unlist(group_filter_list))

#lapply(group_filter_list, length)
#712 + 286 + 796 + 370 + 534 + 249
#length(unique(unlist(group_filter_list)))


# def groupfilter(group, msdata_blankfilterg, analysis_params):
#     """
#     Filters ions in msdata_blankfilterg based on presence in group.

#     Parameters:
#         group (str): Group to filter by.
#         msdata_blankfilterg (pandas.DataFrame): Dataframe to filter.
#         analysis_params (object): Analysis parameters.

#     Returns:
#         list: List of filtered ions.
#     """
#     # Filters ions based on group presence > 100 ions/abundance
#     msdata_blankfilterg = msdata_blankfilterg[msdata_blankfilterg[group] > analysis_params.blankfilthresh]
#     msdata_blankfilterg = msdata_blankfilterg.reset_index()
#     # Returns filtered dataframe index as list
#     groupfilterlist = msdata_blankfilterg.iloc[0:,0].tolist()
#     return groupfilterlist
      
# def parsionlists(analysis_params):
#    """
#    Parses lists of ions present in each group.
#
#    Parameters:
#        analysis_params (object): Analysis parameters.
#
#    Returns:
#        dict: Dictionary of group ion lists.
#    """
#    msdata_blankfilterg = pd.read_csv(analysis_params.outputdir / (analysis_params.filename.stem + '_groupaverages.csv'), sep = ',', header = [0], index_col = [0, 1, 2, 3]).unstack()
#    msdata_blankfilterg.columns = msdata_blankfilterg.columns.droplevel()
#    msdata_blankfilterg['max'] = msdata_blankfilterg.max(numeric_only=True, axis=1)
#    for column in msdata_blankfilterg:
#        msdata_blankfilterg[column] = msdata_blankfilterg[column] / msdata_blankfilterg['max']
#    msdata_blankfilterg = msdata_blankfilterg.drop(columns=['max'])
#    biolgroups = msdata_blankfilterg.columns.tolist()
    # groupionlists = {}
    # for group in biolgroups:
    #     groupionlists[group] = groupfilter(group, msdata_blankfilterg, analysis_params)
    # return groupionlists

# if analysis_params.grpave:
#        stats.groupave(analysis_params)
#        print('Parsing ion lists')



#         groupionlists = filter.parsionlists(analysis_params)
# if analysis_params.blnkfltr:
#         msdata = pd.read_csv(analysis_params.outputdir / (analysis_params.filename.stem + '_formatted.csv'), sep=',', header=[0, 1, 2], index_col=None)
#         msdata = filter.listfilter(msdata, groupionlists[analysis_params.blnkgrp], False)
#         msdata = msdata.drop(analysis_params.blnkgrp, axis=1, level=0)
#         msdata.to_csv(analysis_params.outputdir / (analysis_params.filename.stem + '_formatted.csv'), header=True, index=False)
#         iondict = pd.read_csv(analysis_params.outputdir / 'iondict.csv', sep=',', header=[0], index_col=0)
#         iondict['pass_blnkfil'] = ~iondict.index.isin(groupionlists[analysis_params.blnkgrp])
#         iondict.to_csv(analysis_params.outputdir / 'iondict.csv', header=True, index=True)