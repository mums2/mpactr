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
# full_meta <- sample_df %>% left_join(meta, by = "Sample_Code") %>%
#  filter(Biological_Group != "NA") %>%
#  select(Injection, Sample_Code, Biological_Group) 

# data_frame <- peak_df_relfil
# metadata <- full_meta

filter_blank <- function(data_frame, metadata) {
  
  # extract feature sample (rows are features, columns are samples)
  # data_frame <- peak_df_relfil
  # metadata <- full_meta
  ft <- as.data.frame(data_frame[ , !(colnames(data_frame) %in% c("mz", "rt", "kmd"))])
  rownames(ft) <- as.character(ft$Compound)
  ft <- ft[ , !(colnames(ft) %in% c("Compound"))]
  
  # transpose so rows are samples and columns are features
  ft_t <- as.data.frame(t(ft)) 
  # ft_t <- tibble::rownames_to_column(ft_t, var = "filename")
  
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

# 05/15/2024
filter_blank_2 <- function(data_frame, metadata) {

  # extract feature sample (rows are features, columns are samples)
  data_frame <- peak_df_relfil
  metadata <- full_meta
  ft <- as.data.frame(data_frame[ , !(colnames(data_frame) %in% c("mz", "rt", "kmd"))])
  rownames(ft) <- as.character(ft$Compound)
  ft <- ft[ , !(colnames(ft) %in% c("Compound"))]
  row_count <- nrow(metadata)
  injections <- metadata$Injection
  injection_sample_code_dict <- vector(mode ="list", length = row_count)
  names(injection_sample_code_dict) <- injections
  injection_biol_group_dict <- vector(mode ="list", length = row_count)
  names(injection_biol_group_dict) <- injections
  for(index in 1:row_count)
  {
    injection_sample_code_dict[[metadata$Injection[[index]]]] <- metadata$Sample_Code[[index]]
    injection_biol_group_dict[[metadata$Injection[[index]]]] <- metadata$Biological_Group[[index]]
  }
  sample_codes <- unique(metadata$Sample_Code)
  column_count <- ncol(ft)
  list_of_rows_for_computation <- vector(mode ="list", length = length(sample_codes))
  names(list_of_rows_for_computation) <- sample_codes
  for(i in 1:column_count)
  {
    injection <- colnames(ft)[i]
    sample_code <- injection_sample_code_dict[[injection]]
    list_of_rows_for_computation[[sample_code]] <- c(list_of_rows_for_computation[[sample_code]], i)
  }
  average <- c()
  biol_RSD <- c()
  biol_n <- c()
  for(i in 1:nrow(ft))
  {
    column_name <-
    ft[i,]
  }


}

# data_frame <- peak_df_relfil
# metadata <- full_meta
filter_blank_3 <- function(data_frame, metadata) {
  # extract feature sample (rows are features, columns are samples)
  # data_frame <- peak_df_relfil
  # metadata <- full_meta
  ft <- as.data.frame(data_frame[ , !(colnames(data_frame) %in% c("mz", "rt", "kmd"))])
  rownames(ft) <- as.character(ft$Compound)
  ft <- ft[ , !(colnames(ft) %in% c("Compound"))]
  
  # transpose so rows are samples and columns are features
  ft_t <- as.data.frame(t(ft)) 
  # ft_t <- tibble::rownames_to_column(ft_t, var = "filename")
  
  # make sure metadata rows in the same order as ft_t rownames
  metadata <- metadata[order(match(metadata$Injection, rownames(ft_t))), ]
  # metadata$Injection <- as.factor(metadata$Injection)
  # metadata$Sample_Code <- as.factor(metadata$Sample_Code)
  # metadata$Biological_Group <- as.factor(metadata$Biological_Group)
  
  groups_stats_list <- vector(mode="list", length = length(unique(metadata$Biological_Group)))
  names(groups_stats_list) <- unique(metadata$Biological_Group)
  lapply(names(groups_stats_list), 
         function(x) {
           dat <- ft_t[which(rownames(ft_t) %in% 
                                metadata$Injection[which(metadata$Biological_Group == x)]), ]

          tech_lvls <- as.factor(metadata$Sample_Code[which(rownames(dat) %in% metadata$Injection)])
          
          groups_stats_list[[x]]$Compound <<- colnames(dat)
          groups_stats_list[[x]]$Biological_Group <<- c(rep(x, ncol(dat)))
          groups_stats_list[[x]]$average <<- unlist(lapply(dat, mean))
          groups_stats_list[[x]]$BiolRSD <<- unlist(lapply(dat, rsd))
          groups_stats_list[[x]]$bioln <<- unlist(lapply(dat, length))
          groups_stats_list[[x]]$techRSD <<- unlist(lapply(lapply(dat, 
                function(x) {
                  by(x, tech_lvls, rsd)}),
            mean))
          groups_stats_list[[x]]$techn <<- unlist(lapply(lapply(dat, 
                function(x) {
                  by(x, tech_lvls, length)}),
            mean))
        
      }
    )

  group_stats <- data.table::rbindlist(lapply(groups_stats_list, as.data.frame))
  group_stats <- group_stats[order(group_stats$Compound, decreasing = FALSE), ]

  return(group_stats)
}
# ## allison - testing a combination of apply funs and by to calculate group statistics
# rsd <- function(values) {
#   ifelse(mean(values) != 0, (sd(values) / mean(values)), NA_real_)
# }

# groups_stats_list <- vector(mode="list", length = length(unique(metadata$Biological_Group)))
# names(groups_stats_list) <- unique(metadata$Biological_Group)
# lapply(names(groups_stats_list), 
#       function(x) {
#         dat <- ft_t[which(rownames(ft_t) %in% 
#                                 metadata$Injection[which(metadata$Biological_Group == x)]), ]

#         tech_lvls <- as.factor(metadata$Sample_Code[which(rownames(dat) == metadata$Injection)])
        
#         groups_stats_list[[x]]$Compound <<- colnames(dat)
#         groups_stats_list[[x]]$Biological_Group <<- c(rep(x, ncol(dat)))
#         groups_stats_list[[x]]$average <<- unlist(lapply(dat, mean))
#         groups_stats_list[[x]]$BiolRSD <<- unlist(lapply(dat, rsd))
#         groups_stats_list[[x]]$bioln <<- unlist(lapply(dat, length))
#         groups_stats_list[[x]]$techRSD <<- unlist(lapply(lapply(dat, 
#               function(x) {
#                 by(x, tech_lvls, rsd)}),
#           mean))
#         groups_stats_list[[x]]$techn <<- unlist(lapply(lapply(dat, 
#               function(x) {
#                 by(x, tech_lvls, length)}),
#           mean))
        
        
        
#       }
#     )

# group_stats <- data.table::rbindlist(lapply(groups_stats_list, as.data.frame))
# group_stats <- group_stats[order(group_stats$Compound, decreasing = FALSE)]



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

apply_group_filter  <- function(data_frame, group_filter_list, group, remove_ions = TRUE) {
  if(isFALSE(remove_ions)) {
    return(data_frame)
  }

  ions <- group_filter_list[[group]]
  data_frame <- data_frame[!(data_frame$Compound %in% ions), ]
  return(data_frame)

}



### CV filter
cv_filter <- function(data_frame, metadata, cv_threshold, cv_param) {
    # if mismatch and merge and group filters (with remove_ions == TRUE) 
    # are selected, this is expecting the filtered df as input

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
  
  cv <- ft_t %>%
    pivot_longer(cols = c(as.character(data_frame$Compound[1]):as.character(tail(data_frame$Compound, 1))), names_to = "Compound") %>%
    group_by(Compound, Sample_Code, Biological_Group) %>%
    summarise(cv = if_else(mean(value) != 0, sd(value) / mean(value), NA_real_)) %>%
    ungroup() %>%
    group_by(Compound) %>%
    summarise(mean_cv = mean(cv, na.rm = TRUE),
              median_cv = median(cv, na.rm = TRUE))

  if (cv_param == "mean") {
    fail_cv <- cv$Compound[!(cv$mean_cv < cv_threshold)]
    return(fail_cv)
  }

  if (cv_param == "median") {
    fail_cv <- cv$Compound[!(cv$median_cv < cv_threshold)]
    return(fail_cv)
  }
}


#insource ion filter
identify_insource_ions <- function(data_frame) {
  data_frame <- peak_df_filtered
  rt_list <- unique(data_frame$rt)
  
  singleslist <- c()
  insourcelist <- c()
  clusterlist <- c()
  clusterdfs <- c()
  filtered_df <- data_frame[ , colnames(data_frame) %in% c("Compound", "mz", "rt")] #msdata_ind
  ftrgrps <- list()
  mergegroups <- list()
  for (rt in rt_list) {
    filtered_df <- data_frame[data_frame$rt == rt , colnames(data_frame) %in% c("Compound", "mz", "rt")]

    if (nrow(filtered_df) == 1) {
      singleslist <- c(singleslist, filtered_df$Compound)
    }
    else{ 
      filtered_df2 <- data_frame[data_frame$rt == rt , !(colnames(data_frame) %in% c("Compound", "mz", "rt", "kmd"))]
      filtered_df2<- as.data.frame(t(filtered_df2))
      colnames(filtered_df2) <- filtered_df$Compound
      corr <- stats::cor(filtered_df2, method = c("pearson")) # why correlate 
      inverse_unique_corr <-  1 - (unname(corr)[ , 1][which(unname(corr)[ , 1] != 1)])
      cluster <- hclust(inverse_unique_corr, method = "complete") # why cluster the correlated intensities and not just the intensities 

    }
  }
}
# R
#            12       761
# 12  1.0000000 0.7899264
# 761 0.7899264 1.0000000
# mpact
#            12       761
# 12  1.000000 0.751101
# 761 0.751101 1.000000

# 0, 0.248899

test_df2 <- data_frame[which(data_frame$rt == data_frame[data_frame$Compound == "12", "rt"] ), !(colnames(data_frame) %in% c("Compound", "mz", "rt", "kmd"))]
test_df2 <- as.data.frame(t(test_df2))
colnames(test_df2) <- test_df2$Compound
corr <- stats::cor(test_df2, method = c("pearson"))
x <- 
inverse_unique_corr <-  1 - (unname(corr)[ , 1][which(unname(corr)[ , 1] != 1)])
cluster <- hclust(inverse_unique_corr, method = "complete")

hclust(test_df2, method = "complete")



# def decon(analysis_params, ionfilters):
#     """
#     Perform peak deconvolution on the data and update the ion filters and dictionary.
    
#     This function reads in formatted data and performs peak deconvolution using the correlation clustering method. It
#     then groups peaks into clusters, generates a list of precursor and fragment ions, and updates the ion filters and
#     dictionary. Finally, the function writes the formatted data back out to disk.
    
#     Parameters:
#         analysis_params (object): Analysis parameters.
#         ionfilters (dict): Dictionary of ion filters.
    
#     Returns:
#         dict: Updated dictionary of ion filters.
#     """
#     # Read in data
#     # msdata_ind = pd.read_csv(analysis_params.outputdir / (analysis_params.filename.stem + '_formatted.csv'),
#     #                          sep=',', header=[2], index_col=[0, 1]).iloc[:, :1]
#     # rtlist = msdata_ind['Retention time (min)'].unique()
#     # msdata_ind = pd.read_csv(analysis_params.outputdir / (analysis_params.filename.stem + '_formatted.csv'),
#     #                          sep=',', header=[2], index_col=[0]).iloc[:, :2]
#     # msdata = pd.read_csv(analysis_params.outputdir / (analysis_params.filename.stem + '_formatted.csv'),
#     #                      sep=',', header=[2], index_col=[0]).iloc[:, 2:]

#     # Initialize variables
#     # singleslist = []
#     # insourcelist = []
#     # clusterlist = []
#     # clusterdfs = []
#     # filtereddf = msdata_ind
#     # ftrgrps = {}
#     # mergegroups = {}

#     # Loop through retention times
#     for elem in rtlist:
#         # filtereddf = msdata_ind.loc[msdata_ind.iloc[:, 1] == elem, :] # Return a list for duplicates

#         # if len(filtereddf.index) == 1:
#         #     singleslist.append(filtereddf.index.to_list()[0])
#         else:
#             # Perform deconvolution
#             # unseplist = filtereddf.index.tolist()
#             # filtereddf2 = msdata.loc[unseplist, :].transpose()
#             # corr = filtereddf2.corr()
#             # corr2 = reformatcorr(corr)



    #         linkage = spc.linkage(corr2, method='complete')
    #         np.clip(linkage, 0, None, linkage)
    #         idx = spc.fcluster(linkage, 1 - analysis_params.deconthresh, 'distance')

    #         # Group peaks by cluster
    #         decongroups = {}
    #         for group in range(1, max(idx) + 1):
    #             decongroups[group] = []
    #         for peak in range(0, idx.shape[0]):
    #             if peak < len(unseplist):
    #                 decongroups[idx[peak]].append(unseplist[peak])

    #         # Append groups to appropriate lists
    #         for group in decongroups:
    #             if len(decongroups[group]) == 1:
    #                 singleslist.append(decongroups[group][0])
    #             else:
    #                 clusterlist.append(decongroups[group])
    #                 tempdf = msdata_ind.loc[decongroups[group]].sort_values(by=['m/z'], ascending=False)
    #                 if len(tempdf.index.to_list()) > 0:
    #                     precursor = tempdf.index.to_list()[0]
    #                     fragments = tempdf.index.to_list()[1:]
    #                     singleslist.append(precursor)
    #                     insourcelist += fragments
    #                     clusterdfs.append(tempdf)
    #                     mergegroups[precursor] = ionmerge(precursor, fragments)


    # # Write out data
    # msdata = pd.read_csv(analysis_params.outputdir / (analysis_params.filename.stem + '_formatted.csv'),
    #                      sep=',', header=[0, 1, 2], index_col=[0])
    # msdata = msdata.loc[singleslist]
    # msdata.to_csv(analysis_params.outputdir / (analysis_params.filename.stem + '_formatted.csv'),
    #               header=True, index=True)

    # # Update ion filters
    # ionfilters['insource'] = ionfilter('', insourcelist)
    # ionfilters['insource'].merge = mergegroups

    # # Update ion dictionary
    # iondict = pd.read_csv(analysis_params.outputdir / 'iondict.csv', sep=',', header=[0], index_col=0)
    # iondict['pass_insource'] = ~iondict.index.isin(ionfilters['insource'].ions)
    # iondict.to_csv(analysis_params.outputdir / 'iondict.csv', header=True, index=True)
    
    # # Convert ion filters to MSP format passes and prints errors
    # try:
    #     mspwriter.convert_to_msp(ionfilters['insource'], analysis_params)
    # except Exception:
    #     print('Error in MSP writer')
    #     pass
    
    # return ionfilters
