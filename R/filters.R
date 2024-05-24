# TODO Generate a class-like structure for the tree-plot



#' @params data_table with first three columns: Compound, mz, rt
solvent_blank_filter <- function(data_table, sample_to_filter)
{
  result <- apply(data_table[ , sample_to_filter, drop = FALSE], 1, function(x) {sum(x) <= 0})
  return(data_table[result, !names(data_table) %in% sample_to_filter])
}



# add params: df
#' @param data_table a data.table 
check_mismatched_peaks <- function(data_table, ringwin, isowin, trwin, max_iso_shift, merge_peaks) {
ion_filter_list <- list()
current_rt <- 0
current_mass <- 0
mass_diff <- 0
kmd_diff <- 0
cut_ions <- c() # list
merge_groups <- list() # dictonary

  data_table <- data_table[order(data_table$mz, decreasing = FALSE), ]
  number_of_rows <- nrow(data_table)
  for(i in seq_along(1:number_of_rows)) {
    current_feature <- data_table[i, c("Compound", "mz", "rt")]
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
        if(data_table$Compound[j] %in% cut_ions) {
          next
        }

        if(current_ion == 1248)
        {
          print("stop")
           if( data_table$Compound[j] == 1250)
          {print("stop")}
        }
        mass_diff <- data_table$mz[j] - current_mass
        kmd_diff <- mass_diff - floor(mass_diff)

        if (abs(mass_diff) > max_iso_shift - 0.4) { # BL  - why 0.4??
          break
        }

        rt_diff <- data_table$rt[j] - current_rt
        ring_band <- floor(abs(mass_diff) * (1/ringwin)) %% (1/ringwin)
        double_band <- kmd_diff - .5004 -
        (floor(mass_diff) * .007) # BL - why 0.500 and 0.007?

        if((abs(rt_diff) <= trwin) && 
           (mass_diff <= max_iso_shift - 0.4) &&
           (ring_band == 0 ||
            double_band < isowin)) {
            cut_ions <- c(cut_ions, data_table$Compound[j])
            merge_groups[[as.character(current_ion)]] <- c(merge_groups[[as.character(current_ion)]], data_table$Compound[j])
        }
      }
    }
  }
  # TODO Look into removing merge groups
  ion_filter_list[["cut_ions"]] <- cut_ions
  ion_filter_list[["merge_groups"]] <- merge_groups
  
  if (isTRUE(merge_peaks)) {
    #data_table_merged <- data_table[which(!(data_table$Compound %in% cut_ions)), ]
    #return(as.data.frame(data_table_merged[order(data_table_merged$Compound, decreasing = FALSE), ]))
    return(merge_ions(data_table, ion_filter_list))
  }
  return(ion_filter_list)

}
# data_table <- data.table(peak_df)
# ion_filter_list <- relfil_ion_list
merge_ions <- function(data_table, ion_filter_list) {  
  dat <- melt(data_table, id.vars = c("Compound", "mz", "rt", "kmd"), variable.name = "sample", value.name = "intensity", variable.factor = FALSE)

  for (ion in names(ion_filter_list$merge_groups)) {
    dat <- dat[Compound %in% c(ion, ion_filter_list$merge_groups[[ion]]), intensity:=sum(intensity), by = .(sample)]
  }
  return( dcast(dtm, Compound + mz + kmd + rt ~ sample, value.var = "intensity")[
    (!Compound %in% ion_filter_list$cut_ions), ])
  }

filter_blank <- function(data_table, metadata) {
  
  # extract feature sample (rows are features, columns are samples)
  # data_table <- peak_df_relfil
  # metadata <- full_meta
  ft <- as.data.frame(data_table[ , !(colnames(data_table) %in% c("mz", "rt", "kmd"))])
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
    pivot_longer(cols = c(as.character(data_table$Compound[1]):as.character(tail(data_table$Compound, 1))), names_to = "Compound") %>%
    group_by(Compound, Biological_Group) %>%
    summarise(average = mean(value),
              biolRSD = if_else(average != 0, sd(value) / mean(value), 0),
              bioln = n()) %>%
    ungroup()
  
  tech_stats <- ft_t %>%
    pivot_longer(cols = c(as.character(data_table$Compound[1]):as.character(tail(data_table$Compound, 1))), names_to = "Compound") %>%
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
filter_blank_2 <- function(data_table, metadata) {

  # extract feature sample (rows are features, columns are samples)
  data_table <- peak_df_relfil
  metadata <- full_meta
  ft <- as.data.frame(data_table[ , !(colnames(data_table) %in% c("mz", "rt", "kmd"))])
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

# data_table <- peak_df_relfil
# metadata <- full_meta
filter_blank_3 <- function(data_table, metadata) {
  # extract feature sample (rows are features, columns are samples)
  ft <- data_table[ , !(colnames(data_table) %in% c("mz", "rt", "kmd"))]
  
  rownames(ft) <- as.character(ft$Compound)
  ft <- ft[ , !(colnames(ft) %in% c("Compound"))]
  
  # transpose so rows are samples and columns are features
  ft_t <- as.data.frame(t(ft)) 
  # ft_t <- tibble::rownames_to_column(ft_t, var = "filename")
  
  # make sure metadata rows in the same order as ft_t rownames
  # metadata <- metadata[order(match(metadata$Injection, rownames(ft_t))), ]
  
  groups_stats_list <- vector(mode="list", length = length(unique(metadata$Biological_Group)))
  names(groups_stats_list) <- unique(metadata$Biological_Group)
  lapply(names(groups_stats_list), 
         function(x) {
           dat <- ft_t[which(rownames(ft_t) %in% 
                                metadata$Injection[which(metadata$Biological_Group == x)]), ]

          tech_lvls <- as.factor(metadata$Sample_Code[which(rownames(dat) %in% metadata$Injection)])
          
          groups_stats_list[[x]]$Compound <<- colnames(dat)
          groups_stats_list[[x]]$Biological_Group <<- rep(x, ncol(dat))
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

# microbenchmark::microbenchmark(filter_blank(peak_df_relfil, full_meta), filter_blank_3(peak_df_relfil, full_meta))
# library(data.table)
# data_table <- data.table(peak_df_relfil)
# metadata <- full_meta

#' @import data.table
filter_blank_4 <- function(data_table, metadata) {
  b <- data.table::melt(data_table, id.vars = c("Compound", "mz", "rt", "kmd"), variable.name = "sample", value.name = "intensity", variable.factor = FALSE)[
    metadata, on = .(sample = Injection)][
      , .(average = mean(intensity), BiolRSD = rsd(intensity), Bioln = length(intensity)), by = .(Compound, Biological_Group)]
  
  t <- melt(data_table, id.vars = c("Compound", "mz", "rt", "kmd"), variable.name = "sample", value.name = "intensity", variable.factor = FALSE)[
    metadata, on = .(sample = Injection)][
      , .(sd = rsd(intensity), n = length(intensity)), by = .(Compound, Biological_Group, Sample_Code)][
        , .(techRSD = mean(sd), techn = mean(n)), by = .(Compound, Biological_Group)]

  group_stats <- b[t, on = .(Compound, Biological_Group)]
  setorder(group_stats, Compound, Biological_Group)
  
  return(group_stats)
}

# filter_blank_4(data.table(peak_df_relfil), full_meta)
# microbenchmark::microbenchmark(filter_blank(data_table, metadata), filter_blank_3(data_table, metadata), filter_blank_4(data_table, metadata), times = 1)

# parse ions by group
parse_ions_by_group <- function(group_stats, group_threshold = 0.01) {
  
  # if we pass group_stats_list instead of group_stats (data_table)
  # we would use lines 284-288 in place of the tidy-version (290-297)
  # avgs <- lapply(groups_stats_list, function(x){ 
  #                                     x[[vars]]}
  #      )
  # avgs_df <- data.table::cbindlist(lapply(avgs, as.data.frame))
  # group_stats_tbl <- t(as.data.frame(do.call(cbind, avgs_df)))

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

apply_group_filter  <- function(data_table, group_filter_list, group, remove_ions = TRUE) {
  if(isFALSE(remove_ions)) {
    return(data_table)
  }

  ions <- group_filter_list[[group]]
  data_table <- data_table[!(data_table$Compound %in% ions), ]
  return(data_table)

}



### CV filter
cv_filter <- function(data_table, metadata, cv_threshold, cv_param) {
    # if mismatch and merge and group filters (with remove_ions == TRUE) 
    # are selected, this is expecting the filtered df as input

  # extract feature sample (rows are features, columns are samples)
  ft <- as.data.frame(data_table[ , !(colnames(data_table) %in% c("mz", "rt", "kmd"))])
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
    pivot_longer(cols = c(as.character(data_table$Compound[1]):as.character(tail(data_table$Compound, 1))), names_to = "Compound") %>%
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

# identify_insource_ions(peak_df_filtered, 0.95)
# data_table <- peak_df_filtered
#insource ion filter
filter_insource_ions <- function(data_table, cluster_threshold) {
  rt_list <- unique(data_table$rt)
  
  insourcelist <- c()
  df_meta <- data_table[ , colnames(data_table) %in% c("Compound", "mz", "rt")] #msdata_ind
  ftrgrps <- list()
  mergegroups <- list()
  for (rt in rt_list) {
    filtered_df <- data_table[data_table$rt == rt , colnames(data_table) %in% c("Compound", "mz", "rt")]

    if (nrow(filtered_df) != 1) {
      filtered_df2 <- data_table[data_table$rt == rt , !(colnames(data_table) %in% c("Compound", "mz", "rt", "kmd"))]
      filtered_df2<- as.data.frame(t(filtered_df2))
      colnames(filtered_df2) <- filtered_df$Compound
      corr <- stats::cor(filtered_df2, method = c("pearson")) # why correlate 
      dist <- dist(corr, method = "euclidian")
      cluster <- hclust(dist, method = "complete")
      groups_idx <- cutree(cluster, h = 1 - cluster_threshold)
      
      decon_groups <- vector("list", max(groups_idx))
      
      for(i in 1:length(decon_groups)){
        decon_groups[[i]] <- names(which(groups_idx == i))
      }

      for(group in 1:length(decon_groups))
      {
        if(length(decon_groups[[i]]) != 1)
        {

          temp_df <- df_meta[df_meta$Compound %in% decon_groups[[i]],]
          temp_df <- temp_df[order(temp_df$mz, decreasing = TRUE), ]
          fragments_to_remove <- temp_df$Compound[2:nrow(temp_df)]
          insourcelist <- c(insourcelist, fragments_to_remove)
          # singleslist <- c(fragments_to_remove)
        }
      }      
    }
  }

  # remove ions in insourcelist from data_table 
  #- ask MB if dropping fragments without merging intensity is ok
  return(data_table[!(data_table$Compound %in% insourcelist), ])
}
