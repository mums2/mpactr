# TODO: Generate a class-like structure for the tree-plot
# TODO: Add filter step to CV filter
# TODO: Outline data flow through CV filter - gdj
# TODO: Outline data flow through blank filter - arm
# TODO: Compare CV filter (with, without solvent blanks, solvent and media blanks as 1 group - need data)
# TODO: Outline data flow through filters (both as we see in mpact and how we envision mpactR)


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
  return( dcast(dat, Compound + mz + kmd + rt ~ sample, value.var = "intensity")[
    (!Compound %in% ion_filter_list$cut_ions), ])
  }



#' @import data.table
filter_blank <- function(data_table, metadata) {
  b <- data.table::melt(data_table, id.vars = c("Compound", "mz", "rt", "kmd"), variable.name = "sample", value.name = "intensity", variable.factor = FALSE)[
    data.table(metadata), on = .(sample = Injection)][
      , .(average = mean(intensity), BiolRSD = rsd(intensity), Bioln = length(intensity)), by = .(Compound, Biological_Group)]

  t <- data.table::melt(data_table, id.vars = c("Compound", "mz", "rt", "kmd"), variable.name = "sample", value.name = "intensity", variable.factor = FALSE)[
    data.table(metadata), on = .(sample = Injection)][
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

  group_avgs <- dcast(group_stats, Compound ~ Biological_Group, value.var = "average")

  max <- group_stats[, .(Compound, Biological_Group, average)][
     , .(max = max(average)), by = Compound]

  group_max <- group_avgs[max, on = .(Compound = Compound)][
    , lapply(.SD, function(x) {x / max}), .SDcols = unique(group_stats$Biological_Group)]

  biol_groups <- as.list(group_max)
  biol_groups <- lapply(biol_groups, setNames, group_avgs[ , Compound])
  group_filter_list <- lapply(biol_groups, \(x){names(x)[which(x > group_threshold)]})
  return(group_filter_list)
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
# data_table <- peak_df_filtered
# metadata <- full_meta
cv_filter <- function(data_table, metadata, cv_threshold, cv_param) {
  cv <- data.table::melt(data_table, id.vars = c("Compound", "mz", "rt", "kmd"), variable.name = "sample", value.name = "intensity", variable.factor = FALSE)[
    data.table(metadata), on = .(sample = Injection)][
      , .(cv = rsd(intensity)), by = .(Compound, Biological_Group, Sample_Code)][
        , .(mean_cv = mean(cv, na.rm = TRUE), median_cv = median(cv, na.rm = TRUE)), by = .(Compound)]

    if (cv_param == "mean") {
      fail_cv <- cv[mean_cv > cv_threshold, Compound]
      return(fail_cv)
    }
    if (cv_param == "median") {
       fail_cv <- cv[median_cv > cv_threshold, Compound]
      return(fail_cv)
    }
}

filter_insouce_ions <- function(data_table, cluster_threshold) {
  data_table_filterd <- data_table[ , cor:= deconvolute_correlation(.SD, 0.95), by = .(rt)][
  cor ==TRUE, ]

  return(data_table_filterd)
}
# filter_insouce_ions_pat(data_table, 0.95)
deconvolute_correlation <- function(group_1, cluster_threshold = 0.95) {
  # return(rep(TRUE, nrow(group_1)))
  if (nrow(group_1) <= 1) {
    return(TRUE)
  }

  dat <- melt(group_1, id.vars = c("Compound", "mz", "kmd"),
              variable.name = "sample", value.name = "intensity", variable.factor = FALSE)[
                                                                                           , c("Compound", "sample", "intensity")]
  data <- dcast(dat, sample ~ Compound, value.var = "intensity")
  data[ , c("sample"):=NULL]
  corr <- stats::cor(data, method = c("pearson"))
  dist <- dist(corr, method = "euclidian")
  cluster <- hclust(dist, method = "complete")
  cut_tree <- cutree(cluster, h = 1 - cluster_threshold)

  x <- as.data.table(cut_tree, keep.rownames = "Compound")[
    , Compound:=as.numeric(Compound)][
      group_1, on = .(Compound = Compound)][
        , keep:=cluster_max(mz), by = .(cut_tree)
      ]

  return(x$keep)
}

cluster_max <- function(mz) {
    keep = rep(FALSE, length(mz))

    keep[which.max(mz)] <- TRUE

    return(keep)
  }
