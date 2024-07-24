####  filter 1: mismatched peaks    ###
filter_pactr$set("public", "check_mismatched_peaks", function(ringwin, isowin, trwin, max_iso_shift, merge_peaks, merge_method = NULL) {
  cli::cli_alert_info("Checking {nrow(self$mpactr_data$get_peak_table())} peaks for mispicked peaks.")

  ion_filter_list <- list()
  cut_ions <- c() # list
  merge_groups <- list() # dictonary

  self$mpactr_data$set_peak_table(self$mpactr_data$get_peak_table()[order(self$mpactr_data$get_peak_table()$mz,
                                                                         decreasing =
    FALSE),])
  number_of_rows <- nrow(self$mpactr_data$get_peak_table())
  for (i in seq_along(1:number_of_rows)) {
    current_feature <- self$mpactr_data$get_peak_table()[i, c("Compound", "mz", "rt")]
    current_rt <- current_feature$rt
    current_mass <- current_feature$mz
    current_ion <- current_feature$Compound
    if (!(current_ion %in% cut_ions))
    {
      for (j in (i + 1):number_of_rows) {
        if (j > number_of_rows)
        {
          break
        }
        if (self$mpactr_data$get_peak_table()$Compound[j] %in% cut_ions) {
          next
        }
        mass_diff <- self$mpactr_data$get_peak_table()$mz[j] - current_mass
        kmd_diff <- mass_diff - floor(mass_diff)

        if (abs(mass_diff) > max_iso_shift - 0.4) { # BL  - why 0.4??
          break
        }

        rt_diff <- self$mpactr_data$get_peak_table()$rt[j] - current_rt
        ring_band <- floor(abs(mass_diff) * (1 / ringwin)) %% (1 / ringwin)
        double_band <- kmd_diff -
          .5004 -
          (floor(mass_diff) * .007) # BL - why 0.500 and 0.007?

        if ((abs(rt_diff) <= trwin) &&
          (mass_diff <= max_iso_shift - 0.4) &&
          (ring_band == 0 ||
            double_band < isowin)) {
          cut_ions <- c(cut_ions, self$mpactr_data$get_peak_table()$Compound[j])
          merge_groups[[as.character(current_ion)]] <- c(merge_groups[[as.character(current_ion)]], self$mpactr_data$get_peak_table()$Compound[j])
        }
      }
    }
  }
  # TODO Look into removing merge groups
  ion_filter_list[["cut_ions"]] <- cut_ions
  ion_filter_list[["merge_groups"]] <- merge_groups
  self$logger[["check_mismatched_peaks"]] <- ion_filter_list

  if (isTRUE(merge_peaks)) {
    cli::cli_alert_info("Argument merge_peaks is: {merge_peaks}. Merging mispicked peaks with method {merge_method}.")
    
    private$merge_ions(ion_filter_list, merge_method)
  }
  else{
    cli::cli_alert_warning("Argument merge_peaks is: {merge_peaks}. Mispicked peaks will not be merged.")
  }
  
  self$logger$list_of_summaries$mispicked <- summary$new(filter = "mispicked", failed_ions = cut_ions,
                              passed_ions = self$mpactr_data$get_peak_table()$Compound)
  
  self$logger$list_of_summaries$mispicked$summarize()
})

filter_pactr$set("private", "merge_ions", function(ion_filter_list, method) {
  
  if (is.null(method)) {
    cli::cli_abort("No method has been supplied for merging peaks. method must be one of: sum")
  }
  
  if (method == "sum") {
      dat <- melt(self$mpactr_data$get_peak_table(), id.vars = c("Compound", "mz", "rt", "kmd"), variable.name = "sample", value.name = "intensity", variable.factor = FALSE)

    for (ion in names(ion_filter_list$merge_groups)) {
      dat <- dat[Compound %in% c(ion, ion_filter_list$merge_groups[[ion]]), intensity := sum(intensity), by = .(sample)]
    }
    self$mpactr_data$set_peak_table(dcast(dat, Compound + mz + kmd + rt ~ sample, value.var = "intensity")[
      (!Compound %in% ion_filter_list$cut_ions),])
  }

})

####  filter 2: group filter    ###
# Calculates statisics for each feature (rsd, n) accross biological groups and technical replicates
filter_pactr$set("public", "filter_blank", function() {
  b <- data.table::melt(self$mpactr_data$get_peak_table(), id.vars = c("Compound", "mz", "rt", "kmd"), variable.name =
    "sample", value.name = "intensity", variable.factor = FALSE)[
    data.table(self$mpactr_data$get_meta_data()), on = .(sample = Injection)][
    , .(average = mean(intensity), BiolRSD = rsd(intensity), Bioln = length(intensity)), by = .(Compound, Biological_Group)]

  t <- data.table::melt(self$mpactr_data$get_peak_table(), id.vars = c("Compound", "mz", "rt", "kmd"), variable.name = "sample", value.name = "intensity", variable.factor = FALSE)[
    data.table(self$mpactr_data$get_meta_data()), on = .(sample = Injection)][
    , .(sd = rsd(intensity), n = length(intensity)), by = .(Compound, Biological_Group, Sample_Code)][
    , .(techRSD = mean(sd), techn = mean(n)), by = .(Compound, Biological_Group)]

  group_stats <- b[t, on = .(Compound, Biological_Group)]
  setorder(group_stats, Compound, Biological_Group)
  self$logger[["group_filter-group_stats"]] <- group_stats
})

# this function determines group ions above the group threshold given group statistics (see filter_blank).
# The result is a list of ions by group whose relative abundance is greater than the threshold.
filter_pactr$set("public", "parse_ions_by_group", function(group_threshold = 0.01) {
  group_avgs <- dcast(self$logger[["group_filter-group_stats"]], Compound ~ Biological_Group, value.var = "average")

  max <- self$logger[["group_filter-group_stats"]][, .(Compound, Biological_Group, average)][
    , .(max = max(average)), by = Compound]

  group_max <- group_avgs[max, on = .(Compound = Compound)][
    , lapply(.SD, function(x) { x / max }),
      .SDcols = unique(self$logger[["group_filter-group_stats"]]$Biological_Group)]

  biol_groups <- as.list(group_max)
  biol_groups <- lapply(biol_groups, setNames, group_avgs[, Compound])
  group_filter_list <- lapply(biol_groups, function(x) { names(x)[which(x > group_threshold)] })
  self$logger[["group_filter-failing_list"]] <- group_filter_list
})

# Given a group name, removes flagged ions from the peak table.
filter_pactr$set("public", "apply_group_filter", function(group, remove_ions = TRUE) {
  cli::cli_alert_info("Parsing {nrow(self$mpactr_data$get_peak_table())} peaks based on the following sample group: {group}.")

  if (isFALSE(remove_ions)) {
    cli::cli_alert_warning("Argument remove_ions is {remove_ions}. Peaks from {group} will not be removed.")
    return()
  }
  
  cli::cli_alert_info("Argument remove_ions is: {remove_ions}. Removing peaks from {group}.")
  
  ions <- self$logger[["group_filter-failing_list"]][[group]]
  self$mpactr_data$set_peak_table(self$mpactr_data$get_peak_table()[!(self$mpactr_data$get_peak_table()$Compound
    %in% ions),])

   self$logger$list_of_summaries[[paste0("group-",group)]] <- summary$new(filter = group,
                                                              failed_ions = as.numeric(ions),
                                                              passed_ions = self$mpactr_data$get_peak_table()$Compound)


  self$logger$list_of_summaries[[paste0("group-",group)]]$summarize()
})

####  filter 3: cv filter    ###
filter_pactr$set("public", "cv_filter", function(cv_threshold = NULL, cv_params) {
  if (is.null(cv_threshold)) {
    cli::cli_abort("{.var cv_threshold} must be supplied.")
  }
  ## abort if an incorrect cv_params argument is supplied.
  if (!(cv_params %in% c("mean", "median"))) {
    cli::cli_abort("{.var cv_params} must be one of mean or median.")
  }
  
  ## abort if there are no technical replicates.
  if (isFALSE(self$mpactr_data$isMultipleTechReps())) {
      cli::cli_abort("There are no technical replicates in the dataset provided. In order to run the replicability filter, technical replicates are required.")
  }
  
  input_ions <- self$mpactr_data$get_peak_table()$Compound
  cli::cli_alert_info("Parsing {length(input_ions)} peaks for replicability across technical replicates.")

  cv <- data.table::melt(self$mpactr_data$get_peak_table(), id.vars = c("Compound", "mz", "rt", "kmd"), variable.name =
    "sample", value.name = "intensity", variable.factor = FALSE)[
    self$mpactr_data$get_meta_data(), on = .(sample = Injection)][
      , .(cv = rsd(intensity)), by = .(Compound, Biological_Group, Sample_Code)][
        , .(mean_cv = mean(cv, na.rm = TRUE), median_cv = median(cv, na.rm = TRUE)), by = .(Compound)]
  
    self$logger[["cv_values"]] <- cv

    if (cv_params == "mean") {
      failed_ions <- cv[mean_cv > cv_threshold, Compound]
    }
    else {
       failed_ions<- cv[median_cv > cv_threshold, Compound]
    }

    self$mpactr_data$set_peak_table(self$mpactr_data$get_peak_table()[Compound %in% setdiff(input_ions,
                                                                      failed_ions), ])

    self$logger$list_of_summaries$replicability <- summary$new(filter = "cv_filter",
                                                     failed_ions = failed_ions,
                                                     passed_ions = self$mpactr_data$get_peak_table()$Compound)
    
    self$logger$list_of_summaries$replicability$summarize()
})

####  filter 4: insource ions   ###

filter_pactr$set("public", "filter_insource_ions", function(cluster_threshold = 0.95) {

  input_ions <- self$mpactr_data$get_peak_table()$Compound
  cli::cli_alert_info("Parsing {length(input_ions)} peaks for insource ions.")

  self$mpactr_data$set_peak_table(self$mpactr_data$get_peak_table()[ ,
    cor:= private$deconvolute_correlation(.SD, cluster_threshold), by = .(rt)][cor ==TRUE, ])

  self$logger$list_of_summaries$insource <- summary$new(filter = "insource",
                                                          failed_ions = setdiff(input_ions,
                                                                            self$mpactr_data$get_peak_table()$Compound),
                                                          passed_ions = self$mpactr_data$get_peak_table()$Compound)

  self$logger$list_of_summaries$insource$summarize()
})

filter_pactr$set("private", "deconvolute_correlation", function(group_1, cluster_threshold = 0.95) {
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
  dist <- stats::dist(corr, method = "euclidian")
  cluster <- stats::hclust(dist, method = "complete")
  cut_tree <- stats::cutree(cluster, h = 1 - cluster_threshold)

  x <- as.data.table(cut_tree, keep.rownames = "Compound")[
    , Compound:=as.numeric(Compound)][
      group_1, on = .(Compound = Compound)][
        , keep:=private$cluster_max(mz), by = .(cut_tree)
      ]

  return(x$keep)
})

# Given a list of mz values, this function will determine which value has the largest mz and modify its "KEEP" Column.
filter_pactr$set("private", "cluster_max", function(mz) {
  keep <- rep(FALSE, length(mz))
  keep[which.max(mz)] <- TRUE
  return(keep)
})