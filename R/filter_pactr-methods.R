filter_pactr$set("public", "check_mismatched_peaks", function(ringwin, isowin, trwin, max_iso_shift, merge_peaks) {
  number_of_input_ions <- nrow(self$mpactr_data$peak_table)

  ion_filter_list <- list()
  cut_ions <- c() # list
  merge_groups <- list() # dictonary

  self$mpactr_data$peak_table <- self$mpactr_data$peak_table[order(self$mpactr_data$peak_table$mz, decreasing = FALSE),]
  number_of_rows <- nrow(self$mpactr_data$peak_table)
  for (i in seq_along(1:number_of_rows)) {
    current_feature <- self$mpactr_data$peak_table[i, c("Compound", "mz", "rt")]
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
        if (self$mpactr_data$peak_table$Compound[j] %in% cut_ions) {
          next
        }
        mass_diff <- self$mpactr_data$peak_table$mz[j] - current_mass
        kmd_diff <- mass_diff - floor(mass_diff)

        if (abs(mass_diff) > max_iso_shift - 0.4) { # BL  - why 0.4??
          break
        }

        rt_diff <- self$mpactr_data$peak_table$rt[j] - current_rt
        ring_band <- floor(abs(mass_diff) * (1 / ringwin)) %% (1 / ringwin)
        double_band <- kmd_diff -
          .5004 -
          (floor(mass_diff) * .007) # BL - why 0.500 and 0.007?

        if ((abs(rt_diff) <= trwin) &&
          (mass_diff <= max_iso_shift - 0.4) &&
          (ring_band == 0 ||
            double_band < isowin)) {
          cut_ions <- c(cut_ions, self$mpactr_data$peak_table$Compound[j])
          merge_groups[[as.character(current_ion)]] <- c(merge_groups[[as.character(current_ion)]], self$mpactr_data$peak_table$Compound[j])
        }
      }
    }
  }
  # TODO Look into removing merge groups
  ion_filter_list[["cut_ions"]] <- cut_ions
  ion_filter_list[["merge_groups"]] <- merge_groups
  self$logger[["check_mismatched_peaks"]] <- ion_filter_list
  if (isTRUE(merge_peaks)) {
    private$merge_ions(ion_filter_list)
  }

  number_of_failing_ions <- number_of_input_ions - nrow(self$mpactr_data$peak_table)
  remaining_ions <- nrow(self$mpactr_data$peak_table)
  merged_peak_string <- paste0("merge peaks was set to: ", merge_peaks, ", ", number_of_failing_ions,
                               " failing peaks were not removed. ", remaining_ions, " remain")
  if (isTRUE(merge_peaks)) {
    merged_peak_string <- paste0("merge peaks was set to: ", merge_peaks, ", ", number_of_failing_ions,
                                 " failing peaks were removed. ", remaining_ions, " remain")
  }
  self$logger[["check_mismtach_peaks_summary"]] <- merged_peak_string
  print(merged_peak_string)
})
filter_pactr$set("public", "filter_blank", function() {
  b <- data.table::melt(self$mpactr_data$peak_table, id.vars = c("Compound", "mz", "rt", "kmd"), variable.name =
    "sample", value.name = "intensity", variable.factor = FALSE)[
    data.table(self$mpactr_data$meta_data), on = .(sample = Injection)][
    , .(average = mean(intensity), BiolRSD = rsd(intensity), Bioln = length(intensity)), by = .(Compound, Biological_Group)]

  t <- data.table::melt(self$mpactr_data$peak_table, id.vars = c("Compound", "mz", "rt", "kmd"), variable.name = "sample", value.name = "intensity", variable.factor = FALSE)[
    data.table(self$mpactr_data$meta_data), on = .(sample = Injection)][
    , .(sd = rsd(intensity), n = length(intensity)), by = .(Compound, Biological_Group, Sample_Code)][
    , .(techRSD = mean(sd), techn = mean(n)), by = .(Compound, Biological_Group)]

  group_stats <- b[t, on = .(Compound, Biological_Group)]
  setorder(group_stats, Compound, Biological_Group)
  self$logger[["group_filter-group_stats"]] <- group_stats
})
filter_pactr$set("public", "parse_ions_by_group", function(group_threshold = 0.01) {
  group_avgs <- dcast(self$logger[["group_filter-group_stats"]], Compound ~ Biological_Group, value.var = "average")

  max <- self$logger[["group_filter-group_stats"]][, .(Compound, Biological_Group, average)][
    , .(max = max(average)), by = Compound]

  group_max <- group_avgs[max, on = .(Compound = Compound)][
    , lapply(.SD, function(x) { x / max }),
      .SDcols = unique(self$logger[["group_filter-group_stats"]]$Biological_Group)]

  biol_groups <- as.list(group_max)
  biol_groups <- lapply(biol_groups, setNames, group_avgs[, Compound])
  group_filter_list <- lapply(biol_groups, \(x) { names(x)[which(x > group_threshold)] })
  self$logger[["group_filter-failing_list"]] <- group_filter_list
})
filter_pactr$set("public", "apply_group_filter", function(group, remove_ions = TRUE) {
  number_of_input_ions <- nrow(self$mpactr_data$peak_table)
  if (isFALSE(remove_ions)) {
    print(paste0("remove_ions  was set to: ", remove_ions, " therefore, no data has been changed"))
    return()
  }
  ions <- self$logger[["group_filter-failing_list"]][[group]]
  self$mpactr_data$peak_table <- self$mpactr_data$peak_table[!(self$mpactr_data$peak_table$Compound %in% ions),]

  number_of_failing_ions <- length(ions)
  remaining_ions <- nrow(self$mpactr_data$peak_table)
  group_peak_string <- paste0("remove_ions was set to: ", remove_ions, ", ", number_of_failing_ions,
                              " failing peaks were removed. ", remaining_ions, " remain")
  self$logger[["group_filter_summary"]] <- group_peak_string
  print(group_peak_string)
})
filter_pactr$set("private", "merge_ions", function(ion_filter_list) {
  dat <- melt(self$mpactr_data$peak_table, id.vars = c("Compound", "mz", "rt", "kmd"), variable.name = "sample", value.name = "intensity", variable.factor = FALSE)

  for (ion in names(ion_filter_list$merge_groups)) {
    dat <- dat[Compound %in% c(ion, ion_filter_list$merge_groups[[ion]]), intensity := sum(intensity), by = .(sample)]
  }
  self$mpactr_data$peak_table <- dcast(dat, Compound + mz + kmd + rt ~ sample, value.var = "intensity")[
    (!Compound %in% ion_filter_list$cut_ions),]
})
