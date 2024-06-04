
filter_pactr <- R6Class("filter_pactr", public = list(
  mpactr_data = NA,
  logger = NA,
  initialize = function(mpactr){
    self$mpactr_data <- mpactr
    self$logger <- new.env(hash = TRUE)
  },
  cv_filter = function(x)
  {
    # etc
  },
  check_mismatched_peaks = function(ringwin, isowin, trwin, max_iso_shift, merge_peaks) {
  ion_filter_list <- list()
  current_rt <- 0
  current_mass <- 0
  mass_diff <- 0
  kmd_diff <- 0
  cut_ions <- c() # list
  merge_groups <- list() # dictonary

   self$mpactr_data$peak_table_csv <-self$mpactr_data$peak_table_csv[order(self$mpactr_data$peak_table_csv$mz, decreasing = FALSE), ]
  number_of_rows <- nrow(self$mpactr_data$peak_table_csv)
  for(i in seq_along(1:number_of_rows)) {
    current_feature <- self$mpactr_data$peak_table_csv[i, c("Compound", "mz", "rt")]
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
        if(self$mpactr_data$peak_table_csv$Compound[j] %in% cut_ions) {
          next
        }
        mass_diff <- self$mpactr_data$peak_table_csv$mz[j] - current_mass
        kmd_diff <- mass_diff - floor(mass_diff)

        if (abs(mass_diff) > max_iso_shift - 0.4) { # BL  - why 0.4??
          break
        }

        rt_diff <- self$mpactr_data$peak_table_csv$rt[j] - current_rt
        ring_band <- floor(abs(mass_diff) * (1/ringwin)) %% (1/ringwin)
        double_band <- kmd_diff - .5004 -
        (floor(mass_diff) * .007) # BL - why 0.500 and 0.007?

        if((abs(rt_diff) <= trwin) &&
           (mass_diff <= max_iso_shift - 0.4) &&
           (ring_band == 0 ||
            double_band < isowin)) {
            cut_ions <- c(cut_ions, self$mpactr_data$peak_table_csv$Compound[j])
            merge_groups[[as.character(current_ion)]] <- c(merge_groups[[as.character(current_ion)]],self$mpactr_data$peak_table_csv$Compound[j])
        }
      }
    }
  }
  # TODO Look into removing merge groups
  ion_filter_list[["cut_ions"]] <- cut_ions
  ion_filter_list[["merge_groups"]] <- merge_groups
  self$logger[["check_mismatched_peaks"]] <- ion_filter_list
  if (isTRUE(merge_peaks)) {
    #data_table_merged <- data_table[which(!(data_table$Compound %in% cut_ions)), ]
    #return(as.data.frame(data_table_merged[order(data_table_merged$Compound, decreasing = FALSE), ]))

    self$merge_ions(ion_filter_list)
  }
  return(ion_filter_list)

},
merge_ions = function(ion_filter_list) {
  dat <- melt(self$mpactr_data$peak_table_csv, id.vars = c("Compound", "mz", "rt", "kmd"), variable.name = "sample", value.name = "intensity", variable.factor = FALSE)

  for (ion in names(ion_filter_list$merge_groups)) {
    dat <- dat[Compound %in% c(ion, ion_filter_list$merge_groups[[ion]]), intensity:=sum(intensity), by = .(sample)]
  }
  self$mpactr_data$peak_table_csv <- dcast(dat, Compound + mz + kmd + rt ~ sample, value.var = "intensity")[
    (!Compound %in% ion_filter_list$cut_ions), ]
  }



  
  # merge_mismtached_peaks = function(mpactr_object, ringwin, isowin, trwin, max_iso_shift, merge_peaks)
  # {
  #   check_mismatched_peaks(mpactr_object, ringwin, isowin, trwin, max_iso_shift, merge_peaks)
  #   # We would pass the mpact object over to the class (as a reference) and work with the reference of it. Therefore,
  #   # there will be no copy and it will be passing by a reference.
  # }
))
