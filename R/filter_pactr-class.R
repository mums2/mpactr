filter_pactr <- R6Class("filter_pactr", public = list(
  mpactr_data = NA,
  logger = NA,
  initialize = function(mpactr) {
    self$mpactr_data <- mpactr
    self$logger <- new.env(hash = TRUE)
    self$logger[["list_of_summaries"]] <- list()
  },
  print = function() {
    print(self$mpactr_data$get_peak_table())
  },
  is_filter_run = function(filter, group = NULL) {
    if (!is.null(group)) {
      filter <- paste(filter, group, sep = "-")
    }

    if ((filter %in% names(self$logger$list_of_summaries))) {
      return(TRUE)
    }
    FALSE
  },
  get_log = function(filter, group = NULL) {
    if (!(filter %in% c("mispicked", "group", "replicability", "insource"))) {
      cli::cli_abort("{.var filter} must be one of mpactr's supported filters:
                     mispicked, group, replicability, insource.")
    }

    if (!is.null(group)) {
      filter <- paste(filter, group, sep = "-")
    }

    if (!(filter %in% names(self$logger$list_of_summaries))) {
      cli::cli_abort("{.var filter} {filter} has not yet been applied to
                     the data. Run the corresponding filter function prior
                     to extracting the summary.")
    }

    list(
      "failed_ions" = self$logger$list_of_summaries[[filter]]$get_failed_ions(),
      "passed_ions" = self$logger$list_of_summaries[[filter]]$get_passed_ions()
    )
  },
  get_mispicked_ions = function() {
    if (!exists("check_mismatched_peaks", self$logger)) {
      cli::cli_abort("The mispicked filter has not yet been applied to
                     the data - run filter_mispicked_ions() first.")
    }

    merge_groups <- self$logger$check_mismatched_peaks$merge_groups

    similar_ions <- data.table(
      "main_ion" = names(merge_groups),
      "similar_ions" = merge_groups
    )

    similar_ions
  },
  get_group_averages = function() {
    # return averages and variations for all ions in filtered table
    b <- data.table::melt(self$mpactr_data$get_peak_table(),
      id.vars = c("compound", "mz", "rt", "kmd"), variable.name =
        "sample", value.name = "intensity", variable.factor = FALSE
    )[
      data.table(self$mpactr_data$get_metadata()),
      on = .(sample = injection)
    ][
      , .(
        average = mean(intensity),
        BiolRSD = rsd(intensity),
        Bioln = length(intensity)
      ),
      by = .(compound, biological_group)
    ]

    t <- data.table::melt(self$mpactr_data$get_peak_table(),
      id.vars = c("compound", "mz", "rt", "kmd"),
      variable.name = "sample",
      value.name = "intensity",
      variable.factor = FALSE
    )[
      data.table(self$mpactr_data$get_metadata()),
      on = .(sample = injection)
    ][
      , .(sd = rsd(intensity), n = length(intensity)),
      by = .(compound, biological_group, sample_code)
    ][
      , .(techRSD = mean(sd), techn = mean(n)),
      by = .(compound, biological_group)
    ]

    group_stats <- b[t, on = .(compound, biological_group)]
    setorder(group_stats, compound, biological_group)

    group_stats
  },
  get_cv = function() {
    if (!exists("cv_values", self$logger)) {
      cli::cli_abort("The cv filter has not yet been applied
                      to the data - run filter_cv() first.")
    }

    self$logger$cv_values
  }
))
