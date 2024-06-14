filter_pactr <- R6Class("filter_pactr", public = list(
  mpactr_data = NA,
  logger = NA, # Might have to add getter/setters
  initialize = function(mpactr) {
    self$mpactr_data <- mpactr
    self$logger <- new.env(hash = TRUE)
    self$logger[["list_of_summaries"]] <- list()
  },
  get_log = function(filter) {
    if(!(filter %in% c("mispicked", "group", "replicability", "insource"))) {
      cli::cli_abort("{.var filter} must be one of mpactR's supported filters: mispicked, group, cv, insource.")
    }
    
    if (!(filter %in% names(self$logger$list_of_summaries))) {
      cli::cli_abort("{.var filter} {filter} has not yet been applied to the data. Run the corresponding filter function prior to extracting the summary.")
    }
    
    return(list("failed_ions" = self$logger$list_of_summaries[[filter]]$get_failed_ions(),
                "passed_ions" = self$logger$list_of_summaries[[filter]]$get_passed_ions()))

  },
  get_mispicked_ions = function(){
    if (!exists("check_mismatched_peaks", self$logger)) {
      cli::cli_abort("The mispicked filter has not yet been applied to the data - run filter_mispicked_ions() first.")
    }
    
    merge_groups <- self$logger$check_mismatched_peaks$merge_groups
    
    similar_ions <- data.table("main_ion" = names(merge_groups),
                               "similar_ions" = merge_groups)
    
    return(similar_ions)
  }
  )
)