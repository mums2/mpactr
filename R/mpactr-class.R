mpactr <- R6Class("mpactr",
  public = list(
    # Properties

    # Constructor
    initialize = function(peak_table, metadata) {
      stopifnot(any(class(peak_table$raw_table) == "data.table"))
      stopifnot(any(class(peak_table$peak_table) == "data.table"))
      stopifnot(any(class(metadata) == "data.table"))
      private$raw_peak_table <- data.table(peak_table$raw_table)
      private$peak_table <- data.table(peak_table$peak_table)
      private$metadata <- data.table(metadata)
    },
    isMultipleTechReps = function() {
      any(private$metadata[, .N, by = sample_code][["N"]] > 1)
    }
  ),
  private = list(
    peak_table = NA,
    metadata = NA,
    raw_peak_table = NA
  )
)
