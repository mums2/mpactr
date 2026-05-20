mpactr$set("public", "setup", function() {
  private$initialize_data()
})
mpactr$set("private", "initialize_data", function() {
  private$peak_table <- private$peak_table[which(rowSums(
    private$peak_table[, .SD, .SDcols = private$metadata$injection]
  ) > 0), ]
  private$set_kmd()
  private$peak_table$compound <- as.character(private$peak_table$compound)
})
mpactr$set("private", "set_kmd", function() {
  private$peak_table[, kmd := mz - floor(mz)]
})
mpactr$set("public", "get_peak_table", function() {
  private$peak_table
})
mpactr$set("public", "set_peak_table", function(peak_table) {
  private$peak_table <- peak_table
})
mpactr$set("public", "get_metadata", function() {
  private$metadata
})
mpactr$set("public", "get_raw_data", function() {
  private$raw_peak_table
})
