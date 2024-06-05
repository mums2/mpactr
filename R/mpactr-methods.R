mpactr$set("public", "setup", function()
{
  private$initialize_data()
})
mpactr$set("private", "initialize_data", function()
{
  setnames(self$peak_table,
           c("m/z", "Retention time (min)"), c("mz", "rt"))
  self$peak_table <- self$peak_table[which(rowSums(
    self$peak_table[, .SD, .SDcols = self$meta_data$Injection]) > 0),]
  private$set_kmd()
})
mpactr$set("private", "set_kmd", function()
{
  self$peak_table[, kmd := mz - floor(mz)]
})
