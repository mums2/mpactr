mpactr$set("public", "setup", function()
{
  private$initialize_data()
})
mpactr$set("private", "initialize_data", function()
{
  setnames(private$peak_table,
           c("m/z", "Retention time (min)"), c("mz", "rt"))
  private$peak_table <- private$peak_table[which(rowSums(
    private$peak_table[, .SD, .SDcols = private$meta_data$Injection]) > 0),]
  private$set_kmd()
})
mpactr$set("private", "set_kmd", function()
{
  private$peak_table[, kmd := mz - floor(mz)]
})
mpactr$set("public", "get_peak_table", function() # make a R facing accessor and make private?
{
  return(private$peak_table)
})
mpactr$set("public", "set_peak_table", function(peak_table) 
{
  private$peak_table <- peak_table
})
mpactr$set("public", "get_meta_data", function()
{
  return(private$meta_data)
})