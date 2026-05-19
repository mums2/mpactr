format_by_type <- function(peak_table_path,
                           type_of_peak_table,
                           sample_names = NULL) {
  if (!(type_of_peak_table
        %in% c("Progenesis", "MzMine", "Metaboscape", "None"))) {
    cli::cli_abort("{.var type_of_peak_table} must be one of Progenesis,
                    MzMine, Metaboscape, or none. See documentation
                    for more details.")
  }
  if (isFALSE(grepl("https://", peak_table_path)) &&
        !file.exists(peak_table_path)) {
    cli::cli_abort("Your peak_table is not a valid file path,
                    please enter a new one")
  }
  result <- list()
  if (type_of_peak_table == "Progenesis") {
    result <- progenesis_formatter(peak_table_path)
  } else if (type_of_peak_table == "MzMine") {
    result <- mz_mine_formatter(peak_table_path)
  } else if (type_of_peak_table == "Metaboscape") {
    result <- metaboscape_formatter(peak_table_path, sample_names)
  } else if (type_of_peak_table == "None") {
    peak_table <- data.frame()
    if (!(any(c("data.table", "data.frame") %in% class(peak_table_path)))) {
      peak_table <- fread(peak_table_path,
                          sep = ",")
    } else {
      peak_table <- peak_table_path
    }
    result <- list(
      "peak_table" = peak_table,
      "raw_table" = peak_table
    )
  } else {

  } # default condition = NULL

  non_injection_columns <-
    which(!(colnames(result$peak_table) %in% sample_names))
  colnames(result$peak_table)[non_injection_columns] <-
    tolower(colnames(result$peak_table)[non_injection_columns])
  result
}

progenesis_formatter <- function(peak_table) {
  if (!(any(c("data.table", "data.frame") %in% class(peak_table)))) {
    peak_table <- fread(peak_table,
      sep = ",",
      skip = 2,
    )
  }
  raw_peak_table <- peak_table
  with(peak_table, setnames(
    peak_table,
    c("m/z", "Retention time (min)"),
    c("mz", "rt")
  ))

  list(
    "peak_table" = peak_table,
    "raw_table" = raw_peak_table
  )
}


mz_mine_formatter <- function(peak_table) {
  if (!(any(c("data.table", "data.frame") %in% class(peak_table)))) {
    peak_table <- fread(peak_table,
      sep = ",",
      skip = 2,
    )
  }
  raw_peak_table <- peak_table

  list(
    "peak_table" = peak_table,
    "raw_table" = raw_peak_table
  )
}

metaboscape_formatter <- function(peak_table, sample_names) {
  if (!(any(c("data.table", "data.frame") %in% class(peak_table)))) {
    peak_table <- fread(peak_table,
      sep = ","
    )
  }
  peak_table_convert <- data.table::copy(peak_table)
  peak_table_convert <- with(peak_table_convert, peak_table_convert[
    , ion := gsub(
      ".*\\[(.+)\\].*", "\\1",
      ADDUCT
    )
  ][
    , charge_string := gsub(".*\\](.+)", "\\1", ADDUCT)
  ][
    charge_string == "+", charge := 1
  ][
    charge_string == "2+", charge := 2
  ][
    charge_string == "3+", charge := 3
  ][
    utils::read.csv(system.file("extdata/ion_masses",
      "DefinedIons.csv",
      package = "mpactr"
    )),
    on = .(ion = IONS),
    nomatch = NULL
  ][
    , mz := (PEPMASS / charge) + MASS
  ])

  with(peak_table_convert, setnames(
    peak_table_convert,
    c("FEATURE_ID", "RT"),
    c("compound", "rt")
  ))
  peak_table_mpactr <- with(peak_table_convert, peak_table_convert[
    , .SD,
    .SDcols = c(
      "compound", "mz", "rt",
      sample_names
    )
  ])

  list(
    "peak_table" = peak_table_mpactr,
    "raw_table" = peak_table_convert
  )
}
