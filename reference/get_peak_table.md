# Return the peak table data.table from the mpactr object.

`get_peak_table()` a wrapper function to return the peak table object of
the given mpactr object.

## Usage

``` r
get_peak_table(mpactr_object)
```

## Arguments

- mpactr_object:

  The mpactr object that is created by calling the import_data()
  function.

## Value

a `data.table`.

## Examples

``` r
data <- import_data(
  example_path("coculture_peak_table.csv"),
  example_path("metadata.csv"),
  format = "Progenesis"
)

peak_table <- get_peak_table(data)
```
