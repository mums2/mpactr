# Get groups averages.

`get_group_averages()` is a wrapper function to return group averages
for the filtered peak table.

## Usage

``` r
get_group_averages(mpactr_object)
```

## Arguments

- mpactr_object:

  The mpactr object that is created by calling the import_data()
  function.

## Value

a `data.table` reporting the average and relative standard deviation
across biological groups and technical replicates within each group.

## Examples

``` r
data <- import_data(example("coculture_peak_table.csv"),
  example("metadata.csv"),
  format = "Progenesis"
)

data_filter <- filter_group(data, group_to_remove = "Blanks")
#> ℹ Parsing 1303 peaks based on the sample group: Blanks.
#> ℹ Argument remove_ions is: TRUE.Removing peaks from Blanks.
#> ✔ 796 ions failed the Blanks filter, 507 ions remain.

group_averages <- get_group_averages(data_filter)
head(group_averages)
#>    Compound  Biological_Group    average    BiolRSD Bioln    techRSD techn
#>      <char>            <char>      <num>      <num> <int>      <num> <num>
#> 1:     1000 ANG18 monoculture      0.000         NA     3         NA     3
#> 2:     1000 ANGDT monoculture   3762.133 0.01990561     3 0.01990561     3
#> 3:     1000            Blanks      0.000         NA     3         NA     3
#> 4:     1000         Coculture      0.000         NA     3         NA     3
#> 5:     1000   JC1 monoculture      0.000         NA     3         NA     3
#> 6:     1000  JC28 monoculture 182681.783 0.01643330     3 0.01643330     3
```
