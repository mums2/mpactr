# Get similar ion groups.

`get_similar_ions()` is a wrapper function to return similar ion groups
determined with the
[`filter_mispicked_ions()`](https://www.mums2.org/mpactr/reference/filter_mispicked_ions.md).

## Usage

``` r
get_similar_ions(mpactr_object)
```

## Arguments

- mpactr_object:

  The mpactr object that is created by calling the import_data()
  function.

## Value

a `data.table` reporting the main ion and those found to be similar with
[`filter_mispicked_ions()`](https://www.mums2.org/mpactr/reference/filter_mispicked_ions.md).

## Examples

``` r
data <- import_data(
  example_path("coculture_peak_table.csv"),
  example_path("metadata.csv"),
  format = "Progenesis"
)

data_filter <- filter_mispicked_ions(data)
#> ℹ Checking 1303 peaks for mispicked peaks.
#> ℹ Argument merge_peaks is: TRUE. Merging mispicked peaks with method sum.
#> ✔ 70 ions failed the mispicked filter, 1233 ions remain.

mispicked_ion_groups <- get_similar_ions(data_filter)
mispicked_ion_groups
#>     main_ion similar_ions
#>       <char>       <list>
#>  1:     1188         1189
#>  2:      939          945
#>  3:      898          896
#>  4:     1214         1210
#>  5:     1253         1249
#>  6:      886          884
#>  7:     1059         1060
#>  8:     1011         1008
#>  9:     1305         1304
#> 10:     1272         1271
#> 11:     1332         1333
#> 12:      944          937
#> 13:     1195         1196
#> 14:     1022         1014
#> 15:      897          895
#> 16:     1233         1228
#> 17:      811          817
#> 18:      851          853
#> 19:     1123         1122
#> 20:      951          938
#> 21:     1146         1145
#> 22:     1019         1021
#> 23:      941          946
#> 24:      878          161
#> 25:     1030    1032,1033
#> 26:     1009          282
#> 27:      153          154
#> 28:      785          779
#> 29:     1256         1254
#> 30:     1216         1224
#> 31:      259          260
#> 32:      263          264
#> 33:      770           18
#> 34:     1004         1012
#> 35:      784          800
#> 36:     1258         1251
#> 37:      903          901
#> 38:      248          969
#> 39:      586         1259
#> 40:     1222         1223
#> 41:      525          528
#> 42:     1103         1095
#> 43:      918          927
#> 44:     1243         1242
#> 45:      411          410
#> 46:     1264         1260
#> 47:      426     427,1157
#> 48:      529         1206
#> 49:      395          393
#> 50:     1247         1252
#> 51:     1203         1204
#> 52:      973          975
#> 53:      583          584
#> 54:     1281         1282
#> 55:     1091         1090
#> 56:     1262         1257
#> 57:      508          507
#> 58:     1089         1087
#> 59:     1155         1153
#> 60:     1263         1261
#> 61:      372      371,373
#> 62:      342          340
#> 63:     1266         1267
#> 64:      545          546
#> 65:      496          495
#> 66:     1288    1289,1290
#>     main_ion similar_ions
```
