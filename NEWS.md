# eager2poseidon 0.3.6

## [0.3.6] - 2023-05-10

### `Added`

### `Fixed`

- `fill_in_janno()`: Now correctly keeps `Relation_*` columns.
- `fill_genetic_sex()`: Should now work with newer R and pacakge versions.
- `read_eager_stats_table()`: Eager tables are now joined with eager TSV to allow for sample names that include `_` before the first `.`.
- `read_eager_stats_table()`: Use `sprintf` instead of `format` for formatting decimals. Should avoid `digits=0` error in newer R versions.
- `eager2poseidon.R`: Added explicit `.cols` statement to `mutate`. needed for newer dplyr versions.

### `Dependencies`

### `Deprecated`

# eager2poseidon 0.3.5

## [0.3.5] - 2023-03-01

### `Added`

### `Fixed`

- `compile_eager_result_tables()`: Corrected column name for number of libraries. `Nr_Libs` -> `Nr_Libraries`. Closes #4 
- `compile_across_lib_results()`: Inference of `Nr_Libraries` corrected to count only unique values. This also fixed inference of `Capture_Type`. Closes #5 
- `compile_across_lib_results()`: Corrected weighted sum calculation to only include unique rows of library-level data. Closes #6 
- `compile_across_lib_results()`: Corrected `Library_Names` field. Now includes only unique library names.

### `Dependencies`

- Added `vctrs` to package dependencies

### `Deprecated`

# eager2poseidon 0.3.2

* Update C14 quickcalibrate backend for compatibility with `poseidonR` v0.9.0.

# eager2poseidon 0.3.1

* x-/y-rates from sexdeterrine are now reported in `Sex_Determination_Note` field, and updated when the Genetic_Sex is updated or the field is missing.
* In cases where more than one set of sex determination results exists in the general stats table, use the first available. Generally this will be at the Sample level, which is a good default. In cases where both ds and ss libraries are used, one of the two is picked up depending on ordering. However, the resuls should not vary much between libraries, so this should not affect results much.

# eager2poseidon 0.3.0

* Added a `NEWS.md` file to track changes to the package.
* `Genetic_Sex` field is now filled based on eager sexdeterrmine results.
* x-/y-rates from sexdeterrine are now added to the `Note` field if that is empty.
