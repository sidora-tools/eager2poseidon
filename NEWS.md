# eager2poseidon 0.4.0

* Update C14 quickcalibrate backend for compatibility with `poseidonR` v0.9.0.

# eager2poseidon 0.3.1

* x-/y-rates from sexdeterrine are now reported in `Sex_Determination_Note` field, and updated when the Genetic_Sex is updated or the field is missing.
* In cases where more than one set of sex determination results exists in the general stats table, use the first available. Generally this will be at the Sample level, which is a good default. In cases where both ds and ss libraries are used, one of the two is picked up depending on ordering. However, the resuls should not vary much between libraries, so this should not affect results much.

# eager2poseidon 0.3.0

* Added a `NEWS.md` file to track changes to the package.
* `Genetic_Sex` field is now filled based on eager sexdeterrmine results.
* x-/y-rates from sexdeterrine are now added to the `Note` field if that is empty.
