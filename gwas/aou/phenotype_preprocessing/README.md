# Phenotype preprocessing for AoU

This script takes as input a phenotype ID and runs preprocessing to create a file with phenotype values and phenotype-specific covariates ready to be used for GWAS or other applications.

Usage example:

```
./aou_phenotype_preprocessing.py --phenotype ALT [--sample samplefile.txt]
```

Usage notes:

* Before you run a new phenotype for the first time, you must edit the variable `pt_info` in `aou_queries.py` to add the following info:
  * phenotype ID (e.g. "ALT")
  * concept ID. You'll have to get this from the AoU cohort browser (TODO: add instructions for constructing a concept ID)
  * Acceptable units for that trait. `BLOOD_UNITS` contains common units as a starting point.
  * The acceptable range of values (min, max) for the trait.
* You can optionally provide a list of samples to restrict to. Typically this would a list of samples that passed upstream sample-level QC info.


Notes:
* get trait-specific covars (LDL)
* Take in list of blessed samples from Tara's script
* Run on all phenotypes and store in shared bucket