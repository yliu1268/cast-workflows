# Phenotype preprocessing for AoU

This script takes as input a phenotype ID and runs preprocessing to create a file with phenotype values and phenotype-specific covariates ready to be used for GWAS or other applications.

Usage example:

```
./aou_phenotype_preprocessing.py \
   --phenotype ALT \
   --concept-id 37047736 \
   --units blood \
   --range 0,250 \
   [--drugexposure-covariate-concept-ids 21601855:statin] \
   [--sample samplefile.txt]
```

* The concept ID must be obtained from the AoU workbench (TODO: add instructions)
* `--units blood` allows units: "IU/L", "No matching concept", "international unit per liter", "no value", "unit per liter". Might change this based on units for other phenotypes
* You can optionally provide a list of samples to restrict to. Typically this would a list of samples that passed upstream sample-level QC info.
* You can optionally provide a list of conceptid:conceptname to use as drug exposure covariates

The output file will have columns: "person_id", "phenotype", "age", plus an additional column for each drug exposure named by the "conceptname" provided for that drug. e.g. the example above will have columns "person_id", "phenotype", "age", "statin".