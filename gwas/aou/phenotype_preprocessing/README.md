# Phenotype preprocessing for AoU

This script takes as input a phenotype ID and runs preprocessing to create a file with phenotype values and phenotype-specific covariates ready to be used for GWAS or other applications.

Usage example:

```
./aou_phenotype_preprocessing.py \
   --phenotype ALT \
   --concept-id 37047736 \
   --units blood \
   --range 0,250 \
   [--sample samplefile.txt]
```

* The concept ID must be obtained from the AoU workbench (TODO: add instructions)
* `--units blood` allows units: "IU/L", "No matching concept", "international unit per liter", "no value", "unit per liter"
* You can optionally provide a list of samples to restrict to. Typically this would a list of samples that passed upstream sample-level QC info.

Notes:
* get trait-specific covars (LDL)
* Take in list of blessed samples from Tara's script
* Run on all phenotypes and store in shared bucket