# Running GWAS in All of Us

**under construction**

**Note: all of this needs to be run on the AoU workbench.**

## Important file locations:

* Samples passing QC: `${WORKSPACE_BUCKET}/samples/passing_samples_v7.csv`
* Preprocessed phenotypes+covariates: `${WORKSPACE_BUCKET}/phenotypes/${phenotype}_phenocovar.csv`
* Manifest of existing preprocessed phenotypes: `phenotypes_manifest.csv` (this repo)
* Manifest of existing GWAS results: `gwas_manifest.csv` (this repo)

## Setup

The GWAS instructions below assume you have:

1. Run sample preprocessing to get a set of samples passing all filters. See the `sample_preprocessing` folder for more details.

We will store the list of samples passing QC at:

```
${WORKSPACE_BUCKET}/samples/passing_samples_v7.csv
```

2. Run phenotype preprocessing on your phenotype of interest. See the `phenotype_preprocessing` folder for more details.

We will store phenotypes+shared covariates at:
```
${WORKSPACE_BUCKET}/phenotypes/${phenotype}_phenocovar.csv
```

## Running a SNP-based GWAS using Hail

TODO - what kind of environment do we have to run this in if using hail?

```
./aou_gwas.py \
  --phenotypes ${phenotype} \
  --method hail
```