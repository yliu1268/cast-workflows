# Running GWAS in All of Us

**under construction**

## Setup

The GWAS instructions below assume you have:

1. Run sample preprocessing to get a set of samples passing all filters (TODO)
2. Run phenotype preprocessing on your phenotype of interest. See the `phenotype_preprocessing` folder for more details.

We will store phenotypes+shared covariates at:
```
${WORKSPACE_BUCKET}/phenotypes/${phenotype}_phenocovar.csv
```

## Running a SNP-based GWAS using Hail

./aou_gwas.py \
  --phenotypes ${phenotype} \
  --method hail