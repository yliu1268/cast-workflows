# Running GWAS in All of Us

**under construction**

**Note: all of this needs to be run on the AoU workbench.**

## Quickstart

Assuming the phenotype you are analyzing has already been preprocessed (see below), you can run GWAS (by default using Hail):

1. Example test on a small region (can be run on default Hail environment)
```
./aou_gwas.py --phenotype ALT --num-pcs 10 --region chr11:119206339-119308149
```
will output `ALT_hail_chr11_119206339_119308149.gwas.tab` 

2. Example full GWAS. (TODO: which Hail compute resources)
```
./aou_gwas.py --phenotype ALT --num-pcs 10
```
will output `ALT_hail.gwas.tab`

3. Small example associaTR run (requires installing TRTools)

```
TRVCF=${WORKSPACE_BUCKET}/cromwell-execution/targetTR/031367d7-36d5-4da7-950f-0b42ae50c91e/call-sort_index/test_10.filtered.sorted.vcf.gz
./aou_gwas.py --phenotype ALT --num-pcs 10 --method associaTR --tr-vcf $TRVCF
```

## Detailed usage

Required arguments:

* `--phenotype <STR>`: phenotype ID to use to name the output file

Optional arguments:

* `--method <STR>`: Which GWAS method to use. Default: hail
* `--samples <FILE>`: csv file with list of samples to restrict to. Needs columns "person_id" and "sex_at_birth_Male". Typically this would a list of samples that passed upstream sample-level QC info. Defaults to `${WORKSPACE_BUCKET}/samples/passing_samples_v7.csv`.
* `--ancestry-pred-path <FILE>`: Path to file with PC info. Defaults to the AoU v7 ancestry predictions.
* `--region <STR>`: Region to restrict to (chrom:start-endd)
* `--num-pcs <INT>`: Number of population PCs to include as covariates. Default: 10
* `--ptcovars <STR>`: Comma-separated list of phenotype-specific covariates. Default: age
* `--sharedcovars <STR>`: Comma-separated list of shared covariates (besides PCs). Default: sex_at_birth_Male
* `--plot`: Output Manhattan and QQ plots

## Important file locations:

* Samples passing QC: `${WORKSPACE_BUCKET}/samples/passing_samples_v7.csv`
* Preprocessed phenotypes+covariates: `${WORKSPACE_BUCKET}/phenotypes/${phenotype}_phenocovar.csv`
* GWAS results: `${WORKSPACE_BUCKET}/gwas/${phenotype}_hail.gwas.tab`
* Manifest of existing preprocessed phenotypes: `phenotypes_manifest.csv` (this repo)

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

