# Setup

This must be run on the All of Us workbench. You must have launched both an analysis and cromshell environment.

```
git clone https://github.com/cast-genomics/cast-workflows
cd cast-workflows/
git checkout gnomix
cd local_ancestry
../utils/configure-cromshell.py
```

# Set up sample batches of 1000 samples each

```
gsutil cp $WORKSPACE_BUCKET/tr_imputation/tr_imputation/sample/aou_sample_list.txt .
mkdir batches
split -l 1000 aou_sample_list.txt batches/batch
counter=1
for file in batches/batch*
do
    mv "$file" "batches/batch$counter"
    ((counter++))
done
```

# Test

```
chrom=21; ./gnomix_launcher.py \
  --name test \
  --vcf ${WGS_ACAF_THRESHOLD_VCF_PATH}/acaf_threshold.chr${chrom}.vcf.bgz \
  --chrom ${chrom} \
  --sample-batches batches/ --max-batches 2
```