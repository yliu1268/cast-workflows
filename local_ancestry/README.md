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
    mv "$file" "bathes/batch$counter"
    ((counter++))
done
```

# Test

```
mkdir batches-small/
head -n 10 aou_sample_list.txt > batches-small/batch1
head -n 20 aou_sample_list.txt | tail -n 10 > batches-small/batch2
chrom=21; ./gnomix_launcher.py \
  --name test \
  --vcf ${WGS_ACAF_THRESHOLD_VCF_PATH}/acaf_threshold.chr${chrom}.vcf.bgz \
  --chrom ${chrom} \
  --sample-batches batches-small/ --max-batches 2 --extra-subset-args "-r chr21:1-5252900"

# another test
head -n 11 AFR_BLACK.csv | tail -n 10 | cut -f 1 -d',' > batches-ancestry-small/batchAFR.txt
head -n 11 EUR_WHITE.csv | tail -n 10 | cut -f 1 -d',' > batches-ancestry-small/batchEUR.txt
chrom=1; ./gnomix_launcher.py \
  --name test-ancestry \
  --vcf ${WGS_ACAF_THRESHOLD_VCF_PATH}/acaf_threshold.chr${chrom}.vcf.bgz \
  --chrom ${chrom} \
  --sample-batches batches-ancestry-small/ --max-batches 2 --extra-subset-args "-r chr1:1-44157423"

```

