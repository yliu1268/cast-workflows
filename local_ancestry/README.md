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
mkdir batches-ancestry-small
gsutil cp ${WORKSPACE_BUCKET}/samples/AFR_BLACK.csv .
gsutil cp ${WORKSPACE_BUCKET}/samples/EUR_WHITE.csv .
head -n 11 AFR_BLACK.csv | tail -n 10 | cut -f 1 -d',' > batches-ancestry-small/batchAFR.txt
head -n 11 EUR_WHITE.csv | tail -n 10 | cut -f 1 -d',' > batches-ancestry-small/batchEUR.txt

chrom=21; ./gnomix_launcher.py \
  --name test-chr21 \
  --vcf ${WGS_ACAF_THRESHOLD_VCF_PATH}/acaf_threshold.chr${chrom}.vcf.bgz \
  --chrom ${chrom} \
  --sample-batches batches-ancestry-small/ --max-batches 2 #--extra-subset-args "-r chr21:1-5252900"

# Bigger test
mkdir batches-ancestry-big
head -n 1001 AFR_BLACK.csv | tail -n 1000 | cut -f 1 -d',' > batches-ancestry-big/batchAFR.txt
head -n 1001 EUR_WHITE.csv | tail -n 1000 | cut -f 1 -d',' > batches-ancestry-big/batchEUR.txt
chrom=21; ./gnomix_launcher.py \
  --name test-chr21-big \
  --vcf ${WGS_ACAF_THRESHOLD_VCF_PATH}/acaf_threshold.chr${chrom}.vcf.bgz \
  --chrom ${chrom} \
  --sample-batches batches-ancestry-big/ --max-batches 2
```

