# Setup

This must be run on the All of Us workbench. You must have launched both an analysis and cromshell environment.

```
git clone https://github.com/cast-genomics/cast-workflows
cd cast-workflows/
git checkout gnomix
cd local_ancestry
../utils/configure-cromshell.py
```

Get list of SNPs for each chrom

```
for chrom in $(seq 1 22)
do
  cat pretrained_gnomix_models/chr${chrom}/chr${chrom}.bim | awk '{print $1 "\t" $4-1 "\t" $4}' > snplist_chr${chrom}.txt
  gsutil cp snplist_chr${chrom}.txt ${WORKSPACE_BUCKET}/gnomix/resources/snplist_chr${chrom}.txt
done
```

# Test

```
chrom=11; ./gnomix_launcher.py \
  --name test-chr${chrom} \
  --vcfdir ${WORKSPACE_BUCKET}/acaf_batches/chr${chrom} \
  --chrom ${chrom} \
  --max-batches 2 
```

