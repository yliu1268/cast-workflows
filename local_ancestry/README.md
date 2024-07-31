# Setup

```
git clone https://github.com/cast-genomics/cast-workflows
cd cast-workflows/
git checkout gnomix
cd local_ancestry
../utils/configure_cromshell.py
```

# Test

```
chrom=21; ./gnomix_launcher.py \
  --name test \
  --vcf ${WGS_ACAF_THRESHOLD_VCF_PATH}/acaf_threshold.chr${chrom}.vcf.bgz \
  --chrom ${chrom} \
  --dryrun \
  --sample-batches ./tests/dummy_sample_batches
```