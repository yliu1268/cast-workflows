# Setup

This must be run on the All of Us workbench. You must have launched both an analysis and cromshell environment. This also requires having run the `subset_vcf` workflow to generate batch VCFs at `${WORKSPACE_BUCKET}/acaf_batches/chr${chrom}`.

```
git clone https://github.com/cast-genomics/cast-workflows
cd cast-workflows/
git checkout gnomix
cd local_ancestry
../utils/configure-cromshell.py
```

# Test

```
chrom=11; ./gnomix_launcher.py \
  --name test-chr${chrom} \
  --vcfdir ${WORKSPACE_BUCKET}/acaf_batches/chr${chrom} \
  --chrom ${chrom} \
  --max-batches 2 
```

# Full chromosome

```
chrom=11; ./gnomix_launcher.py \
  --name gnomix-chr${chrom} \
  --vcfdir ${WORKSPACE_BUCKET}/acaf_batches/chr${chrom} \
  --chrom ${chrom}
```
