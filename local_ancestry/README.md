# Setup

This must be run on the All of Us workbench. You must have launched both an analysis and cromshell environment.

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
  --name test-chr11 \
  --vcfdir ${WORKSPACE_BUCKET}/acaf_batches/chr${chrom} \
  --chrom ${chrom} \
  --max-batches 2 
```

