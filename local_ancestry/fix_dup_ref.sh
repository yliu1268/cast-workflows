#!/bin/bash

# Fix 1kg ref panel for chroms that had
# duplicate markers

chrom=17

gsutil cp ${WORKSPACE_BUCKET}/gnomix/resources/ALL.chr${chrom}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz .

bcftools norm -d all -Oz \
	-o ALL.chr${chrom}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes-dedup.vcf.gz \
	ALL.chr${chrom}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
tabix -p vcf ALL.chr${chrom}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes-dedup.vcf.gz