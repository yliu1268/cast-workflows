#!/bin/bash

vcf="$WORKSPACE_BUCKET/tr_imputation/tr_imputation/ensemble_chr21_filtered.vcf.gz"
#local_vcf="data/ensemble_chr21_filtered.vcf.gz"
python imputation_aou.py \
	--name tr_imputation \
	--vcf $vcf \
	--ref-panel $vcf
