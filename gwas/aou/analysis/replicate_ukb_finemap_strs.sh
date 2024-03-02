#!/bin/bash

# TODO - add ability to do outprefix 

# TODO update with full UKB finemapped results
VCF=${WORKSPACE_BUCKET}/cromwell-execution/targetTR/4519d903-dde6-4c8e-b1ee-bcc2d7cd6dd7/call-sort_index/CBL_test.filtered.sorted.vcf.gz

while IFS= read -r line
do
	phenotype=$(echo $line | cut -d',' -f 1)
  	if [[ $phenotype == "phenotype" ]]; then
     	continue
  	fi
  	echo $phenotype
  	cd ..
  	for cohort in passing_samples_v7 EUR_WHITE AFR_BLACK
  	do
  		echo "  ${cohort}"
  		./aou_gwas.py --phenotype ${phenotype} --num-pcs 10 \
  			--method associaTR --tr-vcf $VCF --norm quantile \
  			--samples ${WORKSPACE_BUCKET}/samples/${cohort}.csv
  	done
  	cd -
  	mv ../${phenotype}_associaTR_*.gwas.tab .
done < ../phenotypes_manifest.csv 