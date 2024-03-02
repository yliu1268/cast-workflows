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

# compile results
cat urea_associaTR_ALL.gwas.tab | grep -v "^#" | head -n 1 | awk '{print "trait\tcohort\t" $0}'> all_traits_associaTR.tab
for resfile in $(ls *gwas.tab)
do
	phenotype=$(echo $resfile | sed 's/_associaTR_ALL.gwas.tab//' | sed 's/_associaTR_AFR_BLACK.gwas.tab//' | sed 's/_associaTR_EUR_WHITE.gwas.tab//')
	cohort=$(echo $resfile | awk -F"_" '{print $NF}' | sed 's/.gwas.tab//')
	cat $resfile | grep -v "^#" | grep -v "^chrom" | awk -v"trait=$phenotype" -v"cohort=$cohort" '{print trait "\t" cohort "\t" $0}' >> all_traits_associaTR.tab
done