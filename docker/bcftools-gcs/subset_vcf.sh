#!/bin/bash
VCF=$1
SAMPFILE=$2
# Get column subset to extract
sampcols=$(echo $(bcftools query -l $VCF | awk '{print 9+NR "\t" $1}' | grep -w -f $SAMPFILE | cut -f 1 ) | tr ' ' ,)
cols=$(echo "1-9,${sampcols}") # always need 1-9. samples start at 10
# Get header
bcftools view -h $VCF | grep "^##" > subset_custom.vcf
bcftools view -h $VCF | grep "^#CHROM" | cut -f${cols} >> subset_custom.vcf
# Get rest of data
zcat $VCF | grep -v "^#" | cut -f${cols} >> subset_custom.vcf
# Use vcf-annotate to update AC/AN. Not necessary for Beagle but
# vcf-validator complains. Maybe can remove this step if it is long
#cat subset_custom.vcf | vcf-annotate --fill-AC-AN | bgzip -c > subset_custom.vcf.gz
# Update: skip annotate
bgzip subset_custom.vcf