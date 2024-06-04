#!/bin/bash
VCF=$1
SAMPFILE=$2
OUTPREFIX=$3

# Get column subset to extract
sampcols=$(echo $(bcftools query -l $VCF | awk '{print 9+NR "\t" $1}' | grep -w -f $SAMPFILE | cut -f 1 ) | tr ' ' ,)
cols=$(echo "1-9,${sampcols}") # always need 1-9. samples start at 10

# Get header
bcftools view -h $VCF | grep "^##" > ${OUTPREFIX}.vcf
bcftools view -h $VCF | grep "^#CHROM" | cut -f${cols} >> ${OUTPREFIX}.vcf

# Get rest of data
# TODO: use curl to get VCF contents
zcat $VCF | grep -v "^#" | cut -f${cols} >> ${OUTPREFIX}.vcf

# Compress and index
bgzip ${OUTPREFIX}.vcf
tabix -p vcf ${OUTPREFIX}.vcf.gz
