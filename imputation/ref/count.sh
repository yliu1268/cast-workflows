#!/bin/bash
ref="ref_panel_summary.txt"
echo "chr,snp,tr" | sed 's/,/\t/g' > ${ref}

for chr in {chr1..chr22}; do
    vcf="${chr}_final_SNP_merged_additional_TRs.vcf.gz"
    count_snp=$(bcftools view -i 'ID!="."' $vcf |grep -v "^#" | wc -l)
    count_tr=$(bcftools view -i 'ID="."' $vcf |grep -v "^#" | wc -l)

    echo -e "${chr}\t${count_snp}\t${count_tr}" >> ${ref}

done