#!/bin/bash

# Define the base URL where files are located
vcf="https://ensemble-tr.s3.us-east-2.amazonaws.com/additional-phased-trs/chr"
gcs="$WORKSPACE_BUCKET/tr_imputation/tr_imputation/ref/"
# Loop through chromosomes 1 to 21
for (( chr=1; chr<=21; chr++ ))
do
    # Construct the URL for the current chromosome
    vcf="${vcf}${chr}_final_SNP_merged_additional_TRs.vcf.gz"
    vcf_index="${vcf}.tbi"
    
    # Download the file using wget
    echo "Downloading ${vcf} ..."
    wget "${vcf}"
    wget "${vcf_index}"

    echo "Move ref vcf files to google cloud"
    gsutil cp "chr${chr}_final_SNP_merged_additional_TRs.vcf.gz" $gcs
done

echo "Download complete."

