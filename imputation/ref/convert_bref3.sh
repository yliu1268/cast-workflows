#!/bin/bash

#set up gc bucket storage"
gcs="$WORKSPACE_BUCKET/tr_imputation/tr_imputation/ref/"

#dowaload bref3.jar
wget https://faculty.washington.edu/browning/beagle/bref3.27May24.118.jar 


#convert vcf to bref3
bref="bref3.27May24.118.jar"
for (( chr=1; chr<=21; chr++ ))
do 
    echo "converting vcf to bref3 format"
    zcat ${chr}_final_SNP_merged_additional_TRs.vcf.gz | java -jar $bref > ${chr}_final_SNP_merged_additional_TRs.bref3

    echo "Done converting"

    gsutil cp ${chr}_final_SNP_merged_additional_TRs.bref3 $gcs
    echo "Successfully transfer bref to google cloud bucket"
done 



