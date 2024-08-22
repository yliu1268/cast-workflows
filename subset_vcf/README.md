# Pre-compute batches of ACAF VCF files for use in future workflows

This set of workflows precomputes small subsets of the ACAF VCF files. The steps below:

Setup work on the terminal:

* Break up samples into batches of 1000
* Sets up the "groups" file needed for subsetting with bcftools split (`aou_sample_groups.txt`) 
* Sets up 10Mb regions to run on one at a time

Workflows launched with cromshell:

* `subset_vcf.wdl`: Launch a WDL that creates 1 smaller VCF per sample batch per region
* `concatenate_batch_vcfs.wdl`: Oragnize all the subset VCF files per batch in a single folder in our workspace for future use

## Setup: Set up group batches

```
# Set up groups file for batches of 1000
gsutil cp $WORKSPACE_BUCKET/tr_imputation/tr_imputation/sample/aou_sample_list.txt .
mkdir batches
split -l 1000 aou_sample_list.txt batches/batch
counter=1
for file in batches/batch*
do
    mv "$file" "batches/batch$counter"
    ((counter++))
done

# Make full group set
for file in batches/batch*
do
	batchname=$(basename $file)
	cat $file | awk -v"prefix=$batchname" '{print $1 "\t-\t"prefix}'
done > aou_sample_groups.txt
gsutil cp aou_sample_groups.txt ${WORKSPACE_BUCKET}/subset_vcf/metadata/aou_sample_groups.txt

# Make small group set
cat aou_sample_groups.txt | awk '($3=="batch1" || $3=="batch2")' > aou_sample_groups_test.txt
gsutil cp aou_sample_groups_test.txt ${WORKSPACE_BUCKET}/subset_vcf/metadata/aou_sample_groups_test.txt
```

## Set up: Set up region batches
```
bedtools makewindows -g hg38.txt -w 10000000 > genome_windows.bed
for chrom in $(seq 1 22)
do
	cat genome_windows.bed | grep -w "chr"${chrom} | awk '{print $1 ":" $2 "-"$3}' > regions_chr${chrom}.txt
	gsutil cp regions_chr${chrom}.txt ${WORKSPACE_BUCKET}/subset_vcf/metadata/regions_chr${chrom}.txt
done
```

## Subset VCF: Launch jobs

```
# Test
./subset_vcf_launcher.py --test --chrom 11
```

```
# Full jobs
for chrom in $(seq 1 22)
do
	./subset_vcf_launcher.py --chrom ${chrom} --name chr${chrom}
done
```

## Concatenation: Organize files for each batch

```
./concatenate_batches_v2.py $jobid chr11
```

TODO - changing this to wdl instead (see concate_batches_v2.py - in progress)

Note assumes we compiled bcftools with ability to read from gcs:
```
# install bcftools
wget -O bcftools-1.20.tar.bz2 https://github.com/samtools/bcftools/releases/download/1.20/bcftools-1.20.tar.bz2
tar -xjvf bcftools-1.20.tar.bz2
cd bcftools-1.20
./configure --enable-gcs --enable-libcurl --enable-s3
make install
make -j

# Set env variables
export GCS_REQUESTER_PAYS_PROJECT=${GOOGLE_PROJECT}
export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
```

```

# Notes for doing this for test on chr11
# But ugh this is kind of slow should we move this to wdl to parallelize?
# Also will need to keep track of the jobid for each chrom from subset_vcf_launcher.py above
chrom=11
cromshell --machine_processable list-outputs -j $jobid > chr${chrom}.json
./concatenate_batches.py chr${chrom}.json chr${chrom} ~/bcftools-1.20/bcftools
```

