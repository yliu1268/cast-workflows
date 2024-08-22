# Pre-compute batches of ACAF VCF files for use in future workflows

This set of workflows precomputes small subsets of the ACAF VCF files. The steps below:

Setup work on the terminal:

* Break up samples into batches of 1000
* Sets up the "groups" file needed for subsetting with bcftools split (`aou_sample_groups.txt`) 
* Sets up 10Mb regions to run on one at a time

Workflows launched with cromshell:

* `subset_vcf.wdl`: Launch a WDL that creates 1 smaller VCF per sample batch per region
* `concatenate_batch_vcfs.wdl`: Oragnize all the subset VCF files per batch in a single folder in our workspace for future use

## Setup

### Set up group batches

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

### Set up region batches
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
# Full jobs
# Keep track of the job ID of each chromosome, needed for next steps below
for chrom in $(seq 1 22)
do
	./subset_vcf_launcher.py --chrom ${chrom} --name chr${chrom}
done
```

## Concatenation: Organize files for each batch

```
# Launch concatenation wdl ( example on chr11)
chrom=chr11
jobid=<jobid from the subset_vcf job for this chromosome>
./concatenate_batches_v2.py ${jobid} ${chrom}

# After completion, rename job outputs to be in ${WORKSPACE_BUCKET}/acaf_batches/${chrom}
# cromshell buries them in ${WORKSPACE_BUCKET}/acaf_batches/${chrom}/concatenate_batch_vcfs/
concatjobid=<job id from concatenation>
cromshell -mc list-outputs -j ${concatjobid} | \
	python -c "import json, sys; data=json.load(sys.stdin); [sys.stdout.write(item+'\n') for item in data['concatenate_batch_vcfs.vcf_outputs']+data['concatenate_batch_vcfs.vcf_indices']]" | \
	xargs -n1 -I% -P1 sh -c "gsutil mv % ${WORKSPACE_BUCKET}/acaf_batches/${chrom}/"
```