# Pre-compute batches of ACAF VCF files for use in future workflows

This workflow precomputes small subsets of the ACAF VCF files. The steps below:

* Break up samples into batches of 1000
* Sets up the "groups" file needed for subsetting with bcftools split (`aou_sample_groups.txt`) 
* Sets up 10Mb regions to run on one at a time
* Launch a WDL that creates 1 smaller VCF per sample batch per region
* Oragnize all the subset VCF files per batch in a single folder in our workspace for future use

TODO:
* maybe can reduce disk space, based on output vcf size * num batches (~70Mx246 batches=~17GB? so maybe need closer to 25GB?)

## Set up group batches

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

## Set up region batches
```
bedtools makewindows -g hg38.txt -w 10000000 > genome_windows.bed
for chrom in $(seq 1 22)
do
	cat genome_windows.bed | grep -w "chr"${chrom} | awk '{print $1 ":" $2 "-"$3}' > regions_chr${chrom}.txt
	gsutil cp regions_chr${chrom}.txt ${WORKSPACE_BUCKET}/subset_vcf/metadata/regions_chr${chrom}.txt
done
```

## Launch jobs

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

## Organize files for each batch - TODO

Might keep in separate regions per batch, jobs can concatenate them downstream?