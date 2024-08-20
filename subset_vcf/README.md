# Pre-compute batches of ACAF VCF files for use in future workflows

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
# TODO - precompute regions per chromosome once we figure out a good size
```

## Launch jobs