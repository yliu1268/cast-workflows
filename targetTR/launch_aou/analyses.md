# TargetTR analysis runs on AoU

## STR GWAS replication analysis

Run HipSTR genotyping on all TRs of interest from [Margoliash et al.](https://www.cell.com/cell-genomics/pdfExtended/S2666-979X(23)00302-6).

```
# Run on one batch
batch=1 # 1..10
./targetTR_launcher_aou.py \
  --tr-bed ../strsets/ukb_finemapped_hg38/ukb_finemapped_hg38_batch${batch}.bed \
  --name UKB_finemapped-batch${batch} \
  --separate-hipstr-runs \
  --action run-batches 
```

# Outputs

You should copy outputs of each job (.vcf.gz and .vcf.gz.tbi) to:
```
${WORKSPACE_BUCKET}/ukb_finemapped_hg38/
```
