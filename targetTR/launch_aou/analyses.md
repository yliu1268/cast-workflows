# TargetTR analysis runs on AoU

## STR GWAS replication analysis

Run HipSTR genotyping on all TRs of interest from [Margoliash et al.](https://www.cell.com/cell-genomics/pdfExtended/S2666-979X(23)00302-6).

```
./targetTR_launcher_aou.py \
  --tr-bed ../strsets/ukb_all_hg38.bed \
  --name UKB_ALL \
  --action run-batches 
```