# TargetTR analysis runs on AoU

## STR GWAS replication analysis

Run HipSTR genotyping on all TRs of interest from [Margoliash et al.](https://www.cell.com/cell-genomics/pdfExtended/S2666-979X(23)00302-6).

```
./targetTR_launcher_aou.py \
  --tr-bed ../strsets/ukb_finemapped_hg38.bed \
  --name UKB_finemapped \
  --separate-hipstr-runs \
  --action run-batches 
```

### Test runs
```
# Testing can remove when done
./targetTR_launcher_aou.py \
  --tr-bed ../strsets/test_10.bed \
  --name UKB_test_allsamples_10TRs \
  --separate-hipstr-runs \
  --action run-batches 

./targetTR_launcher_aou.py \
  --tr-bed ../strsets/ukb_all_hg38.bed \
  --name UKB_test_2batches_allTRs \
  --separate-hipstr-runs \
  --batch-num 2 \
  --action run-batches
```