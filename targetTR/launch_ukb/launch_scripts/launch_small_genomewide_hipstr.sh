#!/bin/bash

# Run HipSTR genome-wide on 20 samples
# file-GfbbqJQJv7B9qQ8Ky7Jz2QkG has hg38.hipstr_reference.bed

../targetTR_launcher_ukb.py \
  --tr-bed file-GfbbqJQJv7B9qQ8Ky7Jz2QkG \
  --name hipstr-gw-small \
  --file-list ../ukb_cram_and_index_files.txt \
  --batch-size 20 \
  --batch-num 1 \
  --max-batches-per-workflow -1 \
  --concurrent