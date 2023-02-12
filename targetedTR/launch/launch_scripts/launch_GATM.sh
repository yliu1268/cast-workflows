#!/bin/bash

../targetTR_launcher_ukb.py \
  --region chr15:45396869-45396892 \
  --period 1 \
  --refcopies 23 \
  --name GATM-Homopolymer \
  --batch-size 3 \
  --batch-num 3 \
  --max-batches-per-workflow 1 \
  --workflow-id workflow-GPYY0p0Jv7B730qb33V6jvG0 \
  --file-list ../ukb_cram_and_index_files.txt