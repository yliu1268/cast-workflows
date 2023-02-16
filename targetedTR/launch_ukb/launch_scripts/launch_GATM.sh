#!/bin/bash

../targetTR_launcher_ukb.py \
  --region chr15:45396869-45396892 \
  --period 1 \
  --refcopies 23 \
  --name GATM-Homopolymer \
  --batch-size 500 \
  --max-batches-per-workflow 10 \
  --concurrent \
  --workflow-id workflow-GPfbXV8Jv7B27kpf6Y50QyQ9 \
  --file-list ../ukb_cram_and_index_files.txt

# TODO remove batch-num

# Notes:
# batch size 500 is good
# up to 10 batches on a single job worked. 25 failed
# for 200K crams: need:
#  --batch-size 500 --max-batches-per-workflow 10 --concurrent
