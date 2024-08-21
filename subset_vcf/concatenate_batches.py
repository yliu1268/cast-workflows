#!/usr/bin/env python3
"""
Usage:

./concanate_batches.py <cromshell_job_output.json> <fileprefix>
"""

import json
import os
import sys

jobdata = json.load(open(sys.argv[1], "r"))
outprefix = sys.argv[2]

# Gather files for each batch
batch_files = {} # batchname -> {"vcf": [], "index": []}
for region in jobdata["subset_vcf.vcf_output_array"]:
	for f in region:
		batchname = os.path.basename(f).split("-")[0]
		if batchname not in batch_files.keys():
			batch_files[batchname] = {}
			batch_files[batchname]["vcf"] = []
			batch_files[batchname]["index"] = []
		batch_files[batchname]["vcf"].append(f)
for region in jobdata["subset_vcf.vcf_index_array"]:
	for f in region:
		batchname = os.path.basename(f).split("-")[0]
		batch_files[batchname]["index"].append(f)

# Process one batch at a time
for batch in batch_files.keys():
	vcf_files = batch_files[batch]["vcf"]
	index_files = batch_files[batch]["index"]
	print("##### Processing %s ######"%batch)
	print(vcf_files)
	print(index_files)