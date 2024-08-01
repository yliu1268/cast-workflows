#!/usr/bin/env python3
"""
Script to launch AOU Gnomix for local ancestry inference

Dryrun:
chrom=21; ./gnomix_launcher.py \
  --name test \
  --vcf ${WGS_ACAF_THRESHOLD_VCF_PATH}/acaf_threshold.chr${chrom}.vcf.bgz \
  --chrom ${chrom} \
  --dryrun \
  --sample-batches ./tests/dummy_sample_batches
"""

import argparse
import glob
import json
import os
import sys

sys.path.append("../utils")
import aou_utils

GNOMIXMODEL = "gs://artifacts.ucsd-medicine-cast.appspot.com/resources/pretrained_gnomix_models.tar.gz"

def main():
	parser = argparse.ArgumentParser(__doc__)
	parser.add_argument("--sample-batches", help="Directory with sample batches", \
		type=str, required=True)
	parser.add_argument("--name", help="Name of the run", type=str, required=True)
	parser.add_argument("--vcf", help="GCP path to the AoU SNP VCF", type=str, required=True)
	parser.add_argument("--chrom", help="Which chromosome to process", type=str, required=True)
	parser.add_argument("--gnomix-model", help="GCP path to pretrained model", type=str, \
		default=GNOMIXMODEL)
	# For debugging
	parser.add_argument("--dryrun", help="Set up job but don't launch", \
		action="store_true")
	parser.add_argument("--max-batches", help="Only run this many batches", \
		type=int, default=-1)
	args = parser.parse_args()

	# Get list of sample files and upload to GCS
	gsprefix = os.getenv("WORKSPACE_BUCKET") + "/gnomix/" + args.name + "/batches"
	sample_file_list = glob.glob(args.sample_batches + "/*")
	if args.max_batches > -1 and len(sample_file_list) > args.max_batches:
		sample_file_list = sample_file_list[:args.max_batches]
	sample_file_list_gcs = []
	for sf in sample_file_list:
		sf_gcs = gsprefix + "/" + os.path.basename(sf)
		UploadGS(sf, sf_gcs)
		sample_file_list_gcs.append(sf_gcs)

	# Set up workflow JSON
	json_dict = {}
	json_dict["local_ancestry.out_prefix"] = args.name
	json_dict["local_ancestry.multi_sample_vcf"] = args.vcf
	json_dict["local_ancestry.samples"] = sample_file_list_gcs
	json_dict["local_ancestry.chrom"] = args.chrom
	json_dict["local_ancestry.GOOGLE_PROJECT"] = os.environ.get("GOOGLE_PROJECT", "")
	json_dict["local_ancestry.gnomix_model"] = args.gnomix_model

	# Convert to JSON and save to a file
	json_file = args.name + ".aou.json"
	with open(json_file, "w") as f:
		json.dump(json_dict, f, indent=4)

	# Set up json options
	json_options_dict = {}
	json_options_file = args.name+".options.aou.json"
	with open(json_options_file, "w") as f:
		json.dump(json_options_dict, f, indent=4)

	# Zip all the WDL depencies
	wdl_dependencies_file = args.name + "-wdl.zip"
	aou_utils.ZipWDL(wdl_dependencies_file, ["wdl/gnomix.wdl"])

	# Run workflow on AoU using cromshell
	aou_utils.RunWorkflow("wdl/local_ancestry.wdl", json_file, json_options_file, \
		wdl_dependencies_file, dryrun=args.dryrun)

if __name__ == "__main__":
	main()