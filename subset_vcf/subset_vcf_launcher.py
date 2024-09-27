#!/usr/bin/env python3
"""
Script to luanch jobs to subset AOU VCFs

To run a small test:
./subset_vcf_launcher.py --test --chrom 11
"""

import argparse
import json
import os
import subprocess

def GetRegions(chrom):
	regions_file = os.environ["WORKSPACE_BUCKET"] + "/subset_vcf/metadata/regions_chr%s.txt"%chrom
	os.system("gsutil cp %s ."%regions_file)
	regions = [item.strip() for item in open(os.path.basename(regions_file), "r").readlines()]
	return regions

def main():
	parser = argparse.ArgumentParser(__doc__)
	parser.add_argument("--chrom", help="Which chromosome to process. Ignored if testing", type=str, required=True)
	parser.add_argument("--test", help="Whether to run on small test sets", action="store_true")
	parser.add_argument("--name", help="Name prefix for intermediate files", type=str, default="test")
	parser.add_argument("--disk", help="Amount of diskspace (in GB) for subset region batches. CURRENTLY IGNORED", default=30, type=int)
	args = parser.parse_args()

	# Set up workflow JSON
	json_dict = {}
	json_dict["subset_vcf.multi_sample_vcf"] = os.environ["WGS_ACAF_THRESHOLD_VCF_PATH"] + "/acaf_threshold.chr%s.vcf.bgz"%args.chrom
	json_dict["subset_vcf.GOOGLE_PROJECT"] = os.environ.get("GOOGLE_PROJECT", "")
	json_dict["subset_vcf.disk"] = args.disk
	if args.test:
		json_dict["subset_vcf.sample_groups"] = os.environ["WORKSPACE_BUCKET"] + "/subset_vcf/metadata/aou_sample_groups.txt"
		json_dict["subset_vcf.regions"] = ["chr%s:0-10000000"%args.chrom, "chr%s:10000000-20000000"%args.chrom]
	else:
		json_dict["subset_vcf.sample_groups"] = os.environ["WORKSPACE_BUCKET"] + "/subset_vcf/metadata/aou_sample_groups.txt"
		json_dict["subset_vcf.regions"] = GetRegions(args.chrom)

	# Convert to JSON and save to a file
	json_file = args.name + ".aou.json"
	with open(json_file, "w") as f:
		json.dump(json_dict, f, indent=4)

	# Run workflow
	cmd = "cromshell submit {wdl} {json}".format(wdl="wdl/subset_vcf.wdl", json=json_file)
	output = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).stdout.read()
	print(output.decode("utf-8"))	

if __name__ == "__main__":
	main()