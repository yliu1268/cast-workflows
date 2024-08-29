#!/usr/bin/env python3
"""
Script to launch AOU imputation use new ref panel 

example code to impute 10 samples at CBL region 

./run_subset_vcf.py \
--name batch_test 
--vcf gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/acaf_threshold_v7.1/vcf/acaf_threshold.chr21.vcf.bgz \
--batch-num 2
"""


import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
import csv
from utils import MSG, ERROR

def GetFileBatches(sample_list, batch_num=None):
    sample_batch = []
    
    # Open txt file
    with open(sample_list, "r") as f:
        if batch_num is None:
            # Read all lines if batch_num is None
            sample_batch = f.readlines()
        else:
            # Read the first batch_num lines
            for _ in range(batch_num):
                line = f.readline().strip()  # Read and strip newline character
                if line:  # Check if the line is not empty
                    sample_batch.append(line)
                else:
                    break  # Exit the loop if there are no more lines
    
    return sample_batch

def RunWorkflow(json_file, json_options_file, dryrun=False):
	"""
	Run workflow on AoU

	Arguments
	---------
	json_file : str
		JSON file path with input arguments
	json_options_file : str
		JSON with additional options for cromshell

	dryrun : bool
		Just print the command, don't actually run cromshell
	"""
	
	cmd = "cromshell submit ../wdl/subset_vcf.wdl {json} -op {options}".format(json=json_file, options=json_options_file)

	if dryrun:
		sys.stderr.write("Run: %s\n"%cmd)
		return
	output = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).stdout.read()
	print(output.decode("utf-8"))


def DownloadGS(filename):
	"""
	Download a GCP path locally

	Arguments
	---------
	filename : str
	   GCP path
	"""
	cmd = "gsutil cp {filename} .".format(filename=filename)
	output = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).stdout.read()
	print(output.decode("utf-8"))	
	
def UploadGS(local_path, gcp_path):
	"""
	Upload a local file to GCP

	Arguments
	---------
	local_path : str
	   Local path
	gcp_path : str
	   GCP path to upload to
	"""
	cmd = "gsutil cp {src} {dest}".format(src=local_path, dest=gcp_path)
	output = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).stdout.read()
	print(output.decode("utf-8"))	


def main():
	# Get token
	token_fetch_command = subprocess.run(['gcloud', 'auth', 'application-default', 'print-access-token'], \
		capture_output=True, check=True, encoding='utf-8')
	token = str.strip(token_fetch_command.stdout)

	# Set up output bucket
	bucket = os.getenv("WORKSPACE_BUCKET")
	project = os.getenv("GOOGLE_PROJECT")
	

	parser = argparse.ArgumentParser(__doc__)
	parser.add_argument("--name", help="Name of the TR job", required=True, type=str)
	parser.add_argument("--vcf", help="Name of the genotype vcf file", required=True, type=str)
	parser.add_argument("--samples", help="List of samples to process ", type=str, required=False, \
					 default=f"{bucket}/tr_imputation/tr_imputation/sample/sample_manifest.txt")
	parser.add_argument("--dryrun", help="Don't actually run the workflow. Just set up", action="store_true")
	parser.add_argument("--batch-num", help="Number of batches. Default: -1 (all)",type=int, required=False, default=None)

	args = parser.parse_args()

	output_bucket = bucket + "/" + args.name

	#Set up sample file list
	if args.samples.startswith("gs://"):
		DownloadGS(args.samples)
		sample_list = os.path.basename(args.samples)
	else: sample_list = args.samples

	# Upload vcf file
	if args.vcf.startswith("gs://"):
		vcf_gcs = args.vcf
	else:
				# Copying the vcf file
		vcf_gcs = output_bucket + "/" + args.name + "/"
		UploadGS(args.vcf, vcf_gcs)
				# Copying the index file
		UploadGS(args.vcf + ".tbi", vcf_gcs)

	# set up batches of file
	sample_batch = GetFileBatches(sample_list,args.batch_num)


	# Set up workflow JSON
	json_dict = {}
	json_dict["run_subset.vcf"] = args.vcf
	json_dict["run_subset.vcf_index"]=args.vcf+".tbi"
	json_dict["run_subset.out_prefix"] = args.name
	json_dict["run_subset.GOOGLE_PROJECT"] = project
	json_dict["run_subset.samples"] = sample_batch 
	

	# Convert to json and save as a file
	json_file = args.name+".aou.json"
	with open(json_file, "w") as f:
		json.dump(json_dict, f, indent=4)


	# Set up json optionsß
	json_options_dict = {}
	json_options_file = args.name+".options.aou.json"
	with open(json_options_file, "w") as f:
		json.dump(json_options_dict, f, indent=4)


	# Run workflow on AoU using cromwell
	RunWorkflow(json_file, json_options_file,dryrun=args.dryrun)


if __name__ == "__main__":
	main()
