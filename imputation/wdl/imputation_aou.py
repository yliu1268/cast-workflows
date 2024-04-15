#!/usr/bin/env python3
"""
Script to launch AOU imputation use new ref panel 
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile


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
	cmd = "cromshell submit ../wdl/beagle.wdl {json} -op {options}".format(json=json_file, options=json_options_file)
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
	cmd = "gsutil -u $GOOGLE_PROJECT cp {filename} .".format(filename=filename)
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
	parser = argparse.ArgumentParser(__doc__)
	parser.add_argument("--name", help="Name of the TR job", required=True, type=str)
	parser.add_argument("--vcf", help="Name of the genotype vcf file", required=True, type=str)
	parser.add_argument("--ref-genome", help="File id of ref genome", type=str, default="https://ensemble-tr.s3.us-east-2.amazonaws.com/additional-phased-trs/chr21_final_SNP_merged_additional_TRs.vcf.gz")
	parser.add_argument("--ref-genome-tbi", help="File id of ref genome index", type=str, default="https://ensemble-tr.s3.us-east-2.amazonaws.com/additional-phased-trs/chr21_final_SNP_merged_additional_TRs.vcf.gz.tbi")
	parser.add_argument("--action", help="Options: create-batches, run-batches, both", type=str, required=True)
	parser.add_argument("--dryrun", help="Don't actually run the workflow. Just set up", action="store_true")

	args = parser.parse_args()


    # Get token
	token_fetch_command = subprocess.run(['gcloud', 'auth', 'application-default', 'print-access-token'], \
		capture_output=True, check=True, encoding='utf-8')
	token = str.strip(token_fetch_command.stdout)

	# Set up output bucket
	bucket = os.getenv("WORKSPACE_BUCKET")
	project = os.getenv("GOOGLE_PROJECT")
	output_bucket = bucket + "/" + args.name

    # Set up workflow JSON
	json_dict = {}
	json_dict["targetTR.genome"] = args.genome_id
	json_dict["targetTR.genome_index"] = args.genome_idx_id
	json_dict["targetTR.outprefix"] = args.name
	json_dict["targetTR.GOOGLE_PROJECT"] = project
	json_dict["targetTR.GCS_OAUTH_TOKEN"] = token


	# Convert to json and save as a file
	json_file = args.name+".aou.json"
	with open(json_file, "w") as f:
		json.dump(json_dict, f, indent=4)

	# Set up json options
	json_options_dict = {}
	json_options_file = args.name+".options.aou.json"
	with open(json_options_file, "w") as f:
		json.dump(json_options_dict, f, indent=4)

	# Run workflow on AoU using cromwell
	RunWorkflow(json_file, json_options_file, dryrun=args.dryrun)


if __name__ == "__main__":
main()