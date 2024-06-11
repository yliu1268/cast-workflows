#!/usr/bin/env python3
"""
Script to launch AOU imputation use new ref panel 

example code to impute 10 samples at CBL region 
./imputation_aou.py \
--name test_imputation 
--vcf gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/acaf_threshold_v7.1/vcf/acaf_threshold.chr21.vcf.bgz \
--ref-panel $WORKSPACE_BUCKET/tr_imputation/tr_imputation/ref/chr21_final_SNP_merged_additional_TRs.bref3 \
--samples-file $WORKSPACE_BUCKET/tr_imputation/tr_imputation/subset_samples/aou_subset_samples_100.txt \
--regions-file $WORKSPACE_BUCKET/tr_imputation/tr_imputation/sub_region5mb.bed \
--beagle_region
--chrom chr21:5101889-6101889
--mem 40
"""


import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
import tqdm
from utils import MSG, ERROR


def GetVCFBatches(file_list, batch_size, batch_num=-1, gsprefix=None):
	with open(file_list, 'r') as f:
		lines = f.readlines()

	num_lines = len(lines)
	batch_num = num_lines // batch_size+ (1 if num_lines % batch_size != 0 else 0)

	for i in range(batch_num):
		start = i * batch_size
		end = min((i + 1) * batch_size, num_lines)
		batch_file = f"sample_batch{i}.txt"
		with open(batch_file, 'w') as f_out:
			f_out.writelines(lines[start:end])
	
	print(f"Split into {batch_num} batches.")
	return [ f"sample_batch{i}.txt" for i in range(batch_num)]



def RunWorkflow(json_file, json_options_file, cromwell, dryrun=False):
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
	if cromwell is False:
		cmd = "cromshell submit ../wdl/imputation.wdl {json} -op {options}".format(json=json_file, options=json_options_file)
	else:
		cmd = "java -jar -Dconfig.file={} ".format("/home/jupyter/cromwell.conf") + \
	  			"cromwell-87.jar run imputation.wdl " + \
	  			"--inputs {} --options {}".format(json_file, json_options_file)
	if dryrun:
		sys.stderr.write("Run: %s\n"%cmd)
		return
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
	parser.add_argument("--ref-panel", help="File id of ref genome", type=str)
	parser.add_argument("--mem", help="Specify run memory ", type=int, required=False, default=32)
	parser.add_argument("--window", help="Specify window size for imputation ", type=int, required=False, default=20)
	parser.add_argument("--samples", help="Name of sub_samples file ", type=str, required=True)
	parser.add_argument("--region", help="Name of chrom position  chr:xxx-xxx", type=str,required=False)
	parser.add_argument("--beagle-region", help="Apply chrom for beagle", action="store_true",required=False)
	parser.add_argument("--subset-region", help="Subsetting region for vcf file", action="store_true",required=False)
	parser.add_argument("--dryrun", help="Don't actually run the workflow. Just set up", action="store_true")
	parser.add_argument("--file-list", help="List of sample to process"
		"Format: person_id", type=str, required=False, \
		default="$WORKSPACE_BUCKET/tr_imputation/tr_imputation/sample/aou_sample_list.txt")
	parser.add_argument("--batch-size", help="Sample batch size", required=False, type=int, default=1000)
	parser.add_argument("--batch-num", help="Number of batches. Default: -1 (all)", required=False, default=-1)
	parser.add_argument("--cromwell", help="Run using cormwell as opposed to the default cromshell",
                            action="store_true", default=False)
	#parser.add_argument("--action", help="Options: create-batches, run-batches, both", type=str, required=True)
	


	args = parser.parse_args()

	
	# Get token
	token_fetch_command = subprocess.run(['gcloud', 'auth', 'application-default', 'print-access-token'], \
		capture_output=True, check=True, encoding='utf-8')
	token = str.strip(token_fetch_command.stdout)

	# Set up output bucket
	bucket = os.getenv("WORKSPACE_BUCKET")
	project = os.getenv("GOOGLE_PROJECT")
	output_bucket = bucket + "/" + args.name


	# Upload vcf file
	if args.vcf.startswith("gs://"):
		vcf_gcs = args.vcf
	else:
				# Copying the vcf file
		vcf_gcs = output_bucket + "/" + args.name + "/"
		UploadGS(args.vcf, vcf_gcs)
				# Copying the index file
		UploadGS(args.vcf + ".tbi", vcf_gcs)

	# Set up batches of files
	file_list = args.file_list
	sample_batches_paths= \
		GetVCFBatches(file_list, int(args.batch_size), int(args.batch_num), \
			gsprefix = bucket + "/" + "tr_imputation" +"/" + + "tr_imputation" +"/"+ str(args.batch_size))
	


	# Upload subset sample file
	if args.samples.startswith("gs://"):
		sample_file = args.samples
	else:
				# Copying the exclude sample file
		sample_file = output_bucket + "/" + args.name + "/"
		UploadGS(args.samples, sample_file)


	if args.subset_region and args.region is None:
		ERROR("Must specify --region for --subset-region")



	# Set up workflow JSON
	json_dict = {}
	json_dict["imputation.vcf"] = args.vcf
	json_dict["imputation.vcf_index"]=args.vcf+".tbi"
	json_dict["imputation.ref_panel"] = args.ref_panel
	json_dict["imputation.out_prefix"] = args.name
	json_dict["imputation.GOOGLE_PROJECT"] = project
	#json_dict["imputation.GCS_OAUTH_TOKEN"] = token
	json_dict["imputation.mem"] = args.mem
	json_dict["imputation.window_size"] = args.window
	json_dict["imputation.sample_file"] = args.samples 
	json_dict["imputation.region"] = args.region
	json_dict["imputation.subset_region"] = args.subset_region 
	json_dict["imputation.file_list"] = args.file_list
	json_dict["imputation.batch_size"] = args.batch_size
	json_dict["imputation.batch_num"] = args.batch_num

	


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
	RunWorkflow(json_file, json_options_file, args.cromwell, dryrun=args.dryrun)


if __name__ == "__main__":
	main()
