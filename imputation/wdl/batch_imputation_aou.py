#!/usr/bin/env python3
"""
Script to launch AOU imputation use new ref panel 

example code to impute 10 samples at CBL region 

./batch_imputation_aou.py \
--name batch_test 
#--samples-batch $WORKSPACE_BUCKETacaf_batches/manifest/chr11_acaf_manifest.txt \
--ref-panel $WORKSPACE_BUCKET/tr_imputation/tr_imputation/ref/chr11_final_SNP_merged_additional_TRs.bref3 \
--mem 40 \
--batch-num 2
--map $WORKSPACE_BUCKET/tr_imputation/tr_imputation/genetic_map/beagle_chr11_b38.map \
--disk 50 \
"""

# need to change GetFileBatches
import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
import csv
from utils import MSG, ERROR

def GetBatchVCFFiles(vcfdir, max_batches):
	cmd = "gsutil ls %s/*.vcf.gz"%(vcfdir)
	output = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).stdout.read().decode("utf-8")
	batch_files = [item.strip() for item in output.split()]
	if max_batches > -1 and max_batches <= len(batch_files):
		batch_files = batch_files[0:max_batches]
	return batch_files

def RunWorkflow(json_file, json_options_file, wdl_dependencies_file, dryrun=False):
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

	cmd = "cromshell submit ../wdl/batch_imputation.wdl {json} -op {options}".format(json=json_file, options=json_options_file)
	
	if wdl_dependencies_file.strip() != "":
		cmd += " -d {otherwdl}".format(otherwdl=wdl_dependencies_file)
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
	
##def UploadGS(local_path, gcp_path):
#	"""
#	Upload a local file to GCP
#
#	Arguments
#	---------
#	local_path : str
#	   Local path
#	gcp_path : str
#	   GCP path to upload to
#	"""
#	cmd = "gsutil cp {src} {dest}".format(src=local_path, dest=gcp_path)
#	output = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).stdout.read()
#	print(output.decode("utf-8"))	

def ZipWDL(wdl_dependencies_file):
	"""
	Put all WDL dependencies into a zip file

	Arguments
	---------
	wdl_dependencies_fie : str
	    Zip file to put other wdls in
	"""
	files = ["beagle.wdl","processTR.wdl","merge_TR_batch.wdl"]
	dirname = tempfile.mkdtemp()
	for f in files:
		shutil.copyfile("../wdl/%s"%f, dirname+"/"+f)
	shutil.make_archive(os.path.splitext(wdl_dependencies_file)[0], "zip", root_dir=dirname)


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
	parser.add_argument("--chrom", help="Which chromosome to process", type=str, required=True)
	parser.add_argument("--mem", help="Specify run memory for beagle ", type=int, required=False, default=50)
	parser.add_argument("--merge-mem", help="Specify run memory for bcftools merge ", type=int, required=False, default=50)
	parser.add_argument("--disk", help="Specify disk memory ", type=int, required=False, default=25)
	parser.add_argument("--window", help="Specify window size for imputation ", type=int, required=False, default=20)
	#parser.add_argument("--sample-batches", help="List of batches samples to process ", type=str, required=False, \
	#				 default=f"{bucket}/tr_imputation/tr_imputation/sample/sample_manifest.txt")
	parser.add_argument("--sample-batches", help="List of batches samples to process ", type=str, required=False, \
					 default=f"{bucket}/acaf_batches/manifest/chr11_acaf_manifest.txt")
	parser.add_argument("--region", help="Name of chrom position  chr:xxx-xxx", type=str,required=False)
	parser.add_argument("--beagle-region", help="Apply chrom for beagle", action="store_true",required=False)
	parser.add_argument("--dryrun", help="Don't actually run the workflow. Just set up", action="store_true")
	parser.add_argument("--batch-num", help="Number of batches. Default: -1 (all)",type=int, required=False, default=None)
	parser.add_argument("--overlap", help="Specify overlap size for imputation ", type=int, required=False, default=2)
	parser.add_argument("--map", help="Specify genetic map for imputation ", type=str, required=True)
				

	args = parser.parse_args()

	batch_vcf_files = GetBatchVCFFiles(args.vcfdir, args.max_batches)
	output_bucket = bucket + "/" + args.name

	#Set up sample file list
	if args.sample_batches.startswith("gs://"):
		DownloadGS(args.sample_batches)
		sample_list = os.path.basename(args.sample_batches)
	else: sample_list = args.sample_batches

	if args.beagle_region and args.region is None:
		ERROR("Must specify --region for --beagle-region")


	# Set up workflow JSON
	json_dict = {}
	json_dict["batch_imputation.ref_panel"] = bucket + "tr_imputation/tr_imputation/ref_panel/ensembletr_refpanel_v3_chr%s.bref3"%args.chrom
	json_dict["batch_imputation.ref_vcf"] = bucket + "tr_imputation/tr_imputation/ref_panel/ensembletr_refpanel_v3_chr%s.vcf.gz"%args.chrom
	json_dict["batch_imputation.ref_index"] = bucket + "tr_imputation/tr_imputation/ref_panel/ensembletr_refpanel_v3_chr%s.vcf.gz.tbi"%args.chrom
	json_dict["batch_imputation.out_prefix"] = args.name
	json_dict["batch_imputation.GOOGLE_PROJECT"] = project
	json_dict["batch_imputation.mem"] = args.mem
	json_dict["batch_imputation.merge_mem"] = args.merge_mem
	json_dict["batch_imputation.disk"] = args.disk
	json_dict["batch_imputation.window_size"] = args.window
	json_dict["batch_imputation.sample_batches"] = batch_vcf_files 
	json_dict["batch_imputation.region"] = args.region 
	json_dict["batch_imputation.beagle_region"] = args.beagle_region
	json_dict["batch_imputation.overlap"] = args.overlap
	json_dict["batch_imputation.map"] = args.map




	# Convert to json and save as a file
	json_file = args.name+".aou.json"
	with open(json_file, "w") as f:
		json.dump(json_dict, f, indent=4)


	# Set up json optionsß
	json_options_dict = {}
	json_options_file = args.name+".options.aou.json"
	with open(json_options_file, "w") as f:
		json.dump(json_options_dict, f, indent=4)

	# Zip all the WDL depencies
	wdl_dependencies_file = args.name + "-wdl.zip"
	ZipWDL(wdl_dependencies_file)


	# Run workflow on AoU using cromwell
	RunWorkflow(json_file, json_options_file, wdl_dependencies_file, dryrun=args.dryrun)


if __name__ == "__main__":
	main()
