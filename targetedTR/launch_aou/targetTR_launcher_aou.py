#!/usr/bin/env python3
"""
Script to launch AOU targeted TR analysis

Desired usage:
./targetTR_launcher_aou.py \
  --region chr21:43776445-43776479 \
  --period 5 \
  --refcopies 7.0 \
  --name CSTB-mini \
  --batch-size 5 \
  --batch-num 3 \
  --file-list manifest.csv

"""

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile

def ParseRegion(str_region):
	"""
	Parse region into chrom, start, end

	Arguments
	---------
	str_region : str
	   Region in the format chr:start-end

	Returns
	-------
	chrom : str
	   Chromosome
	start : int
	   start coordinate
	end : int
	   end coordinate
	"""
	chrom = str_region.split(":")[0]
	start = int(str_region.split(":")[1].split("-")[0])
	end = int(str_region.split(":")[1].split("-")[1])
	return chrom, start, end

def GetFileBatches(file_list, batch_size, \
	batch_num=-1, gsprefix=None, action="both"):
	"""
	Given a file containing a list of file IDs for the 
	crams and their indices, divide into batches of 
	size batch_size. 

	Arguments
	---------
	file_list : str
	   File containing a list of IDs for crams and indices
	   This is the manifest.csv file from AoU
	batch_size : int
	   Number of cram files to include in each batch
	batch_num : int
	   Number of batches to create. If -1, make batches
	   out of everything. Primarily used for debugging
	   to make small test sets
	gsprefix : str
	   Prefix to upload to google cloud
	action : str
	   If set to "run-batches", we don't actually
	   need to create any files, just return file names

	Returns
	-------
	cram_batches_paths : list of str (GCP paths)
	   Each path links to a google cloud file listing the
	   crams to be processed
	cram_idx_batches_paths : list of str (GCP paths)
	   Each path links to a google cloud file listing the
	   indices of crams to be processed
	"""

	# Keep track of batches
	cram_batches = []
	cram_idx_batches = []

	# Crams/indices in current batch
	curr_batch_crams = []
	curr_batch_indices = []
	with open(file_list, "r") as f:
		for line in f:
			if line.startswith("person"): continue # header
			# If we reached the batch size, add batch to 
			# the list and start a new one
			if len(curr_batch_crams) == batch_size:
				cram_batches.append(curr_batch_crams)
				cram_idx_batches.append(curr_batch_indices)
				# Quit if we reached our desired number of batches
				# If -1 we just keep going until we run out of files
				if len(cram_batches) == batch_num: break 
				curr_batch_crams = []
				curr_batch_indices = []
			# Add this cram to the current batch
			cram_id, idx_id = line.strip().split(",")[1:3]
			curr_batch_crams.append(cram_id)
			curr_batch_indices.append(idx_id)
	assert(len(cram_batches) == len(cram_idx_batches))
	cram_batches_paths = []
	cram_idx_batches_paths = []
	for i in range(len(cram_batches)):
		cram_batch_fname = "batch-%s-crams.txt"%i
		cram_index_batch_fname = "batch-%s-crams-indices.txt"%i
		##### only actually generate batch files
		##### if action is either "create-batches"
		##### or "both"
		if action in ["both", "create-batches"]:
			with open(cram_batch_fname, "w") as f:
				for item in cram_batches[i]:
					f.write(item.strip()+"\n")
			with open(cram_index_batch_fname, "w") as f:
				for item in cram_idx_batches[i]:
					f.write(item.strip()+"\n")
			if gsprefix is not None:
				UploadGS(cram_batch_fname, gsprefix + "/" + cram_batch_fname)
				UploadGS(cram_index_batch_fname, gsprefix + "/" + cram_index_batch_fname)
		##########################
		cram_batches_paths.append(gsprefix + "/" + cram_batch_fname)
		cram_idx_batches_paths.append(gsprefix + "/" + cram_index_batch_fname)
	return cram_batches_paths, cram_idx_batches_paths

def RunWorkflow(json_file, json_options_file, wdl_dependencies_file="", dryrun=False):
	"""
	Run workflow on AoU

	Arguments
	---------
	json_file : str
	    JSON file path with input arguments
	json_options_file : str
	    JSON with additional options for cromshell
	wdl_dependencies_file : str (optional)
	    Zip file with other WDLs that are imported
	dryrun : bool
	    Just print the command, don't actually run cromshell
	"""
	cmd = "cromshell submit ../wdl/targetTR.wdl {json} -op {options}".format(json=json_file, options=json_options_file)
	if wdl_dependencies_file.strip() != "":
		cmd += " -d {otherwdl}".format(otherwdl=wdl_dependencies_file)
	if dryrun:
		sys.stderr.write("Run: %s\n"%cmd)
		return
	output = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).stdout.read()
	print(output.decode("utf-8"))

def UploadGS(tr_bedfile, tr_bedfile_gcs):
	cmd = "gsutil cp {src} {dest}".format(src=tr_bedfile, dest=tr_bedfile_gcs)
	output = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).stdout.read()
	print(output.decode("utf-8"))	

def WriteTRBed(region, period, refcopies, name, filename):
	chrom, start, end = ParseRegion(region)
	with open(filename, "w") as f:
		f.write("\t".join([chrom, str(start), str(end), str(period), str(refcopies), name])+"\n")

def ZipWDL(wdl_dependencies_file):
	"""
	Put all WDL dependencies into a zip file

	Arguments
	---------
	wdl_dependencies_fie : str
	    Zip file to put other wdls in
	"""
	files = ["hipstr_multi.wdl", "merge_hipstr.wdl", "dumpstr.wdl"]
	dirname = tempfile.mkdtemp()
	for f in files:
		shutil.copyfile("../wdl/%s"%f, dirname+"/"+f)
	shutil.make_archive(wdl_dependencies_file.strip(".zip"), "zip", root_dir=dirname)

def main():
	parser = argparse.ArgumentParser(__doc__)
	parser.add_argument("--region", help="chrom:start-end of TR region", required=True, type=str)
	parser.add_argument("--period", help="Repeat unit length (bp) of TR", required=True, type=int)
	parser.add_argument("--refcopies", help="Ref num. of copies of the TR", required=True, type=float)
	parser.add_argument("--name", help="Name of the TR job", required=True, type=str)
	parser.add_argument("--batch-size", help="HipSTR batch size", required=False, type=int, default=1000)
	parser.add_argument("--batch-num", help="Number of batches. Default: -1 (all)", required=False, default=-1)
	parser.add_argument("--file-list", help="List of crams and indices to process (manifest.csv)"
		"Format: person_id,cram_uri,cram_index_uri", type=str, required=True)
	parser.add_argument("--genome-id", help="File id of ref genome", type=str, default="gs://genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.fasta")
	parser.add_argument("--genome-idx-id", help="File id of ref genome index", type=str, default="gs://genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.fasta.fai")
	parser.add_argument("--action", help="Options: create-batches, run-batches, both", type=str, default="both")
	parser.add_argument("--dryrun", help="Don't actually run the workflow. Just set up", action="store_true")
	args = parser.parse_args()

	# Check if action is valid
	if args.action not in ["create-batches", "run-batches", "both"]:
		sys.stderr.write("Invalid action: %s\n"%args.action)
		sys.exit(1)

	# Get token
	token_fetch_command = subprocess.run(['gcloud', 'auth', 'print-access-token', env['user_email']], \
		capture_output=True, check=True, encoding='utf-8')
	token = str.strip(token_fetch_command.stdout)

	# Set up output bucket
	bucket = os.getenv("WORKSPACE_BUCKET")
	project = os.getenv("GOOGLE_PROJECT")
	output_bucket = bucket + "/" + args.name

	# Set up batches of files
	cram_batches_paths, cram_idx_batches_paths = \
		GetFileBatches(args.file_list, int(args.batch_size), int(args.batch_num), \
			gsprefix = output_bucket + "/" + args.name +"/" + str(args.batch_size), action=args.action)
	if args.action == "create-batches":
		# We're done! quit before running jobs
		sys.exit(1)

	# Generate TR Bed file
	tr_bedfile = args.name+".bed"
	tr_bedfile_gcs = output_bucket + "/" + args.name + "/" + tr_bedfile
	WriteTRBed(args.region, args.period, args.refcopies, args.name, tr_bedfile)
	UploadGS(tr_bedfile, tr_bedfile_gcs)

	# Set up workflow JSON
	json_dict = {}
	json_dict["targetTR.genome"] = args.genome_id
	json_dict["targetTR.genome_index"] = args.genome_idx_id
	json_dict["targetTR.tr_bed"] = tr_bedfile_gcs
	json_dict["targetTR.str_name"] = args.name
	json_dict["targetTR.using_aou"] = True
	json_dict["targetTR.ukb_names"] = False
	json_dict["targetTR.GOOGLE_PROJECT"] = project
	json_dict["targetTR.GCS_OAUTH_TOKEN"] = token
	json_dict["targetTR.cram_file_batches_str"] = cram_batches_paths

	# Convert to json and save as a file
	json_file = args.name+".aou.json"
	with open(json_file, "w") as f:
		json.dump(json_dict, f, indent=4)

	# Set up json options
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
