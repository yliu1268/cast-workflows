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
import subprocess

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

def GetFileBatches(file_list, batch_size, batch_num=-1):
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

	Returns
	-------
	cram_batches : list of lists
	   Each internal list is a list of file IDs for crams
	cram_idx_batches : list of lists
	   Each internal list is a list of file IDs for cram indices
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
	return cram_batches, cram_idx_batches

def RunWorkflow(json_file, json_options_file):
	"""
	Run workflow on AoU

	Arguments
	---------
	json_file : str
	    JSON file path with input arguments
	name : str
	    Used to determine where to store output
	"""
	conf_file = "/home/jupyter/cromwell.conf"
	cmd = "java -jar -Dconfig.file={conf} cromwell-77.jar run ../wdl/workflows/targetTR.wdl -i {json} -o {options}".format(conf=conf_file, json=json_file, options=json_options_file)
	output = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).stdout.read()
	print(output.decode("utf-8"))

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
	args = parser.parse_args()

	# Set up workflow JSON
	json_dict = {}
	json_dict["targetTR.genome"] = args.genome_id
	json_dict["targetTR.genome_index"] = args.genome_idx_id
	json_dict["targetTR.str_name"] = args.name
	json_dict["targetTR.num_copies"] = args.refcopies
	json_dict["targetTR.motif_len"] = args.period
	
	# Parser region and set in json
	chrom, start, end = ParseRegion(args.region)
	json_dict["targetTR.chrom"] = chrom
	json_dict["targetTR.str_start"] = start
	json_dict["targetTR.str_end"] = end

	# Set up batches of files
	cram_batches, cram_idx_batches = GetFileBatches(args.file_list, int(args.batch_size), int(args.batch_num))
	json_dict["targetTR.cram_file_batches"] = cram_batches
	json_dict["targetTR.cram_index_batches"] = cram_idx_batches

	# Convert to json and save as a file
	json_file = args.name+".aou.json"
	with open(json_file, "w") as f:
		json.dump(json_dict, f, indent=4)

	# Set up json options
	bucket = os.getenv("WORKSPACE_BUCKET")
	project = os.getenv("GOOGLE_PROJECT")
	output_bucket = bucket + "/" + args.name
	json_options_dict = {"jes_gcs_root": output_bucket}
	json_options_file = args.name+".options.aou.json"
	with open(json_options_file, "w") as f:
		json.dump(json_options_dict, f, indent=4)

	# Run workflow on AoU using cromwell
	RunWorkflow(json_file, json_options_file)

if __name__ == "__main__":
	main()
