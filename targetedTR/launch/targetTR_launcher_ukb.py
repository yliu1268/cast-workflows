#!/usr/bin/env python3
"""
Script to launch UKB targeted TR analysis

Desired usage:
./targetTR_launcher_ukb.py \
  --region chr21:43776445-43776479 \
  --period 5 \
  --refcopies 7.0 \
  --name CSTB-mini \
  --batch-size 5 \
  --batch-num 3 \
  --workflow-id workflow-GP8jzvQJv7B3K97ypgvBBqxq \
  --file-list ukb_cram_and_index_files.txt

"""

import argparse
import json
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
	   File must have one row per CRAM/index, and the file
	   IDs must be separated by whitespace
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
	   To play nicely with dxCompiler json, each item
	   in the internal list has format {"$dnanexus_link": cram_id}
	cram_idx_batches : list of lists
	   Each internal list is a list of file IDs for cram indices
	   To play nicely with dxCompiler json, each item
	   in the internal list has format {"$dnanexus_link": cram_id}
	"""

	# Keep track of batches
	cram_batches = []
	cram_idx_batches = []

	# Crams/indices in current batch
	curr_batch_crams = []
	curr_batch_indices = []
	with open(file_list, "r") as f:
		for line in f:
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
			cram_id, idx_id = line.strip().split()
			curr_batch_crams.append({"$dnanexus_link": cram_id})
			curr_batch_indices.append({"$dnanexus_link": idx_id})
	assert(len(cram_batches) == len(cram_idx_batches))
	return cram_batches, cram_idx_batches

def RunWorkflow(json_file, workflow_id, name):
	"""
	Run DNA Nexus workflow

	Arguments
	---------
	json_file : str
	    JSON file path with input arguments
	workflow_id : str
	    ID of the DNA Nexus worfklow
	name : str
	    Used to determine where to store output
	"""
	cmd = 'dx run {workflow} -y -f {json} --destination "TargetedSTR/results/{name}"'.format(workflow=workflow_id, json=json_file, name=name)
	output = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).stdout.read()
	print(output.decode("utf-8"))

def UploadDNANexus(fname, name):
	"""
	Upload a file to DNA Nexus and return
	its file ID

	Arguments
	---------
	fname : str
	   Path to file to upload
	name : str
	   Name of subdirectory to store in

	Returns
	-------
	file_id : str
	   ID of the file on DNA Nexus
	"""
	cmd = 'dx upload {fname} --brief ‑‑destination "TargetedSTR/results/{name}"'.format(fname, name=name)
	output = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).stdout.read()
	return output

def WriteTRBed(region, period, refcopies, name, filename):
	chrom, start, end = ParseRegion(region)
	with open(filename, "w") as f:
		f.write("\t".join([chrom, str(start), str(end), str(period), str(refcopies), name])+"\n")

def main():
	parser = argparse.ArgumentParser(__doc__)
	parser.add_argument("--region", help="chrom:start-end of TR region", required=True, type=str)
	parser.add_argument("--period", help="Repeat unit length (bp) of TR", required=True, type=int)
	parser.add_argument("--refcopies", help="Ref num. of copies of the TR", required=True, type=float)
	parser.add_argument("--name", help="Name of the TR job", required=True, type=str)
	parser.add_argument("--batch-size", help="HipSTR batch size", required=False, type=int, default=1000)
	parser.add_argument("--batch-num", help="Number of batches. Default: -1 (all)", required=False, default=-1)
	parser.add_argument("--workflow-id", help="DNA Nexus workflow ID", required=False, default="workflow-GP8jzvQJv7B3K97ypgvBBqxq")
	parser.add_argument("--file-list", help="List of crams and indices to process"
		"Format of each line: cram-file-id cram-index-id", type=str, required=True)
	parser.add_argument("--genome-id", help="File id of ref genome", type=str, default="file-GGJ1z28JbVqbpqB93YbPqbzz")
	parser.add_argument("--genome-idx-id", help="File id of ref genome index", type=str, default="file-GGJ94JQJv7BGFYq8BGp62xPV")
	args = parser.parse_args()

	# Set up workflow JSON
	json_dict = {}
	json_dict["stage-common.genome"] = {}
	json_dict["stage-common.genome_index"] = {}
	json_dict["stage-common.genome"]["$dnanexus_link"] = args.genome_id
	json_dict["stage-common.genome_index"]["$dnanexus_link"] = args.genome_idx_id
	json_dict["stage-common.str_name"] = args.name
	
	# Make bed file
	tr_bedfile = args.name+".bed"
	WriteTRBed(args.region, args.period, args.refcopies, args.name, tr_bedfile)
	json_dict["stage-common.tr_bed"] = UploadDNANexus(tr_bedfile, args.name)

	# Set up batches of files
	cram_batches, cram_idx_batches = GetFileBatches(args.file_list, int(args.batch_size), int(args.batch_num))
	json_dict["stage-common.cram_file_batches"] = cram_batches
	json_dict["stage-common.cram_index_batches"] = cram_idx_batches

	# Convert to json and save as a file
	json_file = args.name+".json"
	with open(json_file, "w") as f:
		json.dump(json_dict, f, indent=4)

	# Run workflow on dna nexus
	RunWorkflow(json_file, args.workflow_id, args.name)

if __name__ == "__main__":
	main()
