#!/usr/bin/env python3
"""
Script to launch UKB targeted TR analysis

Desired usage:
./targetTR_launcher_ukb.py \
  --region chr21:43776445-43776479 \
  --period 5 \
  --refcopies 7.0 \
  --name CSTB-mini \
  --batch-size 100 \
  --batch-num 3 \
  --workflow-id workflow-GP8jzvQJv7B3K97ypgvBBqxq \
  --file-list ukb_cram_and_index_files.txt

"""

import argparse
import json

def ParseRegion(str_region):
	chrom = str_region.split(":")[0]
	start = int(str_region.split(":")[1].split("-")[0])
	end = int(str_region.split(":")[1].split("-")[1])
	return chrom, start, end

def GetFileBatches(file_list, batch_size, batch_num):
	cram_batches = []
	cram_idx_batches = []

	curr_batch_crams = []
	curr_batch_indices = []
	with open(file_list, "r") as f:
		for line in f:
			# If we reached the batch size, add to the list and start a new one
			if len(curr_batch_crams) == batch_size:
				cram_batches.append(curr_batch_crams)
				cram_idx_batches.append(curr_batch_indices)
				# Quit if we reached our desired number of batches
				if len(cram_batches) == batch_num: break 
				curr_batch_crams = []
				curr_batch_indices = []
			# Add this cram to the current batch
			cram_id, idx_id = line.strip().split()
			curr_batch_crams.append({"$dnanexus_link": cram_id})
			curr_batch_indices.append({"$dnanexus_link": idx_id})

	return cram_batches, cram_idx_batches

def main():
	parser = argparse.ArgumentParser(__doc__)
	parser.add_argument("--region", help="chrom:start-end of TR region", required=True, type=str)
	parser.add_argument("--period", help="Repeat unit length (bp) of TR", required=True, type=int)
	parser.add_argument("--refcopies", help="Ref num. of copies of the TR", required=True, type=float)
	parser.add_argument("--name", help="Name of the TR job", required=True, type=str)
	parser.add_argument("--batch-size", help="HipSTR batch size", required=False, type=int, default=1000)
	parser.add_argument("--batch-num", help="Number of batches. Default: all", required=False, default="all")
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
	json_dict["stage-common.num_copies"] = args.refcopies
	json_dict["stage-common.motif_len"] = args.period
	
	# Parser region and set in json
	chrom, start, end = ParseRegion(args.region)
	json_dict["stage-common.chrom"] = chrom
	json_dict["stage-common.str_start"] = start
	json_dict["stage-common.str_end"] = end

	# Set up batches of files
	cram_batches, cram_idx_batches = GetFileBatches(args.file_list, int(args.batch_size), int(args.batch_num))
	json_dict["stage-common.cram_file_batches"] = cram_batches
	json_dict["stage-common.cram_index_batches"] = cram_idx_batches

	# Convert to json
	json_string = json.dumps(json_dict)

	# Run workflow on dna nexus - TODO
	print(json_string)

if __name__ == "__main__":
	main()
