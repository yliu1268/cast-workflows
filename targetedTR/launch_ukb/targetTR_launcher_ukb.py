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
  --workflow-id workflow-GPfbXV8Jv7B27kpf6Y50QyQ9  \
  --file-list ukb_cram_and_index_files.txt
"""

import argparse
import dxpy
import json
import sys

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
			curr_batch_crams.append(dxpy.dxlink(cram_id))
			curr_batch_indices.append(dxpy.dxlink(idx_id))
	assert(len(cram_batches) == len(cram_idx_batches))
	return cram_batches, cram_idx_batches

def RunWorkflow(json_dict, workflow_id, name, depends=[]):
	"""
	Run DNA Nexus workflow

	Arguments
	---------
	json_dict : dict
		Input arguments for the workflow
	workflow_id : str
		ID of the DNA Nexus worfklow
	name : str
		Used to determine where to store output

	Returns
	-------
	analysis : DXAnalysis object
	"""
	workflow = dxpy.dxworkflow.DXWorkflow(dxid=workflow_id)
	analysis = workflow.run(json_dict, \
		depends_on=[item.get_id() for item in depends], \
		folder="/TargetedSTR/results/{name}".format(name=name))
	sys.stderr.write("Started analysis %s (%s)\n"%(analysis.get_id(), name))
	return analysis

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
	file : dxlink
	   {"$dnanexus_link": file_id}
	"""
	folder = "/TargetedSTR/results/{name}".format(name=name)
	dxfile = dxpy.upload_local_file(fname, folder=folder, parents=True)
	return dxpy.dxlink(dxfile)

def WriteTRBed(region, period, refcopies, name, filename):
	chrom, start, end = ParseRegion(region)
	with open(filename, "w") as f:
		f.write("\t".join([chrom, str(start), str(end), str(period), str(refcopies), name])+"\n")

def GetJBOR(analysis, filename):
	return {
		"$dnanexus_link": {
			"analysis": analysis.get_id(),
			"field": filename
		}
	}

def main():
	parser = argparse.ArgumentParser(__doc__)
	parser.add_argument("--str-ref", help="str reference bed file, containing the tab separated columns chrom, start, end, period, refcopies and name. Must either specify --str-ref (for many STRs) or all of --region, --period and --refcopies (for a single STR)", required=False, type=str)
	parser.add_argument("--region", help="chrom:start-end of TR region, start and end are bed file style (they are zero-indexed and start is inclusive but end is exclusive)", required=False, type=str)
	parser.add_argument("--period", help="Repeat unit length (bp) of TR", required=False, type=int)
	parser.add_argument("--refcopies", help="Ref num. of copies of the TR", required=False, type=float)
	parser.add_argument("--name", help="Name of the TR job", required=True, type=str)
	parser.add_argument("--batch-size", help="HipSTR batch size", required=False, type=int, default=500)
	parser.add_argument("--batch-num", help="Number of batches. Default: -1 (all)", required=False, default=-1)
	parser.add_argument("--workflow-id", help="DNA Nexus workflow ID", required=False, default="workflow-GPfbXV8Jv7B27kpf6Y50QyQ9 ")
	parser.add_argument("--file-list", help="List of crams and indices to process"
		"Format of each line: cram-file-id cram-index-id", type=str, required=True)
	parser.add_argument("--genome-id", help="File id of ref genome", type=str, default="file-GGJ1z28JbVqbpqB93YbPqbzz")
	parser.add_argument("--genome-idx-id", help="File id of ref genome index", type=str, default="file-GGJ94JQJv7BGFYq8BGp62xPV")
	# Options for multi-batches
	parser.add_argument("--merge-workflow-id", help="DNA Nexus workflow ID for merging", required=False, default="workflow-GPfv050Jv7B5ZQqpZb0jG5k6")
	parser.add_argument("--max-batches-per-workflow", help="Maximum number of batches to launch at once. Default: -1 (all)", required=False, default=-1, type=int)
	parser.add_argument("--concurrent", help="Launch all batches at once", action="store_true")
	args = parser.parse_args()

	assert bool(args.region) == bool(args.period) == bool(args.refcopies)
	assert bool(args.region) != bool(args.str_ref)

	# Set up workflow JSON
	json_dict = {}
	json_dict["stage-common.genome"] = {}
	json_dict["stage-common.genome_index"] = {}
	json_dict["stage-common.genome"] = dxpy.dxlink(args.genome_id)
	json_dict["stage-common.genome_index"] = dxpy.dxlink(args.genome_idx_id)
	json_dict["stage-common.str_name"] = args.name
	json_dict["stage-common.ukb_names"] = True
	
	# Make bed file
	sys.stderr.write("Preparing bed file...\n")
	if args.region:
		tr_bedfile = args.name+".bed"
		WriteTRBed(args.region, args.period, args.refcopies, args.name, tr_bedfile)
	else:
		tr_bedfile = args.str_ref
	json_dict["stage-common.tr_bed"] = UploadDNANexus(tr_bedfile, args.name)

	# Set up batches of files
	sys.stderr.write("Setting up batches...\n")
	cram_batches, cram_idx_batches = GetFileBatches(args.file_list, int(args.batch_size), int(args.batch_num))

	# Run batches
	final_vcf = None
	final_vcf_idx = None
	if args.max_batches_per_workflow == -1:
		# Run all at once
		json_dict["stage-common.cram_file_batches"] = cram_batches
		json_dict["stage-common.cram_index_batches"] = cram_idx_batches
		analysis = RunWorkflow(json_dict, args.workflow_id, args.name)
	else:
		# Run in chunks across multiple workflows
		depends = []
		batch_num = 0
		curr_idx = 0
		curr_cram_batches = []
		curr_idx_batches = []
		while curr_idx < len(cram_batches):
			if len(curr_cram_batches) == args.max_batches_per_workflow:
				batch_name = args.name + "-CHUNK%s"%batch_num
				batch_dict = json_dict.copy()
				batch_dict["stage-common.str_name"] = batch_name
				batch_dict["stage-common.cram_file_batches"] = curr_cram_batches
				batch_dict["stage-common.cram_index_batches"] = curr_idx_batches
				if args.concurrent:
					use_dep = []
				else: use_dep = depends
				analysis = RunWorkflow(batch_dict, args.workflow_id, \
					args.name + "/" + args.name+"_%s"%batch_num, depends=use_dep)
				depends.append(analysis)
				batch_num += 1
				curr_cram_batches = []
				curr_idx_batches = []
			curr_cram_batches.append(cram_batches[curr_idx])
			curr_idx_batches.append(cram_idx_batches[curr_idx])
			curr_idx += 1
		# Run last batch if any are left
		if len(curr_cram_batches) > 0:
			batch_name = args.name + "-CHUNK%s"%batch_num
			batch_dict = json_dict.copy()
			batch_dict["stage-common.str_name"] = batch_name
			batch_dict["stage-common.cram_file_batches"] = curr_cram_batches
			batch_dict["stage-common.cram_index_batches"] = curr_idx_batches
			if args.concurrent:
				use_dep = []
			else: use_dep = depends
			analysis = RunWorkflow(batch_dict, args.workflow_id, \
				args.name + "/" + args.name+"_%s"%batch_num, depends=use_dep)
			depends.append(analysis)
		# Run a final job to merge all the meta-batches
		merge_vcfs = []
		merge_vcfs_idx = []
		for analysis in depends:
			#analysis.wait_on_done()
			#vcf = analysis.describe()["output"]["stage-outputs.finalvcf"]
			#vcf_idx = analysis.describe()["output"]["stage-outputs.finalvcf_index"]
			vcf = GetJBOR(analysis, "stage-outputs.finalvcf")
			vcf_idx = GetJBOR(analysis, "stage-outputs.finalvcf_index")
			merge_vcfs.append(vcf)
			merge_vcfs_idx.append(vcf_idx)
		sys.stderr.write("Setting up merge from meta-batches...\n")
		merge_dict = {}
		merge_dict["stage-common.out_name"] = args.name
		merge_dict["stage-common.vcf_files"] = merge_vcfs
		merge_dict["stage-common.vcf_idxs"] = merge_vcfs_idx
		analysis = RunWorkflow(merge_dict, args.merge_workflow_id, args.name, depends=depends)

if __name__ == "__main__":
	main()
