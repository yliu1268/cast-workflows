#!/usr/bin/env python3
"""
Usage:

./concanate_batches_v2.py <cromshell_job_id> <chrom>
"""

DEBUG = True

import json
import os
import subprocess
import sys

cromshell_job_id = sys.argv[1]
chrom = sys.argv[2]

def SortByCoordinate(vcf_files):
	""" Sort files by coordinates in filename """
	# Filename syntax: bucket/batch99-chr11-50000000-60000000.vcf.gz
	files = [] # (startcoord, file)
	for f in vcf_files:
		coord = int(os.path.basename(f).split("-")[2])
		files.append((coord, f))
	files = sorted(files, key = lambda x: x[0])
	return [item[1] for item in files]

# Load jobdata
jobdata_fetch_command = subprocess.run(['cromshell', '-mc', 'list-outputs', '-j', cromshell_job_id], \
			capture_output=True, check=True, encoding='utf-8')
jobdata = json.loads(str.strip(jobdata_fetch_command.stdout))

# Gather files for each batch
batch_files = {} # batchname -> [list of files]
for region in jobdata["subset_vcf.vcf_output_array"]:
	for f in region:
		batchname = os.path.basename(f).split("-")[0]
		if batchname not in batch_files.keys():
			batch_files[batchname] = []
		batch_files[batchname].append(f)
batch_names = list(batch_files.keys())

# Set up workflow json
json_dict = {}
json_dict["concatenate_batch_vcfs.outprefix"] = chrom
json_dict["concatenate_batch_vcfs.GOOGLE_PROJECT"] = os.environ.get("GOOGLE_PROJECT", "")
json_dict["concatenate_batch_vcfs.batch_names"] = batch_names
json_dict["concatenate_batch_vcfs.batch_vcf_files"] = [batch_files[bn] for bn in batch_names]
json_file = chrom + "concat.aou.json"
print(json_dict)
with open(json_file, "w") as f:
	json.dump(json_dict, f, indent=4)

# Set up options json
options_file = chrom + "concat.aou.options.json"
json_options_dict = {}
json_options_dict["final_workflow_outputs_dir"] = os.environ["WORKSPACE_BUCKET"] + "/acaf_batches/" + chrom
with open(options_file, "w") as f:
	json.dump(options_dict, f, indent=4)

# Run workflow
cmd = "cromshell submit {wdl} {json} -op {options}".format(wdl="wdl/concatenate_batch_vcfs.wdl", json=json_file, options=options_file)
print(cmd)
#output = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).stdout.read()
#print(output.decode("utf-8"))	


