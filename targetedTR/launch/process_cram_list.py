#!/usr/bin/env python3
"""
./process_cram_list.py ukb_cram_files_long.txt > ukb_cram_and_index_files.txt

Output file of:
<file ID cram> <file ID index>
"""

import sys

fdict = {} # sample -> {cram, idx}

with open(sys.argv[1], "r") as f:
	for line in f:
		if not line.startswith("closed"): continue
		items = line.strip().split()
		fname = items[5]
		fid = items[6].replace("(","").replace(")","")
		sample = fname.split(".")[0]
		if sample not in fdict:
			fdict[sample] = {}
		if fname.endswith(".crai"):
			fdict[sample]["idx"] = fid
		else:
			fdict[sample]["cram"] = fid

for sample in fdict.keys():
	cramid = fdict[sample].get("cram", None)
	idxid = fdict[sample].get("idx", None)
	if cramid is None or idxid is None:
		sys.stderr.write("Couldn't find both files for %s"%sample)
	else:
		sys.stdout.write(" ".join([cramid, idxid])+"\n")
