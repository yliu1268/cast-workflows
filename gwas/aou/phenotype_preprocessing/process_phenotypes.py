#!/usr/bin/env python3
"""
Run phenotype preprocessing on a list of phenotypes

Usage:
./process_phenotypes.py <phenotype_list.csv>

Phenotype list has columns:
* `phenotype`
* `concept_id`
* `units`
* `min`
* `max`
* `drugcovars`
"""

import hashlib
import os
import pandas as pd
import subprocess
import sys

try:
    ptfile = sys.argv[1]
except:
    sys.stderr.write(__doc__)
    sys.exit(1)

def md5(fname):
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()

def RunCmd(cmd):
    res = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    res.wait()
    retcode = res.returncode
    if retcode != 0:
    	sys.stderr.write("Command failed: %s"%cmd)
    	sys.exit(1)

def UploadGCS(outfile, gcsfile):
	cmd = "gsutil cp {outfile} {gcsfile}".format(outfile=outfile, gcsfile=gcsfile)
	RunCmd(cmd)

manifest = open(ptfile.replace(".csv",".manifest.csv"), "w")
ptdata = pd.read_csv(ptfile)
for index, row in ptdata.iterrows():
    cmd = """./aou_phenotype_preprocessing.py \
           --phenotype {phenotype} \
           --concept-id {concept} \
           --units {units} \
           --range {minval},{maxval}""".format(phenotype=row["phenotype"], \
           	concept=row["concept_id"], units=row["units"], \
           	minval=row["min"], maxval=row["max"])
    if str(row["drugcovars"]) != "nan":
        cmd += " --drugexposure-covariate-concept-ids {drug}".format(drug=row["drugcovars"])
    outfile = row["phenotype"]+"_phenocovar.csv"
    gcsfile = os.path.join(os.environ["WORKSPACE_BUCKET"], "phenotypes", outfile)
    RunCmd(cmd)
    UploadGCS(outfile, gcsfile)
    manifest.write(",".join([row["phenotype"], row["concept_id"], \
    	row["units"], row["min"], row["max"], row["drugcovars"], \
    	gcsfile, md5(outfile)])+"\n")

manifest.close()