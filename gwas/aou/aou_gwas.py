#!/usr/bin/env python3

"""
Run GWAS on All of Us
"""

import argparse
import os
import re
from utils import MSG, ERROR

GWAS_METHODS = ["hail"]

def GetPTCovarPath(phenotype):
	return os.path.join(os.getenv('WORKSPACE_BUCKET'), \
		"phenotypes", "%s_phenocovar.csv"%phenotype)

def CheckRegion(region):
	if region is None: return True
	return re.match(r"\w+:\d+-\d+", query_string) is not None

def GetCovars(covar_str):
	if covar_str == "":
		return []
	else:
		return covar_str.split(",")

def GetOutPath(phenotype, method, region):
	outprefix = "%s_%s"%(phenotype, method)
	if region is not None:
		outprefix += "_%s"%(region.replace(":", "_").replace("-","_"))
	return outprefix += ".gwas.tab"

def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--phenotype", help="Phenotypes file path, or phenotype name", type=str, required=True)
    parser.add_argument("--method", help="GWAS method. Options: %s"",".join(GWAS_METHODS), type=str, default="hail")
    parser.add_argument("--region", help="chr:start-end to restrict to", type=str)
    parser.add_argument("--num-pcs", help="Number of PCs to use as covariates", type=int, default=10)
    parser.add_argument("--covars", help="Comma-separated list of phenotype-specific covariates. Default: age", type=str, default="age")
    parser.add_argument("--no-sex", help="Do not include sex as a covariate", action="store_true")
	args = parser.parse_args()

	# Set up paths
	if args.phenotype.endswith(".csv"):
		ptcovar_path = args.phenotype
	else:
		ptcovar_path = GetPTCovarPath(args.phenotype)

	# Check options
	if args.method not in GWAS_METHODS:
		ERROR("--method must be one of: %s"%",".join(GWAS_METHODS))
	if not CheckRegion(args.region):
		ERROR("Invalid region %s"%args.region)

	# TODO - deal with which cohort to do
	# e.g. all, EUR, AFR, etc.
	# pass sample list to gwas below

	# Set up GWAS method
	if args.method == "hail":
		runner = HailRunner(ptcovar_path, region=args.region, \
			num_pcs=args.num_pcs, ptcovars=GetCovars(args.covars), \
			no_sex=args.no_sex, \
			out_path=GetOutPath(args.phenotype, args.method, args.region))
	else:
		ERROR("GWAS method %s not implemented")

	# Run GWAS
	runner.RunGWAS()

	# Plot Manhattan and QQ - TODO

if __name__ == "__main__":
    main()