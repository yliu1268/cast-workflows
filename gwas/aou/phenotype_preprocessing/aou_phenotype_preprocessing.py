#!/usr/bin/env python3

"""
Preprocess AoU phenotypes.
Output cleaned phenotype and covariates files
for the specified phenotype

Requires these environment variables to be set:
WORKSPACE_CDR
"""

import argparse
import aou_queries
import sys

def MSG(msg_str):
	sys.stderr.write("[aou_phenotype_preprocessing]: %s\n"%msg_str.strip())

def ERROR(msg_str):
	MSG("ERROR: " + msg_str)
	sys.exit(1)

def main():
	parser = argparse.ArgumentParser(__doc__)
	parser.add_argument("--phenotype", help="Phenotype ID", type=str, required=True)
	args = parser.parse_args()
	MSG("Processing %s"%args.phenotype)

	# Pull out sql queries
	if args.phenotype not in aou_queries.pt_queries.keys():
		ERROR("Could not find sql query for %s"%args.phenotype)
	demog_sql = aou_queries.demographics_sql
	pt_sql = aou_queries.pt_queries[args.phenotype]



if __name__ == "__main__":
	main()