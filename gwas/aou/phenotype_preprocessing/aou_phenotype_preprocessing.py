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
import os
import pandas as pd
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

	# Check if we have info for that phenotype
	if args.phenotype not in aou_queries.pt_queries.keys():
		ERROR("Could not find sql query for %s"%args.phenotype)

	# Set up dataframes
	demog = pd.read_gbq(
    	aou_queries.demographics_sql,
    	dialect="standard",
    	use_bqstorage_api=("BIGQUERY_STORAGE_API_ENABLED" in os.environ),
    	progress_bar_type="tqdm_notebook")
	ptdata = pd.read_gbq(
    	aou_queries.pt_queries[args.phenotype],
    	dialect="standard",
    	use_bqstorage_api=("BIGQUERY_STORAGE_API_ENABLED" in os.environ),
    	progress_bar_type="tqdm_notebook")
	data = pd.merge(ptdata, demog, on="person_id", how="inner")
	print(data.head())

if __name__ == "__main__":
	main()