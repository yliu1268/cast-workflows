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
import numpy as np
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
    parser.add_argument("--samples", help="List of sample IDs to keep", type=str)
    args = parser.parse_args()
    MSG("Processing %s"%args.phenotype)

    # Check if we have info for that phenotype
    if args.phenotype not in aou_queries.AVAILABLE_PHENOTYPES:
        ERROR("Could not find concept_id for %s. Please add it to aou_queries.py"%args.phenotype)

    # Set up dataframes
    demog = pd.read_gbq(
        aou_queries.demographics_sql,
        dialect="standard",
        use_bqstorage_api=("BIGQUERY_STORAGE_API_ENABLED" in os.environ),
        progress_bar_type="tqdm_notebook")
    ptdata = pd.read_gbq(
        aou_queries.ConstructTraitSQL(args.phenotype),
        dialect="standard",
        use_bqstorage_api=("BIGQUERY_STORAGE_API_ENABLED" in os.environ),
        progress_bar_type="tqdm_notebook")
    data = pd.merge(ptdata, demog, on="person_id", how="inner")
    MSG("After merge, have %s data points"%data.shape[0])

    # Restrict to samples we want to keep
    if args.samples is not None:
        ERROR("Sample subsetting not implemented yet") # TODO
    MSG("After filter samples, have %s data points"%data.shape[0])

    # Filtering
    data.dropna(axis=0, subset=['value_as_number'],inplace=True)
    MSG("After filter NA, have %s data points"%data.shape[0])
    data = data[data["unit_concept_name"].isin(aou_queries.GetUnits(args.phenotype))]
    MSG("After filter units, have %s data points"%data.shape[0])
    minval, maxval = aou_queries.GetPhenotypeRange(args.phenotype)
    if minval is None or maxval is None:
        ERROR("No minval or maxval specified")
    data = data[(data["value_as_number"]>minval) & (data["value_as_number"]<maxval)]
    MSG("After filter range, have %s data points"%data.shape[0])

    # Determine a single representative value per person
    data['Year'] = data['measurement_datetime'].dt.strftime('%Y')
    median_per_year = data.groupby(['person_id','Year']).agg(median_year=('value_as_number', np.median)).reset_index()
    median_of_medians = median_per_year.groupby(['person_id']).agg(median_median=('median_year', np.median)).reset_index()
    median_of_medians.rename({"median_median": "value_as_number"}, inplace=True, axis=1)

    # Merge back to whole dataframe to only keep median of median value per person
    filtered = pd.merge(data, median_of_medians, on=["person_id", "value_as_number"])

    # De-duplicate to keep one entry per person
    filtered.sort_values("measurement_datetime").drop_duplicates(subset=["person_id"], keep="last", inplace=True)

    # TODO:
    # get age
    # output final file with person_id, phenotype, age at measurement
    print(filtered.head())

if __name__ == "__main__":
    main()