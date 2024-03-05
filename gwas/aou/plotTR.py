#!/usr/bin/env python3

"""
Make a plot summarizing TR/trait association and allele freqs

Example:
./plotTR.py --phenotype platelet_count --tr-vcf ${WORKSPACE_BUCKET}/cromwell-execution/targetTR/4519d903-dde6-4c8e-b1ee-bcc2d7cd6dd7/call-sort_index/CBL_test.filtered.sorted.vcf.gz --region chr11:119206290-119206323 --outprefix CBL_platelet_count
"""

import argparse
import os
import pandas as pd
import trtools.utils.tr_harmonizer as trh
import trtools.utils.utils as utils
from utils import MSG, ERROR

SAMPLEFILE = os.path.join(os.environ["WORKSPACE_BUCKET"], "samples", \
    "passing_samples_v7.csv")

def GetPTCovarPath(phenotype):
    return os.path.join(os.getenv('WORKSPACE_BUCKET'), \
        "phenotypes", "%s_phenocovar.csv"%phenotype)

def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--phenotype", help="Phenotypes file path, or phenotype name", type=str, required=True)
    parser.add_argument("--samples", help="List of sample IDs, sex to keep", type=str, default=SAMPLEFILE)
    parser.add_argument("--tr-vcf", help="VCF file with TR genotypes.", type=str, required=True)
    parser.add_argument("--region", help="chr:start-end of target TR", type=str, required=True)
    parser.add_argument("--outprefix", help="Output prefix", type=str, required=True)
    args = parser.parse_args()

    # Load phenotype data
    if args.phenotype.endswith(".csv"):
        ptcovar_path = args.phenotype
    else:
        ptcovar_path = GetPTCovarPath(args.phenotype)
    data = pd.read_csv(ptcovar_path)
    data["person_id"] = data["person_id"].apply(str)

    # Load TR genotypes
    invcf = utils.LoadSingleReader(args.tr_vcf, checkgz=True)
    samples = invcf.samples
    region = invcf(args.region)
    nrecords = 0
    for record in region:
        trrecord = trh.HarmonizeRecord(vcftype, record)
        afreqs = trrecord.GetAlleleFreqs()
        genotypes = trrecord.GetLengthGenotypes()
        dosages = [sum(item)/2 for item in genotypes]
        trdf = pd.DataFame({"person_id": samples, "tr_dosage": dosages})
        nrecords += 1
    if nrecords == 0:
        ERROR("No matching TR records found")
    if nrecords > 1:
        ERROR("Multiple matching TR records found")

    # Merge phenotype and TR dosages
    df = pd.merge(data, trdf, on=["person_id"])
    print(df.head())

    # Plot - TODO
if __name__ == "__main__":
    main()