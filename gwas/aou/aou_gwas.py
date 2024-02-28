#!/usr/bin/env python3

"""
Run GWAS on All of Us

Example:
./aou_gwas.py --phenotype ALT --num-pcs 10 --region chr11:119206339-119308149
"""

import argparse
from gwas_plotter import PlotManhattan, PlotQQ
import os
import pandas as pd
import re
import sys
from utils import MSG, ERROR

GWAS_METHODS = ["hail", "associaTR"]
ANCESTRY_PRED_PATH = "gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv"
SAMPLEFILE = os.path.join(os.environ["WORKSPACE_BUCKET"], "samples", \
    "passing_samples_v7.csv")
DEFAULTSAMPLECSV = "passing_samples_v7.csv"

def GetPTCovarPath(phenotype):
    return os.path.join(os.getenv('WORKSPACE_BUCKET'), \
        "phenotypes", "%s_phenocovar.csv"%phenotype)

def CheckRegion(region):
    if region is None: return True
    return re.match(r"\w+:\d+-\d+", region) is not None

def GetOutPath(phenotype, method, region, samplefile):
    cohort = "ALL" #if samplefile == DEFAULTSAMPLECSV
    if samplefile != DEFAULTSAMPLECSV:
        cohort = samplefile[:-4] #chop off .csv at the end
    outprefix = "%s_%s_%s"%(phenotype, method, cohort)
    if region is not None:
        outprefix += "_%s"%(region.replace(":", "_").replace("-","_"))
    return outprefix + ".gwas"

def GetFloatFromPC(x):
    x = x.replace("[","").replace("]","")
    return float(x)

def LoadAncestry(ancestry_pred_path):
    if ancestry_pred_path.startswith("gs://"):
        if not os.path.isfile("ancestry_preds.tsv"):
            os.system("gsutil -u ${GOOGLE_PROJECT} cp %s ."%(ancestry_pred_path))
        ancestry_pred_path = "ancestry_preds.tsv"
    ancestry = pd.read_csv(ancestry_pred_path, sep="\t")
    ancestry.rename({"research_id": "person_id"}, axis=1, inplace=True)
    num_pcs = len(ancestry["pca_features"].values[0].split(","))
    pcols = ["PC_%s"%i for i in range(num_pcs)]
    ancestry[pcols] = ancestry["pca_features"].str.split(",", expand=True)
    for p in pcols:
        ancestry[p] = ancestry[p].apply(lambda x: GetFloatFromPC(x), 1)
    return ancestry

def WriteGWAS(gwas, outpath, covars):
    # Ouptut header with command used
    f = open(outpath, "w")
    f.write("#" + " ".join(sys.argv) + "\n")
    f.write("# covars: %s\n"%",".join(covars))
    f.close()
    # Append gwas results
    gwas[["chrom","pos","beta","standard_error","-log10pvalue"]].to_csv(outpath, sep="\t", mode="a", index=False)

def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--phenotype", help="Phenotypes file path, or phenotype name", type=str, required=True)
    parser.add_argument("--method", help="GWAS method. Options: %s"",".join(GWAS_METHODS), type=str, default="hail")
    parser.add_argument("--samples", help="List of sample IDs, sex to keep", type=str, default=SAMPLEFILE)
    parser.add_argument("--ancestry-pred-path", help="Path to ancestry predictions", default=ANCESTRY_PRED_PATH)
    parser.add_argument("--region", help="chr:start-end to restrict to. Default is genome-wide", type=str)
    parser.add_argument("--num-pcs", help="Number of PCs to use as covariates", type=int, default=10)
    parser.add_argument("--ptcovars", help="Comma-separated list of phenotype-specific covariates. Default: age", type=str, default="age")
    parser.add_argument("--sharedcovars", help="Comma-separated list of shared covariates (besides PCs). Default: sex_at_birth_Male", type=str, default="sex_at_birth_Male")
    parser.add_argument("--tr-vcf", help="VCF file with TR genotypes. Required if running associaTR", type=str)
    parser.add_argument("--plot", help="Make a Manhattan plot", action="store_true")
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
    if args.num_pcs > 10:
        ERROR("Specify a maximum of 10 PCs")
    if args.method == "associaTR" and args.tr_vcf is None:
        ERROR("Must specify --tr-vcf for associTR")

    # Get covarlist
    pcols = ["PC_%s"%i for i in range(1, args.num_pcs+1)]
    shared_covars = [item for item in args.sharedcovars.split(",") if item != ""]
    pt_covars = [item for item in args.ptcovars.split(",") if item != ""]
    covars = pcols + pt_covars + shared_covars

    # Set up data frame with phenotype and covars
    data = pd.read_csv(ptcovar_path)
    ancestry = LoadAncestry(args.ancestry_pred_path)
    data = pd.merge(data, ancestry[["person_id"]+pcols], on=["person_id"])
    data["person_id"] = data["person_id"].apply(str)

    # Add shared covars
    sampfile = args.samples
    if sampfile.startswith("gs://"):
        sampfile = sampfile.split("/")[-1]
        if not os.path.isfile(sampfile):
            os.system("gsutil -u ${GOOGLE_PROJECT} cp %s ."%(args.samples))
    samples = pd.read_csv(sampfile)
    samples["person_id"] = samples["person_id"].apply(str)
    data = pd.merge(data, samples)

    # Check we have all covars
    req_cols = ["phenotype"] + covars
    for item in req_cols:
        if item not in data.columns:
            ERROR("Required column %s not found"%item)

    # Set up GWAS method
    if args.method == "hail":
        from hail_runner import HailRunner
        runner = HailRunner(data, region=args.region, covars=covars)
    elif args.method == "associaTR":
        from associatr_runner import AssociaTRRunner
        runner = AssociaTRRunner(data, args.tr_vcf, region=args.region, covars=covars)
    else:
        ERROR("GWAS method %s not implemented")

    # Run GWAS
    outpath = GetOutPath(args.phenotype, args.method, args.region, sampfile)
    runner.RunGWAS()
    WriteGWAS(runner.gwas, outpath+".tab", covars)

    # Plot Manhattan
    if args.plot:
        PlotManhattan(runner.gwas, outpath+".manhattan.png")
        PlotQQ(runner.gwas, outpath+".qq.png")

if __name__ == "__main__":
    main()
