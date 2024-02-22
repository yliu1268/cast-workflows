#!/usr/bin/env python3

"""
Run GWAS on All of Us

Example:
./aou_gwas.py --phenotype ALT --num-pcs 10 --region chr11:119206339-119308149
"""

import argparse
from gwas_runners import HailRunner
from gwas_plotter import PlotManhattan, PlotQQ
import os
import pandas as pd
import re
import sys
from utils import MSG, ERROR
import numpy as np
import scipy.stats as stats
import warnings

GWAS_METHODS = ["hail"]
ANCESTRY_PRED_PATH = "gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv"
SAMPLEFILE = os.path.join(os.environ["WORKSPACE_BUCKET"], "samples", \
    "passing_samples_v7.csv")

def GetPTCovarPath(phenotype):
    return os.path.join(os.getenv('WORKSPACE_BUCKET'), \
        "phenotypes", "%s_phenocovar.csv"%phenotype)

def CheckRegion(region):
    if region is None: return True
    return re.match(r"\w+:\d+-\d+", region) is not None

def GetOutPath(phenotype, method, region):
    outprefix = "%s_%s"%(phenotype, method)
    if region is not None:
        outprefix += "_%s"%(region.replace(":", "_").replace("-","_"))
    return outprefix + ".gwas"

def GetFloatFromPC(x):
    x = x.replace("[","").replace("]","")
    return float(x)

def LoadAncestry():
    if not os.path.isfile("ancestry_preds.tsv"):
        os.system("gsutil -u ${GOOGLE_PROJECT} cp %s ."%(ANCESTRY_PRED_PATH))
    ancestry = pd.read_csv("ancestry_preds.tsv", sep="\t")
    ancestry.rename({"research_id": "person_id"}, axis=1, inplace=True)
    num_pcs = len(ancestry["pca_features"].values[0].split(","))
    pcols = ["PC_%s"%i for i in range(num_pcs)]
    ancestry[pcols] = ancestry["pca_features"].str.split(",", expand=True)
    for p in pcols:
        ancestry[p] = ancestry[p].apply(lambda x: GetFloatFromPC(x), 1)
    return ancestry

#add normalization
def Inverse_Quantile_Normalization(M):   
    R = stats.mstats.rankdata(M,axis=1)  # ties are averaged
    Q = stats.norm.ppf(R/(M.shape[1]+1))
    return Q


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
    parser.add_argument("--region", help="chr:start-end to restrict to. Default is genome-wide", type=str)
    parser.add_argument("--num-pcs", help="Number of PCs to use as covariates", type=int, default=10)
    parser.add_argument("--ptcovars", help="Comma-separated list of phenotype-specific covariates. Default: age", type=str, default="age")
    parser.add_argument("--sharedcovars", help="Comma-separated list of shared covariates (besides PCs). Default: sex_at_birth_Male", type=str, default="sex_at_birth_Male")
    parser.add_argument("--plot", help="Make a Manhattan plot", action="store_true")
    parser.add_argument("--normalization", help="normalize phenotype either quantile or z-score",type=str,default="quantile")
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

    # Get covarlist
    pcols = ["PC_%s"%i for i in range(1, args.num_pcs+1)]
    shared_covars = [item for item in args.sharedcovars.split(",") if item != ""]
    pt_covars = [item for item in args.ptcovars.split(",") if item != ""]
    covars = pcols + pt_covars + shared_covars

    # Set up data frame with phenotype and covars
    data = pd.read_csv(ptcovar_path)
    ancestry = LoadAncestry()
    data = pd.merge(data, ancestry[["person_id"]+pcols], on=["person_id"])
    data["person_id"] = data["person_id"].apply(str)

    # Add normalization quantile
    normalization = args.normalization
    if normalization == "quantile":
        normalize = Inverse_Quantile_Normalization(data[["phenotype"]].transpose()).transpose()
        data["normalized_value"] = normalize.tolist()
        data["phenotype"] = data["normalized_value"].apply(lambda x: ','.join(map(str, x)))
        data["phenotype"] = data["phenotype"].astype(float)

    print(data["phenotype"].head())
    #add normalization z score

  


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
        runner = HailRunner(data, region=args.region, covars=covars)
    else:
        ERROR("GWAS method %s not implemented")

    # Run GWAS
    outpath = GetOutPath(args.phenotype, args.method, args.region)
    runner.RunGWAS()
    WriteGWAS(runner.gwas, outpath+".tab", covars)

    # Plot Manhattan
    if args.plot:
        PlotManhattan(runner.gwas, outpath+".manhattan.png")
        PlotQQ(runner.gwas, outpath+".qq.png")

if __name__ == "__main__":
    main()