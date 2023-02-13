#!/usr/bin/env python3
"""
Visualize TR vs. phenotype
"""

import argparse
import numpy as np
import pandas as pd
from statsmodels.regression.linear_model import OLS
import trtools.utils.tr_harmonizer as trh
import trtools.utils.utils as utils

def RegressOutCovars(df, phenotype, covars):
	Y = np.array(df[phenotype])
	X = df[covars].to_numpy()
	model = OLS(Y, X, missing="drop")
	res = model.fit()
	pred = res.get_prediction(X).predicted_mean
	resid = [Y[i]-pred[i] for i in range(len(Y))]
	return resid

def main():
	parser = argparse.ArgumentParser(__doc__)
	parser.add_argument("--phenotype", help="Name of phenotype", \
		type=str, required=True)
	parser.add_argument("--covars", help="Comma-sep list of covars", \
		type=str, required=False)
	parser.add_argument("--ptfile", help="Pheno covars file", \
		type=str, required=True)
	parser.add_argument("--out", help="Output file", \
		type=str, required=True)
	parser.add_argument("--trvcf", help="TR VCF file", \
		type=str, required=True)
	args = parser.parse_args()

	# Load phenotypes
	df = pd.read_csv(args.ptfile, sep="\t")

	# Process only the first TR locus
	invcf = utils.LoadSingleReader(args.trvcf)
	vcftype = trh.InferVCFType(invcf)
	samples = invcf.samples
	for record in invcf:
		trrecord = trh.HarmonizeRecord(vcftype, record)
		gts = trrecord.GetLengthGenotypes()
		gtsum = [sum(item) for item in gts]
		gt1 = [item[0] for item in gts]
		gt2 = [item[1] for item in gts]
		gtdf = pd.DataFrame({"IID": samples, "GTLEN": gtsum,
			"GT1": gt1, "GT2": gt2})
		break # only process one record!
	df["IID"] = df["IID"].apply(str)
	gtdf["IID"] = df["IID"].apply(str)
	df = pd.merge(df, gtdf, on="IID")

    # Get regressed out phenotype
	df["%s.resid"%args.phenotype] = \
		RegressOutCovars(df, args.phenotype, args.covars.split(","))
	print(df)

	# Visualization of raw pt and adjusted pt - TODO
	with open(args.out, "w") as f:
		f.write("TEST")

if __name__ == "__main__":
	main()
