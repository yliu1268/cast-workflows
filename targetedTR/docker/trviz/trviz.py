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
import matplotlib.pyplot as plt

def RegressOutCovars(df, phenotype, covars):
	Y = np.array(df[phenotype])
	X = df[covars].to_numpy()
	model = OLS(Y, X, missing="drop")
	res = model.fit()
	pred = res.get_prediction(X).predicted_mean
	resid = [Y[i]-pred[i] for i in range(len(Y))]
	return resid

def PlotAssoc(df, ax, ptcol):
	aggdata_mean = df.groupby(["GTLEN"], as_index=False).agg({ptcol: np.mean})
	aggdata_mean = aggdata_mean.rename({ptcol: "mean"}, axis=1)
	aggdata_sd = df.groupby(["GTLEN"], as_index=False).agg({ptcol: lambda x: np.sqrt(np.var(x))})
	aggdata_sd = aggdata_sd.rename({ptcol: "sd"}, axis=1)
	pltdata = pd.merge(aggdata_mean, aggdata_sd, on="GTLEN").sort_values("GTLEN")
	ax.errorbar(pltdata["GTLEN"], pltdata["mean"], yerr=pltdata["sd"], \
		marker="o")
	ax.set_xlabel("Sum of rpt. lens")
	ax.set_ylabel(ptcol)

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

	# Visualization of raw pt and adjusted pt - TODO
	fig = plt.figure()
	ax = fig.add_subplot(1, 2, 1)
	PlotAssoc(df, ax, args.phenotype)
	ax = fig.add_subplot(1, 2, 2)
	PlotAssoc(df, ax, "%s.resid"%args.phenotype)
	fig.tight_layout()
	fig.savefig(args.out)

if __name__ == "__main__":
	main()
