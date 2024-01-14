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

def PlotAssocStratified(df, ax, ptcol):
	alleles = set(list(df["GT1"]) + list(df["GT2"]))
	for a in alleles:
		if a < 0: continue
		adata = df[(df["GT1"]==a) | (df["GT2"]==a)].copy()
		if adata.shape[0] < 500: continue
		adata["other_allele"] = adata.apply(lambda x: x["GT2"] if x["GT1"]==a else x["GT1"], 1)
		aggdata = adata.groupby(["other_allele"], as_index=False).agg({ptcol: np.mean})
		aggdata_num = adata.groupby(["other_allele"], as_index=False).agg({ptcol: len})
		aggdata_num = aggdata_num.rename({ptcol: "num.samples"}, axis=1)
		aggdata = pd.merge(aggdata, aggdata_num, on="other_allele")
		aggdata = aggdata[aggdata["num.samples"]>100].sort_values("other_allele")
		ax.plot(aggdata["other_allele"], aggdata[ptcol], label="A1=%s"%a)
	ax.set_xlabel("rpt. len")
	ax.set_ylabel(ptcol)
	ax.legend()

def PlotAssoc(df, ax, ptcol):
	aggdata_mean = df.groupby(["GTLEN"], as_index=False).agg({ptcol: np.mean})
	aggdata_mean = aggdata_mean.rename({ptcol: "mean"}, axis=1)
	aggdata_se = df.groupby(["GTLEN"], as_index=False).agg({ptcol: lambda x: np.sqrt(np.var(x)/len(x))})
	aggdata_se = aggdata_se.rename({ptcol: "se"}, axis=1)
	aggdata_num = df.groupby(["GTLEN"], as_index=False).agg({ptcol: len})
	aggdata_num = aggdata_num.rename({ptcol: "num.samples"}, axis=1)
	pltdata = pd.merge(aggdata_mean, aggdata_se, on="GTLEN").sort_values("GTLEN")
	pltdata = pd.merge(pltdata, aggdata_num, on="GTLEN").sort_values("GTLEN")
	pltdata = pltdata[pltdata["GTLEN"]> 0]
	pltdata = pltdata[pltdata["num.samples"]>=100]
	ax.errorbar(pltdata["GTLEN"], pltdata["mean"], yerr=pltdata["se"], \
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

	# Visualization of raw pt and adjusted pt
	fig = plt.figure()
	ax = fig.add_subplot(2, 2, 1)
	PlotAssoc(df, ax, args.phenotype)
	ax = fig.add_subplot(2, 2, 2)
	PlotAssoc(df, ax, "%s.resid"%args.phenotype)
	ax = fig.add_subplot(2, 2, 3)
	PlotAssocStratified(df, ax, args.phenotype)
	ax = fig.add_subplot(2, 2, 4)
	PlotAssocStratified(df, ax, "%s.resid"%args.phenotype)
	fig.tight_layout()
	fig.savefig(args.out)

if __name__ == "__main__":
	main()
