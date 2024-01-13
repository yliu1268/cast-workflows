#!/usr/bin/env python3
"""
Usage: ./wrangle_ptcovars.py <covarlist> <phenotype> <ptfile>
"""

import pandas as pd
import numpy as np
import sys

try:
	covarlist = sys.argv[1]
	phenotype = sys.argv[2]
	ptfile = sys.argv[3]
except:
	sys.stderr.write(__doc__)
	sys.exit(1)\

covarlist = [item.strip() for item in covarlist.split(",")]
df = pd.read_csv(ptfile, sep="\t")
np.save(phenotype+"_phenotype.npy", \
	df[["IID", phenotype]].to_numpy())
np.save(phenotype+"_covars.npy", \
	df[["IID"]+covarlist].to_numpy())

sys.exit(0)