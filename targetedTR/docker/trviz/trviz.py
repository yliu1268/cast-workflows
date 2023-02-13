#!/usr/bin/env python3
"""
Visualize TR vs. phenotype
"""

import argparse
import pandas as pd

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

	with open(args.out, "w") as f:
		f.write("TEST")

if __name__ == "__main__":
	main()
