"""
Classes for performing GWAS
"""

import numpy as np
import os
import pandas as pd
import subprocess
import tempfile

MT_WGS_PATH = 'gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/acaf_threshold/multiMT/hail.mt' 
SMALLNUM = 10e-400

class AssociaTRRunner:
    import trtools.associaTR # only require if we are running TRTools
    def __init__(self, ptcovar, trvcf, region=None, covars=[]):
        self.ptcovar = ptcovar
        self.pt_npy = None
        self.covar_npy = None
        self.trvcf = trvcf
        self.region = region
        self.covars = covars
        self.gwas = None
        self.method = "associaTR"
        self.setup()

    def setup(self):
        # Set up phenotypes and covariates
        self.pt_npy = tempfile.NamedTemporaryFile(suffix='.npy')
        self.covar_npy = tempfile.NamedTemporaryFile(suffix='.npy')
        np.save(self.pt_npy.name, \
            self.ptcovar[["person_id", "phenotype"]].to_numpy())
        np.save(self.covar_npy.name, \
            self.ptcovar[["person_id"]+self.covars].to_numpy())

        # Fetch genotypes file if remote
        if self.trvcf.startswith("gs://"):
            trvcf = sampfile.split("/")[-1]
            if not os.path.isfile(trvcf):
                os.system("gsutil -u ${GOOGLE_PROJECT} cp %s ."%self.trvcf)
            self.trvcf = trvcf

    def RunGWAS(self):
        # Set associaTR options
        vcftype = None
        same_samples = False
        sample_list = None
        non_major_cutoff = 20
        beagle_dosages = False
        plotting_phenotype = None
        paired_genotype_plot = False
        plot_genotype_residuals = False
        plotting_ci_alphas = False
        imputed_ukb_strs_paper_period_check = False
        outfile = tempfile.NamedTemporaryFile(suffix='.tsv')
        trtools.associaTR.perform_gwas(
            outfile.name,
            self.trvcf,
            "phenotype",
            [self.pt_npy, self.covar_npy],
            vcftype,
            same_samples,
            sample_list,
            self.region,
            non_major_cutoff,
            beagle_dosages,
            plotting_phenotypes,
            paired_genotype_plot,
            plot_genotype_residuals,
            plotting_ci_alphas,
            imputed_ukb_strs_paper_period_check
        )
        gwas = pd.read_csv(outfile.name, sep="\t")
        gwas.rename({"p_phenotype": "p_value"}, inplace=True, axis=1)
        gwas["-log10pvalue"] = -np.log10(gwas["p_value"])
        self.gwas = gwas

class HailRunner:
    import hail as hl # only import hail if we have to
    def __init__(self, ptcovar, region=None, covars=[]):
        self.ptcovar = ptcovar
        self.region = region
        self.covars = covars
        self.gwas = None
        self.data = None
        self.method = "hail"
        self.setup()

    def setup(self):
    	# Set up hail
        self.hl.init(default_reference = "GRCh38")

        # Load genotypes
        mt = self.hl.read_matrix_table(MT_WGS_PATH)
        if self.region is not None:
            mt = self.hl.filter_intervals(mt, [self.hl.parse_locus_interval(self.region,)])

        # Load phenotype and covariates
        ptcovar = self.hl.Table.from_pandas(self.ptcovar, key="person_id")
        data = mt.annotate_cols(ptcovar = ptcovar[mt.s])

        # Genotype QC
        data = data.annotate_entries(FT = self.hl.coalesce(data.FT,'PASS'))
        data = data.filter_entries(data.FT =='PASS')
        data = data.filter_entries(data.GQ >= 20)

         # Keep track of data
        self.data = data

    def RunGWAS(self):
        linear_r = self.hl.linear_regression_rows(
            y= self.data.ptcovar.phenotype,
            x= self.data.GT.n_alt_alleles(),
            covariates = [1.0] + [self.data.ptcovar[item] \
            	for item in self.covars]
        )
        gwas = linear_r.annotate(p_value_str= self.hl.str(linear_r.p_value)).to_pandas()
        gwas["chrom"] = gwas["locus"].apply(lambda x: str(x).split(":")[0])
        gwas["pos"] = gwas["locus"].apply(lambda x: int(str(x).split(":")[1]))
        gwas["p_value"] = gwas.apply(lambda x: float(x["p_value_str"])+SMALLNUM, 1)
        gwas["-log10pvalue"] = -np.log10(gwas["p_value"])
        self.gwas = gwas