"""
Classes for performing GWAS
"""

import numpy as np

MT_WGS_PATH = 'gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/acaf_threshold/multiMT/hail.mt' 
SMALLNUM = 10e-400

class HailRunner:
    import hail as hl # only import hail if we have to
    def __init__(self, ptcovar, region=None, covars=None):
        self.ptcovar = ptcovar
        self.region = region
        self.covars = covars
        self.gwas = None
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

        # Locus and Sample QC
        data = self.hl.variant_qc(data)
        data = self.hl.sample_qc(data)
        data = data.filter_cols(data.sample_qc.call_rate >= .90, keep = True)
        data = data.filter_rows(data.variant_qc.call_rate >= .90, keep = True)
        data = data.filter_rows(self.hl.min(data.variant_qc.AF) > 0.01, keep = True)
        data = data.filter_rows(data.variant_qc.p_value_hwe > 1e-15, keep = True)


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