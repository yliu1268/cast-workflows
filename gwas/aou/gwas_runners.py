"""
Classes for performing GWAS
"""

import hail as hl

MT_WGS_PATH = 'gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/acaf_threshold/multiMT/hail.mt' 

class HailRunner:
    def __init__(self, ptcovar, region=None, covars=None):
        self.ptcovar = ptcovar
        self.region = region
        self.gwas_results = None
        self.setup()

    def setup(self):
        hl.init(default_reference = "GRCh38")

        # Load genotypes
        mt = hl.read_matrix_table(MT_WGS_PATH)
        if region is not None:
            mt = h.filter_intervals(mt, [hl.parse_locus_interval(region,)])

        # Load phenotype and covariates
        ptcovar = hl.Table.from_pandas(self.ptcovar, key="persion_id")
        data = mt.annotate_cols(ptcovar = ptcovar[mt.s])

        # Genotype QC
        data = data.annotate_entries(FT = hl.coalesce(data.FT,'PASS'))
        data = data.filter_entries(data.FT =='PASS')
        data = data.filter_entries(data.GQ >= 20)

         # Keep track of data
        self.data = data

    def RunGWAS(self):
        linear_r = hl.linear_regression_rows(
            y= data.ptcovar.phenotype,
            x= data.GT.n_alt_alleles(),
            covariates = self.covars
        )
        summary = linear_r.flatten()
        summary.export(self.out_path)

        # Convert to data frame with required columns
        #gwas = linear_r.annotate(p_value_str= hl.str(linear_r.p_value))
        #gwas_pd = gwas.to_pandas()
        #gwas_pd["chrom"] = gwas_pd["locus"].apply(lambda x: str(x).split(":")[0])
        #gwas_pd["pos"] = gwas_pd["locus"].apply(lambda x: int(str(x).split(":")[1]))
        #gwas_pd['p_value']= gwas_pd.apply(lambda x: change_p(float(x['p_value_str'])),axis=1)
        #gwas_pd['-log10pvalue'] = -np.log10(gwas_pd.p_value)
        #self.gwas_results = gwas_pd

