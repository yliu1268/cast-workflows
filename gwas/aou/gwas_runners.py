"""
Classes for performing GWAS
"""

MT_WGS_PATH = 'gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/acaf_threshold/multiMT/hail.mt' 

class HailRunner:
    import hail as hl # only import hail if we have to
    def __init__(self, ptcovar, region=None, covars=None):
        self.ptcovar = ptcovar
        self.region = region
        self.covars = covars
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
            covariates = [1.0 + self.data.ptcovar[item] \
            	for item in self.covars]
        )
        gwas = linear_r.annotate(p_value_str= self.hl.str(linear_r.p_value)).to_pandas()
        print(gwas.head())

        # TODO - export to dataframe and save to the class

        # Convert to data frame with required columns
        #gwas_pd["chrom"] = gwas_pd["locus"].apply(lambda x: str(x).split(":")[0])
        #gwas_pd["pos"] = gwas_pd["locus"].apply(lambda x: int(str(x).split(":")[1]))
        #gwas_pd['p_value']= gwas_pd.apply(lambda x: change_p(float(x['p_value_str'])),axis=1)
        #gwas_pd['-log10pvalue'] = -np.log10(gwas_pd.p_value)
        #self.gwas_results = gwas_pd

