"""
Classes for performing GWAS
"""

MT_WGS_PATH = 'gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/acaf_threshold/multiMT/hail.mt' 

class HailRunner:
    import hail as hl # only import hail if we have to
    def __init__(self, ptcovar, region=None, covars=None):
        self.ptcovar = ptcovar
        self.region = region
        self.gwas_results = None
        self.setup()

    def setup(self):
    	# Set up hail
        self.hl.init(default_reference = "GRCh38")

        # Load genotypes
        mt = self.hl.read_matrix_table(MT_WGS_PATH)
        if self.region is not None:
            mt = hl.filter_intervals(mt, [hl.parse_locus_interval(region,)])

        # Load phenotype and covariates
        ptcovar = self.hl.Table.from_pandas(self.ptcovar, key="person_id")
        data = mt.annotate_cols(ptcovar = ptcovar[mt.s])

        # Genotype QC
        data = data.annotate_entries(FT = self.hl.coalesce(data.FT,'PASS'))
        data = data.filter_entries(data.FT =='PASS')
        data = data.filter_entries(data.GQ >= 20)

         # Keep track of data
        self.data = data

        print(self.data.describe()) # TODO remove

    def RunGWAS(self):
    	# TODO how to access covariates from ptcovar by name
        linear_r = self.hl.linear_regression_rows(
            y= data.ptcovar.phenotype,
            x= data.GT.n_alt_alleles(),
            covariates = [1.0 + data.ptcovar[item] \
            	for item in self.covars]
        )
        summary = linear_r.flatten()
        summary.export(self.out_path)
        # TODO - export to dataframe and save to the class

        # Convert to data frame with required columns
        #gwas = linear_r.annotate(p_value_str= hl.str(linear_r.p_value))
        #gwas_pd = gwas.to_pandas()
        #gwas_pd["chrom"] = gwas_pd["locus"].apply(lambda x: str(x).split(":")[0])
        #gwas_pd["pos"] = gwas_pd["locus"].apply(lambda x: int(str(x).split(":")[1]))
        #gwas_pd['p_value']= gwas_pd.apply(lambda x: change_p(float(x['p_value_str'])),axis=1)
        #gwas_pd['-log10pvalue'] = -np.log10(gwas_pd.p_value)
        #self.gwas_results = gwas_pd

