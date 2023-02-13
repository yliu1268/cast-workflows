version 1.0

workflow trait_association {
	input {
		File trvcf
		File trvcf_idx
		String phenotype
		Array[String] covars
		File phenocovars
	}

	### Set up npy files ###
	call setup_associatr_npy {
		input :
			phenotype=phenotype,
			covars=covars,
			phenocovars=phenocovars
	}

	### Output files ####
	output {
		# TODO
		File ptfile = setup_associatr_npy.pt_file
		File covarfile = setup_associatr_npy.covar_file
	}
}

task setup_associatr_npy {
	input {
		String phenotype
		Array[String] covars
		File phenocovars
	}

	command {
		echo '~{sep="\n" covars}' > covar.list
		python <<CODE
		import pandas as pd
		import numpy as np
		covarlist = [item.strip() for item in open("covar.list").readlines()]
		print('${phenotype}')
		print(str(covarlist))
		df = pd.read_csv('${phenocovars}', sep="\t")
		np.save('${phenotype}'+"_phenotype.npy", \
			df[["IID", '${phenotype}']].to_numpy())
		np.save('${phenotype}'+"_covars.npy", \
			df[["IID"]+covarlist].to_numpy())
		CODE
	}

	runtime {
      docker: "quay.io/biocontainers/pandas:1.4.3"
  	}

	output {
		File pt_file = "${phenotype}_phenotype.npy"
		File covar_file = "${phenotype}_covars.npy"
	}

	meta {
      description: "Convert phenotype/covars to associatr npy format"
    }
}