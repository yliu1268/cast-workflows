version 1.0

workflow trait_association {
	input {
		File trvcf
		File trvcf_idx
		String phenotype
		Array[String] covars
		File phenocovars
		String name
	}

	### Set up npy files ###
	call setup_associatr_npy {
		input :
			phenotype=phenotype,
			covars=covars,
			phenocovars=phenocovars
	}

	### Run associatr ###
	call associatr {
		input :
			ptfile=setup_associatr_npy.pt_file,
			covarfile=setup_associatr_npy.covar_file,
			trvcf=trvcf,
			trvcf_idx=trvcf_idx,
			phenotype=phenotype,
			name=name
	}

	### Visualize genotype vs. phenotype ###
	call viz_assoc {
		input :
			phenocovars=phenocovars,
			covars=covars,
			trvcf=trvcf,
			trvcf_idx=trvcf_idx,
			phenotype=phenotype,
			name=name
	}

	### Output files ####
	output {
		File assocfile = associatr.assocfile
		File visfile = viz_assoc.visfile
	}
}

task setup_associatr_npy {
	input {
		String phenotype
		Array[String] covars
		File phenocovars
	}

	command <<<
		echo '~{sep="\n" covars}' > covar.list
		python <<CODE
		import pandas as pd
		import numpy as np
		import sys
		covarlist = [item.strip() for item in open("covar.list").readlines()]
		print('${phenotype}')
		print(str(covarlist))
		df = pd.read_csv('${phenocovars}', sep="\t")
		np.save('${phenotype}'+"_phenotype.npy", \
			df[["IID", '${phenotype}']].to_numpy())
		np.save('${phenotype}'+"_covars.npy", \
			df[["IID"]+covarlist].to_numpy())
		sys.exit(0)
		CODE
	>>>

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

task associatr {
	input {
		File trvcf
		File trvcf_idx
		String phenotype
		File ptfile
		File covarfile
		String name
	}

	command <<<
	associaTR \
		~{phenotype}_~{name}_associations.tsv \
		~{trvcf} \
		~{phenotype} \
		~{ptfile} ~{covarfile}
	>>>

	runtime {
      docker: "gcr.io/ucsd-medicine-cast/trtools-5.0.1:latest"
  	}

  	output {
  		File assocfile = "${phenotype}_${name}_associations.tsv"
  	}

	meta {
		description: "Run associatr on a phenotype"
	}
}

task viz_assoc {
	input {
		String phenotype
		Array[String] covars
		File phenocovars
		File trvcf
		File trvcf_idx
		String name
	}

	command <<<
	trviz.py \
		--phenotype ~{phenotype} \
		--covars ~{sep="," covars} \
		--ptfile ~{phenocovars} \
		--trvcf ~{trvcf} \
		--out ~{name}_~{phenotype}_viz.png
	>>>

	output {
		File visfile = "${name}_${phenotype}_viz.png"
	}

	runtime {
      docker: "gcr.io/ucsd-medicine-cast/trviz:latest"
  	}

	meta {
		description: "Visualize TR length vs. phenotype" 
	}
}