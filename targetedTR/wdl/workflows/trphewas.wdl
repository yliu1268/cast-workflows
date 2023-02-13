version 1.0

# TODO
# Run associatr to test for association of a phenotype
# Visualize TR vs. phenotype

workflow trphewas {
	input {
		File trvcf
		File trvcf_idx
		File phenocovars
		Array[Array[String]] covars
		Array[String] phenotypes
		String outname
	}

	### Run association testing/vis on each phenotype ###
	scatter(i in range(length(phenotypes))) {
		String use_pt = phenotypes[i]
		Array[String] use_covars = covars[i]
		### Run associaTR ###
		call associaTR {
			input :
				trvcf = trvcf,
				trvcf_idx = trvcf_idx,
				phenotype = use_pt,
				covars = use_covars,
				phenocovars = phenocovars
		}
		### Run visualization ###
		call visTR {
			input :
				trvcf = trvcf,
				trvcf_idx = trvcf_idx,
				phenotype = use_pt,
				covars = use_covars,
				phenocovars = phenocovars
		}
	}

	### Gather outputs ###
	# TODO

	### Output files ####
	output {
		# TODO
	}
}