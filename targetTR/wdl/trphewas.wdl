version 1.0

import "../tasks/associatr.wdl" as associatr_t

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
		### Run associaTR and visualization ###
		call associatr_t.trait_association as trait_association {
			input :
				trvcf = trvcf,
				trvcf_idx = trvcf_idx,
				phenotype = use_pt,
				covars = use_covars,
				phenocovars = phenocovars,
				name = outname
		}
	}

	### Gather outputs ###
	# TODO

	### Output files ####
	output {
		# TODO
	}
}