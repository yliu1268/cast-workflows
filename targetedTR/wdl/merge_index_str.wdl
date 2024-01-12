version 1.0

import "../tasks/merge_hipstr.wdl" as merge_t

workflow merge_index {
	input {
		String out_name
		Array[File] vcf_files
		Array[File] vcf_idxs
	}

	### Merge VCF files ###
	call merge_t.merge_hipstr as merge_hipstr {
		input :
			vcfs=vcf_files,
			vcf_indexes=vcf_idxs,
			out_prefix=out_name
	}

	### Zip and index the VCF ###
	call sort_index {
		input :
			vcf=merge_hipstr.outfile
	}

	### Output files ####
	output {
		File finalvcf = sort_index.outvcf
		File finalvcf_index = sort_index.outvcf_index
	}

}

task sort_index {
	input {
		File vcf
	}

	String basename = basename(vcf, ".vcf")

	command <<<
		vcf-sort ~{vcf} | bgzip -c > ~{basename}.sorted.vcf.gz && tabix -p vcf ~{basename}.sorted.vcf.gz
	>>>

	runtime {
        docker:"mgymrek/vcfutils:latest"
    }

	output {
		File outvcf = "${basename}.sorted.vcf.gz"
		File outvcf_index = "${basename}.sorted.vcf.gz.tbi"
	}
}
