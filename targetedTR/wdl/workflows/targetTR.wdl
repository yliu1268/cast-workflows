version 1.0

# TODO
# Hipstr mergefix broke with only one record
#  commenting for now, need to bring back
# Adding phewas
# Need: obtain phenotypes
# Need: run associaTR
# Need: Make html report 
#   Length distribution
# TODO
# Workflow improvements:
# Take stutter input, allele set, from 
# multi-sample hipstr calls from UKB subset

import "../tasks/hipstr_multi.wdl" as hipstr_multi_t
import "../tasks/merge_hipstr.wdl" as merge_t
import "../tasks/dumpstr.wdl" as dumpstr_t

workflow targetTR {
	input {
		String str_name
		Array[Array[File]] cram_file_batches
		Array[Array[File]] cram_index_batches
		File genome
		File genome_index
		File tr_bed
	}

	### Call HipSTR on batches of samples ###
	scatter(i in range(length(cram_file_batches))) {
		Array[File] crams = cram_file_batches[i]
		Array[File] cram_indices = cram_index_batches[i]
		call hipstr_multi_t.run_hipstr as run_hipstr {
			input :
				bams=crams,
				bam_indices=cram_indices,
				genome=genome,
				genome_index=genome_index,
				str_ref=tr_bed,
				out_prefix=str_name+".BATCH"+i
		}
	}

	### Merge HipSTR across all samples ###
	call merge_t.merge_hipstr as merge_hipstr {
		input :
			vcfs=run_hipstr.outfile,
			vcf_indexes=run_hipstr.outfile_index,
			out_prefix=str_name
	}

	### DumpSTR on merged VCF ###
	call dumpstr_t.run_dumpstr as dumpstr {
		input :
			vcf=merge_hipstr.outfile,
			out_prefix=str_name+".filtered"
	}

	### Zip and index the VCF ###
	call sort_index {
		input :
			vcf=dumpstr.outfile
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