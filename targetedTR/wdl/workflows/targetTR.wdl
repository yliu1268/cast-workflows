version 1.0

# TODO
# Implementation:
# Get list of UKB files
# Change to use GCR rather than dockerhub

# TODO
# Workflow improvements:
# Take stutter input, allele set, from 
# multi-sample hipstr calls from UKB subset
# Just take in hipstr ID

import "../tasks/hipstr.wdl" as hipstr_t
import "../tasks/merge_hipstr.wdl" as merge_t

workflow targetTR {
	input {
		String chrom
		Int str_start
		Int str_end
		Int motif_len
		Float num_copies
		String str_name
		Array[File] cram_files
		Array[File] cram_index_files
		File genome
	}

	### Generate HipSTR BED file ###
	call makeBed {
		input :
			chrom=chrom,
			str_start=str_start,
			str_end=str_end,
			motif_len=motif_len,
			num_copies=num_copies,
			str_name=str_name
	}

	### Call HipSTR on all samples ###
	scatter (i in range(length(cram_files))) {
		File cram = cram_files[i]
		File cram_index = cram_index_files[i]
		call hipstr_t.run_hipstr as run_hipstr {
			input :
				bam=cram,
				bam_index=cram_index,
				genome=genome,
				str_ref=makeBed.tr_bed,
				out_prefix=str_name
		}
	}

	### Merge HipSTR across all samples ###
	call merge_t.merge_hipstr as merge_hipstr {
		input :
			vcfs=run_hipstr.outfile,
			out_prefix=str_name
	}

	### DumpSTR and zip on merged VCF ###
	# TODO

	### Output files ####
	# TODO - merged, filtered VCF
	output {
		File outbed = makeBed.tr_bed
		File mergevcf = merge_hipstr.outfile
	}
}

task makeBed {
	input {
		String chrom
		Int str_start
		Int str_end
		Int motif_len
		Float num_copies
		String str_name
	}

	# TODO: couldn't figure out how to get tabs to be escaped properly in
	# the command <<< >>> 
	String command = "echo '~{chrom},~{str_start},~{str_end},~{motif_len},~{num_copies},~{str_name}' | sed 's/,/\t/g' "
	command <<<
		sh -c "~{command}" > ~{str_name}_strref.bed
	>>>

	output {
		File tr_bed = "${str_name}_strref.bed"
	}
}