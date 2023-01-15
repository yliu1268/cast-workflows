version 1.0

# TODO
# Take stutter input, allele set, from 
# multi-sample hipstr calls from UKB subset
# Just take in hipstr ID

import "../tasks/hipstr.wdl" as hipstr_t

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
		call hipstr_t.run_hipstr {
			input :
				bam=cram,
				bam_index=cram_index,
				genome=genome,
				str_ref=makeBed.tr_bed,
				out_prefix=str_name
		}
	}

	### Merge HipSTR across all samples ###
	# TODO

	### DumpSTR on merged VCF ###
	# TODO

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

	command <<<
		echo "~{chrom}\t~{str_start}\t~{str_end}\t~{motif_len}\t~{num_copies}\t~{num_copies}" \
			> ~{str_name}_strref.bed
	>>>

	output {
		File tr_bed = "${str_name}_strref.bed"
	}
}