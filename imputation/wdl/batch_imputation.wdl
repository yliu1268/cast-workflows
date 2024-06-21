version 1.0

import "imputation.wdl" as imputation_t
import "processTR.wdl" as processTR_t
import "processSNP.wdl" as processSNP_t
import "merge_batch.wdl" as merge_batch_t


workflow batch_imputation {
        input {
                String vcf 
                String vcf_index 
                File ref_panel
                File ref
                File ref_index
                String out_prefix
                String GOOGLE_PROJECT = ""
                Int? mem 
                Int? window_size 
                String? region
                Boolean subset_region = false
                Boolean beagle_region = false
                Array[File] samples = []
                File header_file

        }
    ### Call subsetting samples with batches ###
        Int num_batches = length(samples)
        scatter (i in range(num_batches)) {
                File sample_batch = samples[i]
                call imputation_t.run_imputation as run_imputation {
                    input:
                        sample=sample_batch,
                        vcf=vcf,
                        vcf_index=vcf_index,
                        ref_panel=ref_panel,
                        region=region,
                        GOOGLE_PROJECT=GOOGLE_PROJECT,
                        subset_region=subset_region,
                        beagle_region=beagle_region,
                        out_prefix=out_prefix+".BATCH"+i,
                        mem=mem,
                        window_size=window_size
                }
                ## extract TR from batches of beagle output
                call processTR_t.processTR as processTR {
                    input:
                        vcf=run_imputation.outfile,
                        vcf_index=run_imputation.outfile_index,
                        ref=ref,
                        ref_index=ref_index,
                        out_prefix=out_prefix,
                        header_file=header_file

                }
                ## extract SNP from batches of beagle output
                call processSNP_t.processSNP as processSNP {
                    input:
                        vcf=run_imputation.outfile,
                        vcf_index=run_imputation.outfile_index,
                        out_prefix=out_prefix

                }
        }

        ## use MergeSTR to merge TR
        call merge_batch_t.merge_batch as merge_batch {
            input:
                vcfs=processTR.outfile,
                vcfs_index=processTR.outfile_index,
                out_prefix=out_prefix
        }

        ## use bcftools to merge SNP
        call merge_SNP {
            input:
                vcfs=processSNP.outfile,
                vcfs_index=processSNP.outfile_index,
                out_prefix=out_prefix

        }

        call sort_index {
            input:
                vcf=merge_SNP.outfile
        }

        output {
            File trvcf = merge_batch.outfile 
            File trvcf_index = merge_batch.outfile_index
            File snpvcf = sort_index.outvcf 
            File snpvcf_index = sort_index.outvcf_index
        }
        meta {
            description: "This workflow run imputation on batches of sample, extract TRs and merge across a single chromosome with default parameters "
    }
                   
}

task merge_SNP {
    input {
        Array [File] vcfs
        Array [File] vcfs_index
        String out_prefix
    }

    command <<<
        bcftools merge ~{sep=',' vcfs} -o ~{out_prefix}_merged_SNP.vcf
    >>>

    runtime {
        docker:"gcr.io/ucsd-medicine-cast/bcftools-gcs:latest"
    }

    output {
        File outfile = "${out_prefix}_merged_SNP.vcf"
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
        docker:"gcr.io/ucsd-medicine-cast/vcfutils:latest"
    }

	output {
		File outvcf = "${basename}.sorted.vcf.gz"
		File outvcf_index = "${basename}.sorted.vcf.gz.tbi"
	}
}
