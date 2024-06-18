version 1.0

import "imputation.wdl" as imputation_t
import "processTR.wdl" as processTR_t
import "merge_batch.wdl" as merge_batch_t

workflow batch_imputation {
        input {
                String vcf 
                String vcf_index 
                File ref_panel
                File ref
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
                call processTR_t.processTR as processTR {
                    input:
                        vcf=run_imputation.outfile,
                        vcf_index=run_imputation.outfile_index,
                        ref=ref,
                        out_prefix=out_prefix,
                        header_file=header_file

                }
        }
        call merge_batch_t.merge_batch as merge_batch {
            input:
                vcfs=processTR.outfile,
                vcfs_index=processTR.outfile_index,
                out_prefix=out_prefix
        }
     
        output {
            File finalvcf = merge_batch.outfile 
            File finalvcf_index = merge_batch.outfile_index
        }
        meta {
            description: "This workflow run imputation on batches of sample, extract TRs and merge across a single chromosome with default parameters "
    }
            
        
}

