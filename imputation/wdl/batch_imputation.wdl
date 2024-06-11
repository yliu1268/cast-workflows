version 1.0

import "imputation.wdl" as imputation_t

workflow batch_imputation {
        input {
                String vcf 
                String vcf_index 
                File ref_panel
                String out_prefix
                String GOOGLE_PROJECT = ""
                Int? mem 
                Int? window_size 
                String? region
                Boolean subset_region = false
                Boolean beagle_region = false
                Array[File] samples = []
                Boolean using_batch_files = false

        }
    ### Call subsetting samples with batches ###

        Boolean using_batch_files = (length(samples)>0)
        Int num_batches = length(samples)
        scatter (i in range(num_batches)) {
                File sample_batch = samples [i]
                call imputation_t.run_imputation as run_imputation {
                    input:
                        sample=sample_batch,
                        vcf=vcf,
                        vcf_index=vcf_index,
                        ref_panel=ref_panel,
                        region=region,
                        GOOGLE_PROJECT=GOOGLE_PROJECT,
                        subset_region=subset_region,
                        using_batch_files=using_batch_files,
                        beagle_region=beagle_region,
                        out_prefix=out_prefix+".BATCH"+i,
                        mem=mem,
                        window_size=window_size
                }
        }
        call merge_vcf {
            input:
                vcfs=run_imputation.outvcf,
                vcfs_index=run_imputation.outvcf,
                out_prefix=out_prefix
        }

        
        output {
            File outfile = merge_vcf.outvcf 
            File outfile_index = merge_vcf.outvcf_index
        }
            
        
}

task merge_vcf {
    input {
        Array[File] vcfs
        Array[File] vcfs_index
        String out_prefix
    }

    command <<<
        bcftools merge ${sep=' ' vcfs} -O z -o ${out_prefix}.merged.vcf.gz
        tabix -p vcf ${out_prefix}.merged.vcf.gz

    >>>

    runtime {
        docker:"gcr.io/ucsd-medicine-cast/bcftools-gcs:latest"
    }
    
    output {
        File outvcf = "${out_prefix}.merged.vcf.gz"
        File outvcf_index = "${out_prefix}.merged.vcf.gz.tbi"
    }
}