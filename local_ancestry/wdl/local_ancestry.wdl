version 1.0

import "gnomix.wdl" as gnomix_t

workflow local_ancestry {
    input {
        String out_prefix
        String multi_sample_vcf 
        Array[File] samples
        File gnomix_model
        String chrom
        String GOOGLE_PROJECT = ""
    }

    # Call Beagle/gnomix on batches
    Int num_batches = length(samples)
    scatter (i in range(num_batches)) {
        File sample_batch = samples[i]
        call gnomix_t.run_gnomix as run_gnomix {
            input:
                samples=sample_batch,
                multi_sample_vcf=multi_sample_vcf,
                model=gnomix_model,
                chrom=chrom,
                out_prefix=out_prefix+".BATCH"+i,
                GOOGLE_PROJECT=GOOGLE_PROJECT
        }
    }

    # Merge the output
    call merge_gnomix {
        input:
            gnomix_outputs_msp=run_gnomix.msp_outfile,
            gnomix_outputs_fb=run_gnomix.fb_outfile,
            out_prefix=out_prefix
    }

    ### Output files ####
    output {
        File msp_outfile = merge_gnomix.msp_outfile
        File fb_outfile = merge_gnomix.fb_outfile
    }
}

task merge_gnomix {
    input {
        Array[File] gnomix_outputs_msp
        Array[File] gnomix_outputs_fb
        String out_prefix
    }

    command <<<
        # TODO remove duplicate columns
        paste ~{gnomix_outputs_msp} > ~{out_prefix}.msp
        paste ~{gnomix_outputs_fb} > ~{out_prefix}.fb
    >>>

    output {
        File msp_outfile = "~{out_prefix}.msp"
        File fb_outfile = "~{out_prefix}.fb"
    }
}