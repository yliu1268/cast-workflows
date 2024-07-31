version 1.0

import "gnomix.wdl" as gnomix_t

# TODO
# do we need to give a reference panel to Beagle?
# Sigh the gnomix models are in hs37! Need to liftover...
# At least can then use 1000G ref panel for Beagle
# Need to unpack gnomix model
# Command to merge gnomix outputs

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
        echo "TODO"
    >>>

    output {
        File msp_outfile = "~{out_prefix}.msp"
        File fb_outfile = "~{out_prefix}.fb"
    }
}