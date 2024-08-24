version 1.0

import "gnomix.wdl" as gnomix_t

workflow local_ancestry {
    input {
        String out_prefix
        Array[File] batch_vcf_files
        File gnomix_model
        String chrom
        File chainfile
        File refpanel
        File refpanel_index
    }

    # Call Beagle/gnomix on batches
    Int num_batches = length(batch_vcf_files)
    scatter (i in range(num_batches)) {
        File batch_vcf = batch_vcf_files[i]
        call gnomix_t.run_gnomix as run_gnomix {
            input:
                vcf=batch_vcf,
                model=gnomix_model,
                chrom=chrom,
                out_prefix=out_prefix+".BATCH"+i,
                chainfile=chainfile,
                refpanel=refpanel,
                refpanel_index=refpanel_index
        }
    }

    # Merge the output
    call gnomix_t.merge_gnomix as merge_gnomix {
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

