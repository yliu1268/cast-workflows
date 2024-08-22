version 1.0

import "gnomix.wdl" as gnomix_t

workflow local_ancestry {
    input {
        String out_prefix
        Array[File] batch_vcf_files
        File gnomix_model
        String chrom
        String GOOGLE_PROJECT = ""
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
                GOOGLE_PROJECT=GOOGLE_PROJECT,
                chainfile=chainfile,
                refpanel=refpanel,
                refpanel_index=refpanel_index
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
        Int total = length(gnomix_outputs_msp)
    }

    command <<<
        # First merge msp files
        MSPFILEARRAY=(~{sep=' ' gnomix_outputs_msp})
        head -n 1 ${MSPFILEARRAY[0]} > ~{out_prefix}.msp
        cat ${MSPFILEARRAY[0]} | grep -v "^#Subpopulation" | cut -f 1-6 -d$'\t'> fixedcols.msp
        for (( c = 0; c < ~{total}; c++ ))
        do
            cat ${MSPFILEARRAY[$c]} | grep -v "^#Subpopulation" | cut -f 1-6 -d$'\t' --complement > data_${c}.msp
        done
        paste fixedcols.msp data*.msp >> ~{out_prefix}.msp

        # Next merge fb files
        FBFILEARRAY=(~{sep=' ' gnomix_outputs_fb})
        head -n 1 ${FBFILEARRAY[0]} > ~{out_prefix}.fb
        cat ${FBFILEARRAY[0]} | grep -v "^#" | cut -f 1-4 -d$'\t' > fixedcols.fb
        for (( c = 0; c < ~{total}; c++ ))
        do
            cat ${FBFILEARRAY[$c]} | grep -v "^#" | cut -f 1-4 -d$'\t' --complement > data_${c}.fb
        done
        paste fixedcols.fb data*.fb >> ~{out_prefix}.fb
    >>>

    runtime {
        docker: "gcr.io/ucsd-medicine-cast/bcftools-gcs:latest"
    }

    output {
        File msp_outfile = "~{out_prefix}.msp"
        File fb_outfile = "~{out_prefix}.fb"
    }
}