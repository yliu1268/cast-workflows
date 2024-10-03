version 1.0

workflow batch_imputation {
    input {
        File ref_vcf 
        File ref_index 
        File ref_panel
        File map
        String out_prefix
        Array[File] batch_vcf_files = []
    }

    ### Call subsetting samples with batches ###
    Int num_batches = length(batch_vcf_files)
    scatter (i in range(num_batches)) {
        File batch = batch_vcf_files[i]
        File batch_index = batch_vcf_files[i]+".tbi"
        call beagle {
            input:
                vcf=batch,
                vcf_index=batch_index,
                ref_panel=ref_panel,
                out_prefix=out_prefix+".BATCH"+i,
                map=map
        }
        ## extract TR from batches of beagle output
        call extract_TR {
            input:
                vcf=beagle.outvcf,
                vcf_index=beagle.outvcf_index,
                out_prefix=out_prefix+".BATCH"+i
        }
    }

    ## use bcftools to merge TRs
    call merge_batch {
        input:
            vcfs=extract_TR.outvcf,
            vcfs_index=extract_TR.outvcf_index,
            out_prefix=out_prefix
    }

    call annotaTR {
        input:
            vcf=merge_batch.outvcf,
            vcf_index=merge_batch.outvcf_index,
            ref_vcf=ref_vcf,
            ref_index=ref_index,
            out_prefix=out_prefix
    }

    output {
        File outfile_pgen = annotaTR.pgen
        File outfile_psam = annotaTR.psam
        File outfile_pvar = annotaTR.pvar
        File outfile_vcf = annotaTR.outvcf
        File outfile_vcfind = annotaTR.outvcfind
    }

    meta {
        description: "This workflow runs TR imputation on batches of sample, extract TRs, merges across a single chromosome and runs annotaTR."
    }                 
}

task beagle {
    input {
        File vcf 
        File vcf_index 
        File ref_panel
        String out_prefix
        File map
    } 

    command <<<
        set -e

        java -Xmx25g -jar /beagle.jar \
            gt=~{vcf} \
            ref=~{ref_panel} \
            ap=true \
            out=~{out_prefix}_output \
            map=~{map} 
        tabix -p vcf ~{out_prefix}_output.vcf.gz        
    >>>
    
    runtime {
        docker:"gcr.io/ucsd-medicine-cast/beagle:latest"
        memory: "25G"
        disks: "local-disk 50 SSD"
        preemptible: 1
    }

    output {
       File outvcf = "${out_prefix}_output.vcf.gz"
       File outvcf_index = "${out_prefix}_output.vcf.gz.tbi"
    }
}

task extract_TR {
    input {
        File vcf
        File vcf_index
        String out_prefix
    }
    command <<<
        set -e
        bcftools view -i 'ID~"EnsTR"' ~{vcf} -Oz -o ~{out_prefix}_TR.vcf.gz
        tabix -p vcf ~{out_prefix}_TR.vcf.gz
    >>>

    runtime {
        docker:"gcr.io/ucsd-medicine-cast/bcftools-gcs:latest"
        preemptible: 1
    }

    output {
        File outvcf = "${out_prefix}_TR.vcf.gz"
        File outvcf_index = "${out_prefix}_TR.vcf.gz.tbi"

    }    
}

task merge_batch {
    input {
        Array[File] vcfs
        Array[File] vcfs_index
        String out_prefix
    }

    command <<<
        bcftools merge ~{sep=' ' vcfs} -Oz -o ~{out_prefix}_TR_merged.vcf.gz 
        tabix -p vcf ~{out_prefix}_TR_merged.vcf.gz        
    >>>
    
    runtime {
        docker: "gcr.io/ucsd-medicine-cast/bcftools-gcs:latest"
        memory: "25G"
        disks: "local-disk 60 SSD"
    }
    output {
        File outvcf = "${out_prefix}_TR_merged.vcf.gz"
        File outvcf_index = "${out_prefix}_TR_merged.vcf.gz.tbi"
    }

}

task annotaTR {
    input {
        File vcf 
        File vcf_index
        File ref_vcf
        File ref_index
        String out_prefix
    }
    
    command <<<
        set -e
        annotaTR --vcf ~{vcf} \
            --ref-panel ~{ref_vcf} \
            --out ~{out_prefix}_annotated \
            --vcftype hipstr \
            --outtype pgen vcf \
            --vcf-outtype z \
            --dosages beagleap_norm \
            --ignore-duplicates \
            --match-refpanel-on locid \
            --warn-on-AP-error \
            --update-ref-alt
        tabix -p vcf ~{out_prefix}_annotated.vcf.gz
    >>>

    runtime {
        docker:"gcr.io/ucsd-medicine-cast/trtools-annotatr:latest"
        disks: "local-disk 120 SSD"
        memory: "30G"
    }

    output {
        File pgen = "${out_prefix}_annotated.pgen"
        File psam = "${out_prefix}_annotated.psam"
        File pvar = "${out_prefix}_annotated.pvar"
        File outvcf = "${out_prefix}_annotated.vcf.gz"
        File outvcfind = "${out_prefix}_annotated.vcf.gz.tbi"
    }
}