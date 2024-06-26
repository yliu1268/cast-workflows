version 1.0

workflow merge_batch {
    input {
        Array[File] vcfs
        Array[File] vcfs_index
        String out_prefix
    }

    call mergeSTR {
        input:
            vcfs=vcfs,
            vcfs_index=vcfs_index,
            out_prefix=out_prefix
    }

   # call index_vcf {
   #     input:
   #         vcf=mergeSTR.outfile
   # }

    output {
        File outfile = mergeSTR.outvcf 
        File outfile_index = mergeSTR.outvcf_index
    }

    meta {
        description: "Merge batches of vcf file across a single chromosome"
    }
}

task mergeSTR {
    input {
        Array[File] vcfs
        Array[File] vcfs_index
        String out_prefix
        }

    command <<<
        mergeSTR --vcfs ~{sep=',' vcfs} --out ~{out_prefix}_merged --vcftype hipstr
        vcf-sort ~{out_prefix}_merged.vcf | bgzip -c > ~{out_prefix}_merged.vcf.gz && tabix -p vcf ~{out_prefix}_merged.vcf.gz
    >>>
    
    runtime {
        docker: "gcr.io/ucsd-medicine-cast/trtools-6.0.1:latest"
        disks: "local-disk 20 SSD"
    }
    output {
        File outvcf = "${out_prefix}_merged.vcf.gz"
        File outvcf_index = "${out_prefix}_merged.vcf.gz.tbi"
    }

}

#task index_vcf {
#    input {
#      File vcf
#    }
#    String basename = basename(vcf, ".vcf")
#
#    command <<<
#        vcf-sort ~{vcf} | bgzip -c > ~{basename}.sorted.vcf.gz && tabix -p vcf ~{basename}.sorted.vcf.gz
#      
#    >>>
#
#   runtime {
#        docker:"gcr.io/ucsd-medicine-cast/vcfutils:latest"
#    }
#
#   output {
#    File outvcf = "${basename}.sorted.vcf.gz"
#    File outvcf_index = "${basename}.sorted.vcf.gz.tbi"
#  }
#}
