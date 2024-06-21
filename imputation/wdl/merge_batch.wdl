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

    call index_vcf {
        input:
            vcf=mergeSTR.outfile
    }

    output {
        File outfile = index_vcf.outvcf 
        File outfile_index = index_vcf.outvcf_index
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

    >>>
    
    runtime {
        docker: "gcr.io/ucsd-medicine-cast/trtools-6.0.1:latest"
    }
    output {
        File outfile = "${out_prefix}_merged.vcf"
    }

}

task index_vcf {
    input {
      File vcf
    }
    String basename = basename(vcf, ".vcf")

    command <<<
        bgzip ~{vcf} > ~{basename}.vcf.gz 
        tabix -p vcf ~{basename}.vcf.gz
      
    >>>

    runtime {
        docker:"gcr.io/ucsd-medicine-cast/vcfutils:latest"
    }

    output {
    File outvcf = "${basename}.vcf.gz"
    File outvcf_index = "${basename}.vcf.gz.tbi"
  }
}
