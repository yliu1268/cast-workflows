version 1.0

workflow merge_SNP_batch {
    input {
        Array[File] vcfs
        Array[File] vcfs_index
        String out_prefix
    }

    call mergeSNP {
        input:
            vcfs=vcfs,
            vcfs_index=vcfs_index,
            out_prefix=out_prefix
    }

    call index_vcf {
        input:
            vcf=mergeSNP.outfile
    }

    output {
        File outfile = index_vcf.outvcf 
        File outfile_index = index_vcf.outvcf_index
    }

    meta {
        description: "Merge batches of vcf file containing SNP across a single chromosome"
    }
}

task mergeSNP  {
    input {
        Array[File] vcfs
        Array[File] vcfs_index
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

task index_vcf {
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
