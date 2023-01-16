version 1.0

workflow run_hipstr {
    input {
        File bam
        File bam_index
        File genome
        File str_ref
        String out_prefix
    }

    call hipstr {
        input : 
          bam=bam, 
          bam_index=bam_index,
          genome=genome, 
          str_ref=str_ref,
          out_prefix=out_prefix
    }

    call sort_index {
      input :
        vcf=hipstr.outfile
    }

    output {
       File outfile = sort_index.outvcf 
       File outfile_index = sort_index.outvcf_index
    }
    meta {
      description: "Run HipSTR on a single sample with default parameters"
    }
}

task hipstr {
    input {
        File bam
        File bam_index
        File genome
        File str_ref
        String out_prefix
    } 

    command <<<
      HipSTR \
          --bams ~{bam} \
          --fasta ~{genome} \
          --regions ~{str_ref} \
          --str-vcf ~{out_prefix}.vcf.gz
    >>>
    
    runtime {
        docker:"mgymrek/hipstr-gymreklab:latest"
    }

    output {
       File outfile = "${out_prefix}.vcf.gz"
    }
}

task sort_index {
  input {
    File vcf
  }

  String basename = basename(vcf, ".vcf")

  command <<<
    vcf-sort ~{vcf} | bgzip -c > ~{basename}.vcf.gz && tabix -p vcf ~{basename}.sorted.vcf.gz
  >>>

  runtime {
        docker:"mgymrek/vcfutils:latest"
    }

  output {
    File outvcf = "${basename}.vcf.gz"
    File outvcf_index = "${basename}.vcf.gz.tbi"
  }
}
