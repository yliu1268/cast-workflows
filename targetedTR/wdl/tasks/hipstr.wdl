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

    output {
       File outfile = hipstr.outfile 
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
        docker:"yli091230/hipstr:small"
    }

    output {
       File outfile = "${out_prefix}.vcf.gz"
    }
}
