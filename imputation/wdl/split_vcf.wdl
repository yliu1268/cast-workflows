version 1.0


workflow split_vcf {
    input {
        File vcf
        File vcf_index
        String out_prefix
        String GOOGLE_PROJECT = ""
        String GCS_OAUTH_TOKEN = ""
    }

    call get_sample_list {
        input : 
          vcf=vcf, 
          vcf_index=vcf_index,
          out_prefix=out_prefix,
          GOOGLE_PROJECT=GOOGLE_PROJECT,
          GCS_OAUTH_TOKEN=GCS_OAUTH_TOKEN
    }

    output {
        File outfile = get_sample_list.outfile 
    }
    meta {
      description: "Extract sample list with default parameters"
    }
}

task get_sample_list {
    input {
        File vcf
        File vcf_index
        String out_prefix
        String GOOGLE_PROJECT = ""
        String GCS_OAUTH_TOKEN = ""
    } 

    command <<<
      export GCS_REQUESTER_PAYS_PROJECT=~{GOOGLE_PROJECT}
      export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
      bcftools query -l ~(vcf) > ~(out_prefix).txt
    >>>
    
  
    runtime {
        docker:"nichole/bcftoold-gcs"
    }

    output {
       File outfile = "${out_prefix}.txt"
    }
}



