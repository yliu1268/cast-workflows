version 1.0


workflow split_vcf {
    input {
        Array[File] vcf = []
        Array[File] vcf_ndex = []
        String out_prefix
        #String GOOGLE_PROJECT = ""
        #String GCS_OAUTH_TOKEN = ""
    }

    call get_sample_list {
        input : 
          vcf=vcf, 
          vcf_index=vcf_index,
          out_prefix=out_prefix,
          #GOOGLE_PROJECT=GOOGLE_PROJECT,
          #GCS_OAUTH_TOKEN=GCS_OAUTH_TOKEN
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
        Array[File] vcf
        Array[File] vcf_index 
        String out_prefix
        #String GOOGLE_PROJECT = ""
        #String GCS_OAUTH_TOKEN = ""
    } 

    command <<<
      vcf_input=~{sep=',' vcf}
      #export GCS_REQUESTER_PAYS_PROJECT=~{GOOGLE_PROJECT}
      #export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
      bcftools query -l ~(vcf_input) > ~(out_prefix).txt
    >>>
    
  
    runtime {
        docker:"gcr.io/ucsd-medicine-cast/bcftools-gcs:latest"
    }

    output {
       File outfile = "${out_prefix}.txt"
    }
}



