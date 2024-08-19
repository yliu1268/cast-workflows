version 1.0

workflow run_subset{
    input {
        String vcf 
        String vcf_index 
        String out_prefix
        String GOOGLE_PROJECT = ""
        Array[File] samples = []
        File ? sample 
    }

        ### Call subsetting samples with batches ###
    Int num_batches = length(samples)
    
    scatter (i in range(num_batches)) {
            File sample_batch = samples[i]  

            call subset_vcf {
                input:
                    vcf=vcf,
                    vcf_index=vcf_index,
                    sample=sample_batch,
                    GOOGLE_PROJECT=GOOGLE_PROJECT,
                    out_prefix=out_prefix+".BATCH"+i
            }       
    }
    output {
            Array[File] outfile = subset_vcf.outvcf 
            Array[File] outfile_index = subset_vcf.outvcf_index
    }
    meta {
        description: "Run Beagle on a subset of samples on a single chromesome with default parameters"
    }
}

task subset_vcf {
    input {
        String vcf
        String vcf_index
        File sample 
        String GOOGLE_PROJECT = ""
        String out_prefix=out_prefix
     
    }

    command <<<
        set -e 
        export GCS_REQUESTER_PAYS_PROJECT=~{GOOGLE_PROJECT}
        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

        bcftools view -S ~{sample} -I ~{vcf} -Oz -o ~{out_prefix}.vcf.gz && tabix -p vcf ~{out_prefix}.vcf.gz
          
    >>>

    runtime {
        docker:"gcr.io/ucsd-medicine-cast/bcftools-gcs:latest"
        maxRetries: 3
        preemptible: 3
    }

    output {
        File outvcf = "${out_prefix}.vcf.gz"
        File outvcf_index = "${out_prefix}.vcf.gz.tbi"
        
    }    
}