version 1.0

workflow run_imputation {
    input {
        String vcf 
        String vcf_index 
        File ref_panel
        String out_prefix
        String GOOGLE_PROJECT = ""
        Int? mem 
        Int? window_size
        String? region
        File sample
        Boolean subset_region = false
        Boolean beagle_region = false
    }

       
    call subset_vcf {
        input:
            vcf=vcf,
            vcf_index=vcf_index,
            sample=sample,
            region=region,
            GOOGLE_PROJECT=GOOGLE_PROJECT,
            subset_region=subset_region,
            out_prefix=out_prefix
            }

    call index_vcf {
        input:
            vcf=subset_vcf.outfile
    }
    
    call beagle {
        input : 
          vcf=index_vcf.outvcf, 
          vcf_index=index_vcf.outvcf_index,
          ref_panel=ref_panel, 
          out_prefix=out_prefix,
          GOOGLE_PROJECT=GOOGLE_PROJECT,
          mem=mem,
          window_size=window_size,
          beagle_region=beagle_region,
          region=region
    }
    call sort_index_beagle {
        input :
            vcf=beagle.outfile
    }
    output {
        File outfile = sort_index_beagle.outvcf 
        File outfile_index = sort_index_beagle.outvcf_index
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
	    String? region
        String GOOGLE_PROJECT = ""
        String out_prefix=out_prefix
        Boolean subset_region = false   
    }

    command <<<
        export GCS_REQUESTER_PAYS_PROJECT=~{GOOGLE_PROJECT}
        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

        # Subsetting region for each chromesome

        if [[ "~{subset_region}" == false ]] ; then
            #bcftools view -S ~{sample} ~{vcf} > ~{out_prefix}.vcf
        bash /usr/bin/subset_vcf.sh ~{vcf} ~{sample} ~{out_prefix} ~{GOOGLE_PROJECT} || exit 1

        else 
            bcftools view -r ~{region} -S ~{sample} ~{vcf} > ~{out_prefix}.vcf
        fi 
    >>>

    runtime {
        docker:"gcr.io/ucsd-medicine-cast/bcftools-gcs:latest"
        maxRetries: 3
        preemptible: 3
    }

    output {
        File outfile = "${out_prefix}.vcf"
    }    
}


task index_vcf {
    input {
      File vcf
    }

    String basename = basename(vcf, ".vcf")

    command <<<
        bgzip -c ~{vcf}> ~{basename}.vcf.gz && tabix -p vcf ~{basename}.vcf.gz
    >>>

    runtime {
        docker:"gcr.io/ucsd-medicine-cast/vcfutils:latest"

    }

    output {
        File outvcf = "${basename}.vcf.gz"
        File outvcf_index = "${basename}.vcf.gz.tbi"
    }
}

task beagle {
    input {
        File vcf 
        File vcf_index 
        File ref_panel
        String out_prefix
        String GOOGLE_PROJECT = ""
        Int? mem 
        Int? window_size
        Boolean beagle_region = false
        String? region
    } 

    command <<<

        export GCS_REQUESTER_PAYS_PROJECT=~{GOOGLE_PROJECT}
        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
        if [[ "~{beagle_region}" == false ]] ; then
        java -Xmx~{mem}g -jar /beagle.jar \
            gt=~{vcf} \
            ref=~{ref_panel} \
            window=~{window_size} \
            ap=true \
            out=~{out_prefix}_output
        else
            java -Xmx~{mem}g -jar /beagle.jar \
            gt=~{vcf} \
            ref=~{ref_panel} \
            window=~{window_size} \
            chrom=~{region} \
            ap=true \
            out=~{out_prefix}_output
        fi
    >>>
    
    #file upto 300mb use mem=25
    runtime {
        docker:"gcr.io/ucsd-medicine-cast/beagle:latest"
	    memory: mem + "GB"
    }

    output {
       File outfile = "${out_prefix}_output.vcf.gz"
    }
}

task sort_index_beagle {
    input {
      File vcf
    }

    String basename = basename(vcf, ".vcf.gz")

    command <<<
        zcat ~{vcf} | vcf-sort | bgzip -c > ~{basename}.sorted.vcf.gz && tabix -p vcf ~{basename}.sorted.vcf.gz
        echo "Number of TRs in the genotyped file"
        bcftools view -i 'ID="."' ~{basename}.sorted.vcf.gz | grep -v "^#" | wc -l
    >>>

    runtime {
        docker:"gcr.io/ucsd-medicine-cast/vcfutils:latest"
    }

    output {
    File outvcf = "${basename}.sorted.vcf.gz"
    File outvcf_index = "${basename}.sorted.vcf.gz.tbi"
  }
}
