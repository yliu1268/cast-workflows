version 1.0

workflow run_gnomix {
    input {
        File samples
        String multi_sample_vcf
        File model
        String chrom
        String out_prefix
        String GOOGLE_PROJECT = ""
    }

    call subset_vcf {
        input:
            multi_sample_vcf=multi_sample_vcf,
            samples=samples,
            out_prefix=out_prefix,
            GOOGLE_PROJECT=GOOGLE_PROJECT
    }

    call beagle {
        input:
            vcf=subset_vcf.outvcf,
            vcf_index=subset_vcf.outvcf_index,
            out_prefix=out_prefix
    }

    call gnomix {
        input:
            vcf=beagle.outvcf,
            vcf_index=beagle.outvcf_index,
            model=model,
            chrom=chrom,
            out_prefix=out_prefix
    }

    output {
        File msp_outfile = gnomix.msp_outfile
        File fb_outfile = gnomix.fb_outfile
    }

    meta {
        description: "Run beagle and gnomix on a specified subset of samples on a single chromosome"
    }
}

task subset_vcf {
    input {
        String multi_sample_vcf
        File samples
        String out_prefix
        String GOOGLE_PROJECT = ""
    }

    command <<<
        # Set up credentials
        export GCS_REQUESTER_PAYS_PROJECT=~{GOOGLE_PROJECT}
        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
        bcftools view -S ~{samples} -I ~{multi_sample_vcf} -Oz -o ~{out_prefix}.vcf.gz
        tabix -p vcf ~{out_prefix}.vcf.gz
    >>>

    runtime {
        docker: "gcr.io/ucsd-medicine-cast/bcftools-gcs:latest"
        maxRetries: 3
        preemptible: 3
    }

    output {
        File outvcf = "${out_prefix}.vcf.gz"
        File outvcf_index = "${out_prefix}.vcf.gz.tbi"
    }    
}

task beagle {
    input {
        File vcf
        File vcf_index
        String out_prefix
    }

    command <<<
    java -Xmx25g -jar /beagle.jar \
        gt=~{vcf} \
        out=~{out_prefix}_phased
    tabix -p vcf ~{out_prefix}_phased.vcf.gz
    >>>

    runtime {
        docker: "gcr.io/ucsd-medicine-cast/beagle:latest"
        memory: "25GB"
    }

    output {
       File outvcf = "${out_prefix}_output.vcf.gz"
       File outvcf_index = "${out_prefix}_output.vcf.gz.tbi"
    }
}

task gnomix {
    input {
        File vcf
        File vcf_index
        File model
        String chrom
        String out_prefix
    }

    command <<<
        python3 gnomix.py ~{vcf} . ~{chrom} False ~{model}
        cp query_results.msp ~{out_prefix}.msp
        cp query_results.fb ~{out_prefix}.fb
    >>>

    runtime {
        docker: "gcr.io/ucsd-medicine-cast/gnomix:latest"
    }

    output {
       File msp_outfile = "~{out_prefix}.msp"
       File fb_outfile = "~{out_prefix}.fb"
    }
}