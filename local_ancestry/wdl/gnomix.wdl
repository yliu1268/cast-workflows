version 1.0

workflow run_gnomix {
    input {
        File samples
        String multi_sample_vcf
        File model
        String chrom
        String out_prefix
        String GOOGLE_PROJECT = ""
        File chainfile
        File refpanel
        File refpanel_index
        String extra_subset_args
    }

    call subset_vcf {
        input:
            multi_sample_vcf=multi_sample_vcf,
            samples=samples,
            out_prefix=out_prefix,
            GOOGLE_PROJECT=GOOGLE_PROJECT,
            chainfile=chainfile,
            extra_subset_args=extra_subset_args,
            chrom=chrom
    }

    call beagle {
        input:
            vcf=subset_vcf.outvcf,
            vcf_index=subset_vcf.outvcf_index,
            out_prefix=out_prefix,
            refpanel=refpanel,
            refpanel_index=refpanel_index
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

# TODO add back retries after debugging
task subset_vcf {
    input {
        String multi_sample_vcf
        File samples
        String out_prefix
        String GOOGLE_PROJECT = ""
        String hg38_ref = "gs://genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.fasta"
        String hg19_ref = "gs://gcp-public-data--broad-references/Homo_sapiens_assembly19_1000genomes_decoy/Homo_sapiens_assembly19_1000genomes_decoy.fasta"
        File chainfile
        String extra_subset_args = ""
        String chrom
    }

    command <<<        
        # Set up credentials
        export GCS_REQUESTER_PAYS_PROJECT=~{GOOGLE_PROJECT}
        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

        # Extract subset
        bcftools view ~{extra_subset_args} -S ~{samples} -I -m2 -M2 --min-af 0.01 ~{multi_sample_vcf} -Oz -o ~{out_prefix}.vcf.gz
        tabix -p vcf ~{out_prefix}.vcf.gz

        # Liftover to hg19
        bcftools +liftover --no-version -Ou ~{out_prefix}.vcf.gz -- \
              -f ~{hg19_ref} \
              -s ~{hg38_ref} \
              -c ~{chainfile} \
              --reject ~{out_prefix}.reject.bcf \
              --reject-type b \
              --write-src --drop-tags FORMAT/AD | \
              bcftools sort -o ~{out_prefix}_hg19.vcf.gz -Oz
        tabix -p vcf ~{out_prefix}_hg19.vcf.gz

        # Restrict to target chromosome
        bcftools view -r ~{chrom} ~{out_prefix}_hg19.vcf.gz -Oz -o ~{out_prefix}.filtered.vcf.gz
        tabix -p vcf ~{out_prefix}.filtered.vcf.gz
    >>>

    runtime {
        docker: "gcr.io/ucsd-medicine-cast/bcftools-gcs-plugins:latest"
    }

    output {
        File outvcf = "${out_prefix}.filtered.vcf.gz"
        File outvcf_index = "${out_prefix}.filtered.vcf.gz.tbi"
    }    
}

task beagle {
    input {
        File vcf
        File vcf_index
        String out_prefix
        File refpanel
        File refpanel_index
    }

    command <<<
    java -Xmx25g -jar /beagle.jar \
        gt=~{vcf} \
        ref=~{refpanel} \
        impute=false \
        out=~{out_prefix}_phased
    tabix -p vcf ~{out_prefix}_phased.vcf.gz
    >>>

    runtime {
        docker: "gcr.io/ucsd-medicine-cast/beagle:latest"
        memory: "25GB"
    }

    output {
       File outvcf = "${out_prefix}_phased.vcf.gz"
       File outvcf_index = "${out_prefix}_phased.vcf.gz.tbi"
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
        tar -xzvf ~{model}
        echo "debugging"
        ls # TODO remove. for debugging
        echo ${PWD}
        cd /gnomix
        echo ${PWD}
        ls
        python3 gnomix.py ~{vcf} . ~{chrom} False pretrained_gnomix_models/chr~{chrom}/model_chm_~{chrom}.pkl
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