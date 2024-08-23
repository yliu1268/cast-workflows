version 1.0

workflow run_gnomix {
    input {
        File vcf
        File model
        String chrom
        String out_prefix
        String GOOGLE_PROJECT = ""
        File chainfile
        File refpanel
        File refpanel_index
    }

    call liftover_vcf {
        input:
            vcf=vcf,
            out_prefix=out_prefix,
            GOOGLE_PROJECT=GOOGLE_PROJECT,
            chainfile=chainfile,
            chrom=chrom
    }

    call beagle {
        input:
            vcf=liftover_vcf.outvcf,
            vcf_index=liftover_vcf.outvcf_index,
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

task liftover_vcf {
    input {
        File vcf
        String out_prefix
        String GOOGLE_PROJECT = ""
        File hg38_ref = "gs://genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.fasta"
        File hg19_ref = "gs://gcp-public-data--broad-references/Homo_sapiens_assembly19_1000genomes_decoy/Homo_sapiens_assembly19_1000genomes_decoy.fasta"
        File chainfile
        String chrom
    }

    command <<<
        set -e # fail if anything doesn't succeed

        # Liftover to hg19
        # First refresh credentials
#        export GCS_REQUESTER_PAYS_PROJECT=~{GOOGLE_PROJECT}
#        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
        echo "Liftover"
        bcftools +liftover --no-version -Ou ~{vcf} -- \
              -f ~{hg19_ref} \
              -s ~{hg38_ref} \
              -c ~{chainfile} \
              --reject ~{out_prefix}.reject.bcf \
              --reject-type b \
              --write-src --drop-tags FORMAT/AD > ~{out_prefix}_hg19.bcf
        echo "sort"
        bcftools sort -o ~{out_prefix}_hg19.vcf.gz -Oz ~{out_prefix}_hg19.bcf
        tabix -p vcf ~{out_prefix}_hg19.vcf.gz

        # Restrict to target chromosome
        echo "restrict chromosome"
        bcftools view -r ~{chrom} ~{out_prefix}_hg19.vcf.gz -Oz -o ~{out_prefix}.filtered.vcf.gz
        tabix -p vcf ~{out_prefix}.filtered.vcf.gz
    >>>

    runtime {
        docker: "gcr.io/ucsd-medicine-cast/bcftools-gcs-plugins:latest"
        memory: "25GB"
        disks: "local-disk 50 SSD"
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
        impute=true \
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
        cd /gnomix
        tar -xzvf ~{model}
        python3 gnomix.py ~{vcf} . ~{chrom} False pretrained_gnomix_models/chr~{chrom}/model_chm_~{chrom}.pkl
        cp query_results.msp /cromwell_root/~{out_prefix}.msp
        cp query_results.fb /cromwell_root/~{out_prefix}.fb
    >>>

    runtime {
        docker: "gcr.io/ucsd-medicine-cast/gnomix:latest"
        memory: "100GB"
    }

    output {
       File msp_outfile = "~{out_prefix}.msp"
       File fb_outfile = "~{out_prefix}.fb"
    }
}