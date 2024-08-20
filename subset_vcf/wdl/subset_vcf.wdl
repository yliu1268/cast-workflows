version 1.0

workflow subset_vcf {
    input {
        String multi_sample_vcf
        String GOOGLE_PROJECT = ""
        File sample_groups
        Array[String] regions
    }

    ### Do subsetting for each region ###
    Int num_regions = length(regions)
    scatter (i in range(num_regions)) {
        String region = regions[i]
        call subset_region_batches {
            input:
                multi_sample_vcf=multi_sample_vcf,
                GOOGLE_PROJECT=GOOGLE_PROJECT,
                sample_groups=sample_groups,
                region=region
        }
    }

    ### Output files ####
    output {
        Array[Array[File]] vcf_output_array = subset_region_batches.vcf_output_array
        Array[Array[File]] vcf_index_array = subset_region_batches.vcf_index_array
    }
}

task subset_region_batches {
    input {
        String multi_sample_vcf
        String GOOGLE_PROJECT = ""
        File sample_groups
        String region
    }

    command <<<
        set -e # fail if anything doesn't succeed

        # Set up credentials
        export GCS_REQUESTER_PAYS_PROJECT=~{GOOGLE_PROJECT}
        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

        # Modify output prefixes to include region info
        cat ~{sample_groups} | awk -v"region=~{region}" '{print $0 "-" region}' | sed 's/:/-/' > renamed_sample_groups.txt

        # Get subsets
        bcftools plugin split ~{multi_sample_vcf} -r ~{region} -G renamed_sample_groups.txt -Oz -o .

        # Index
        for f in *.vcf.gz
        do
            tabix -p vcf $f
        done
    >>>

    runtime {
        docker: "gcr.io/ucsd-medicine-cast/bcftools-gcs:latest"
    }

    output {
        Array[File] vcf_output_array = glob("*.vcf.gz")
        Array[File] vcf_index_array = glob("*.vcf.gz.tbi")
    }
}