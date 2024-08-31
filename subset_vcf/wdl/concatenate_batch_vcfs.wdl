version 1.0

workflow concatenate_batch_vcfs {
    input {
        Array[Array[String]] batch_vcf_files
        Array[String] batch_names
        String outprefix
        String GOOGLE_PROJECT = ""
    }

    #### Run separate job for each batch ####
    Int num_batches = length(batch_names)
    scatter (i in range(num_batches)) {
        String batch_prefix = batch_names[i]
        Array[String] batch_files = batch_vcf_files[i]
        call concat_batch {
            input:
                batch_files=batch_files,
                GOOGLE_PROJECT=GOOGLE_PROJECT,
                batch_prefix=batch_prefix,
                outprefix=outprefix
        }
    }

    output {
        Array[File] vcf_outputs = concat_batch.vcf_output
        Array[File] vcf_indices = concat_batch.vcf_indices
    }
}

task concat_batch {
    input {
        Array[String] batch_files
        String GOOGLE_PROJECT = ""
        String batch_prefix
        String outprefix
    }

    command <<<
        set -e

        # Set up credentials
        export GCS_REQUESTER_PAYS_PROJECT=~{GOOGLE_PROJECT}
        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

        # Concatenate and index
        bcftools concat ~{sep=' ' batch_files} -Oz -o ~{outprefix}_~{batch_prefix}.vcf.gz

        echo "Sorting..."
        bcftools sort -Oz -o ~{outprefix}_~{batch_prefix}.sorted.vcf.gz ~{outprefix}_~{batch_prefix}.vcf.gz
        tabix -p vcf ~{outprefix}_~{batch_prefix}.sorted.vcf.gz
    >>>

    runtime {
        docker: "gcr.io/ucsd-medicine-cast/bcftools-gcs:latest"
        disks: "local-disk 25 SSD"
    }

    output {
        File vcf_output = "~{outprefix}_~{batch_prefix}.sorted.vcf.gz"
        File vcf_indices = "~{outprefix}_~{batch_prefix}.sorted.vcf.gz.tbi"
    }
}