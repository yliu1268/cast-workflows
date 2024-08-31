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
        Int total = length(batch_files)
    }

    command <<<
        set -e

        # Set up credentials
        export GCS_REQUESTER_PAYS_PROJECT=~{GOOGLE_PROJECT}
        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

        # Fix regions (comment and revert after we change subset jobs to fix this)
        newlist=""
        VCFARRAY=(~{sep=" " batch_files}) # Process one vcf file at a time.
        for (( c = 0; c < ~{total}; c++ )) # bash array are 0-indexed ;)
        do
            vcf=${VCFARRAY[$c]}
            echo "Processing $vcf..."
            export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
            region=$(basename $vcf | cut -d'-' -f 2-4 | sed 's/-/:/' | sed 's/.vcf.gz//')
            outf=$(echo $vcf | sed 's/.vcf.gz/.fixed.vcf.gz/')
            bcftools view $vcf -r $region --regions-overlap 2 -Oz -o ${outf}
            tabix -p vcf ${outf}
            newlist="$newlist $outf"
        done        
        echo $newlist

        # Concatenate and index
        bcftools concat ${newlist} -Oz -o ~{outprefix}_~{batch_prefix}.vcf.gz
        tabix -p vcf ~{outprefix}_~{batch_prefix}.vcf.gz
    >>>

    runtime {
        docker: "gcr.io/ucsd-medicine-cast/bcftools-gcs:latest"
    }

    output {
        File vcf_output = "~{outprefix}_~{batch_prefix}.vcf.gz"
        File vcf_indices = "~{outprefix}_~{batch_prefix}.vcf.gz.tbi"
    }
}