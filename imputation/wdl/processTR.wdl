version 1.0

workflow processTR {
    input {
        File vcf
        File vcf_index
        String out_prefix

    }

    call extract_TR {
        input:
            vcf=vcf,
            vcf_index=vcf_index,      
            out_prefix=out_prefix

    }


    output {
        File outfile = extract_TR.outvcf
        File outfile_index = extract_TR.outvcf_index
    }

    meta {
        description: "Extract TR from beagle output "
    }
}


task extract_TR {
    input {
        File vcf
        File vcf_index
        String out_prefix
    }
    command <<<
        bcftools view -i 'ID~"EnsTR"' ~{out_prefix}.vcf.gz -Oz -o ~{out_prefix}_TR.vcf.gz && tabix -p vcf ~{out_prefix}_TR.vcf.gz
 
    >>>

    runtime {
        docker:"gcr.io/ucsd-medicine-cast/bcftools-gcs:latest"
    }

    output {
        File outvcf = "${out_prefix}_TR.vcf.gz"
        File outvcf_index = "${out_prefix}_TR.vcf.gz.tbi"

    }    
}
