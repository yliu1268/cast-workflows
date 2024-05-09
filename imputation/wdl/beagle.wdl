version 1.0

workflow beagle {
    input {
        File vcf 
        File vcf_index 
        File ref_panel
        File ref_panel_index
        String out_prefix
        String GOOGLE_PROJECT = ""
        String GCS_OAUTH_TOKEN = ""
        Int? mem = 32
        Int? window_size = 20
        File? excludesamples = excludesamples
    }

    call beagle {
        input : 
          vcf=vcf, 
          vcf_index=vcf_index,
          ref_panel=ref_panel, 
          ref_panel_index=ref_panel_index,
          out_prefix=out_prefix,
          GOOGLE_PROJECT=GOOGLE_PROJECT,
          GCS_OAUTH_TOKEN=GCS_OAUTH_TOKEN,
          mem=mem,
          window_size=window_size,
          excludesamples=excludesamples
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
      description: "Run Beagle on a single chromesome with default parameters"
    }
}

task beagle {
    input {
        File vcf 
        File vcf_index 
        File ref_panel
        File ref_panel_index
        String out_prefix
        String GOOGLE_PROJECT = ""
        String GCS_OAUTH_TOKEN = ""
        Int? mem = 32
        Int? window_size = 20
        File? excludesamples = excludesamples
    } 

    command <<<
        #set -e
    # We need at least 1 GB of available memory outside of the Java heap in order to execute native code, thus, limit
    # Java's memory by the total memory minus 1 GB. We need to compute the total memory as it might differ from
    # memory_size_gb because of Cromwell's retry with more memory feature.
    # Note: In the future this should be done using Cromwell's ${MEM_SIZE} and ${MEM_UNIT} environment variables,
    #       which do not rely on the output format of the `free` command.
    #    available_memory_mb=$(free -m | awk '/^Mem/ {print $2}')
    #   let java_memory_size_mb=available_memory_mb-3072
    #    echo Total available memory: ${available_memory_mb} MB >&2
    #    echo Memory reserved for Java: ${java_memory_size_mb} MB >&2
    #    #export GCS_REQUESTER_PAYS_PROJECT=~{GOOGLE_PROJECT}
    #   #export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

    #    java -Xmx${java_memory_size_mb}m -Xms${java_memory_size_mb}m -jar /beagle.jar  \
    #            gt=~{vcf} \
    #            ref=~{ref_panel} \
    #            out=~{out_prefix}
        export GCS_REQUESTER_PAYS_PROJECT=~{GOOGLE_PROJECT}
        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
        echo "Java max heap size"
        java -XX:+PrintFlagsFinal -version | grep HeapSize
        echo "system memory"
        grep MemTotal /proc/meminfo | awk '{print $2}'
  
        java -Xmx~{mem}g -jar /beagle.jar \
            gt=~{vcf} \
            ref=~{ref_panel} \
            window=~{window_size} \
            excludesamples=~{excludesamples} \
            out=~{out_prefix}





    >>>
    
    #file upto 300mb use mem=25
    runtime {
        docker:"gcr.io/ucsd-medicine-cast/beagle:latest"
	    memory: mem + "GB"
    }

    output {
       File outfile = "${out_prefix}.vcf.gz"
    }
}

task sort_index_beagle {
    input {
      File vcf
    }

    String basename = basename(vcf, ".vcf.gz")

    command <<<
        zcat ~{vcf} | vcf-sort | bgzip -c > ~{basename}.sorted.vcf.gz && tabix -p vcf ~{basename}.sorted.vcf.gz
    >>>

    runtime {
        docker:"gcr.io/ucsd-medicine-cast/vcfutils:latest"
    }

    output {
    File outvcf = "${basename}.sorted.vcf.gz"
    File outvcf_index = "${basename}.sorted.vcf.gz.tbi"
  }
}
