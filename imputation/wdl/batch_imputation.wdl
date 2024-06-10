version 1.0

import "imputation.wdl" as imputation_single

workflow batch_imputation {
	input {
            String vcf 
            String vcf_index 
            File ref_panel
            String out_prefix
            String GOOGLE_PROJECT = ""
            Int? mem 
            Int? window_size 
            File sample_file 
	        String? region
            Boolean subset_region = false
            Boolean beagle_region = false
            Array[String] samples = []

	}

	### Call subsetting samples with batches ###
    scatter (i in range(length(samples))) {
        String sample = samples [i]
        call imputation_single.subset_vcf as imputation_batch_sample {
            input:
                sample_file=sample,
                region=region,
                vcf=vcf,
                vcf_index=vcf_index,
                GOOGLE_PROJECT=GOOGLE_PROJECT,
                subset_region=subset_region,
                out_prefix=out_prefix,
        }
    

        call imputation_single.index_vcf as index_vcf{
            input:
                vcf=imputation_batch_sample.outfile
        }

        call imputation_single.beagle as beagle{
            input: 
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

        call imputation_single.sort_index_beagle as sort_index_beagle {
            input:
                vcf=beagle.outfile
        }
    }
    call merge_vcf {
        input:
            vcf=sort_index_beagle.outvcf,
            out_prefix=out_prefix
    }

        
    output {
        File outfile = merge_vcf.outvcf 
        File outfile_index = merge_vcf.outvcf_index
    }
    meta {
    description: "Run Beagle on batches of samples on a single chromesome with default parameters"
    }
    
}


task merge_vcf {
	input {
        Array[File] vcf
        String out_prefix
	}

	command <<<
		bcftools merge ${sep=' ' vcf} -O z -o ${out_prefix}.merged.vcf.gz
	>>>

	runtime {
        docker:"gcr.io/ucsd-medicine-cast/bcftools-gcs:latest"
    }

	output {
		File outvcf = "${out_prefix}.merged.vcf.gz"
		File outvcf_index = "${out_prefix}.merged.vcf.gz.tbi"
	}
}