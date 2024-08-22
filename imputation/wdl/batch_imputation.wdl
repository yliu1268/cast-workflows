version 1.0

import "beagle.wdl" as imputation_t
import "processTR.wdl" as processTR_t
#import "processSNP.wdl" as processSNP_t
import "merge_TR_batch.wdl" as merge_TR_batch_t
#import "merge_SNP_batch.wdl" as merge_SNP_batch_t


workflow batch_imputation {
        input {
                String vcf 
                String vcf_index 
                File ref_panel
                String out_prefix
                String GOOGLE_PROJECT = ""
                Int? mem 
                Int? window_size 
                String? region
                Boolean beagle_region = false
                Array[File] sample_batches = []
                Array[File] sample_index_batches = []
                Int? disk
                Int? overlap
                File map
                

        }
    ### Call subsetting samples with batches ###
        
        Int num_batches = length(sample_batches)
        scatter (i in range(num_batches)) {
                File batch = sample_batches[i]
                File batch_index = sample_index_batches[i]
                call imputation_t.beagle as run_imputation {
                    input:
                        vcf=batch,
                        vcf_index= batch_index,
                        ref_panel=ref_panel,
                        region=region,
                        GOOGLE_PROJECT=GOOGLE_PROJECT,
                        beagle_region=beagle_region,
                        out_prefix=out_prefix+".BATCH"+i,
                        mem=mem,
                        window_size=window_size,
                        disk=disk,
                        overlap=overlap,
                        map=map
                }
                ## extract TR from batches of beagle output
                call processTR_t.processTR as processTR {
                    input:
                        vcf=run_imputation.outfile,
                        vcf_index=run_imputation.outfile_index,
                        out_prefix=out_prefix+".BATCH"+i
     
                }
                ## extract SNP from batches of beagle output
                #call processSNP_t.processSNP as processSNP {
                #    input:
                #        vcf=run_imputation.outfile,
                #        vcf_index=run_imputation.outfile_index,
                #        out_prefix=out_prefix+".BATCH"+i

                #}
        }

        ## use MergeSTR to merge TR
        call merge_TR_batch_t.merge_batch as merge_TR_batch {
            input:
                vcfs=processTR.outfile,
                vcfs_index=processTR.outfile_index,
                out_prefix=out_prefix,
                disk=disk
        }

        ## use bcftools to merge SNP
        #call merge_SNP_batch_t.merge_SNP_batch as merge_SNP_batch {
        #    input:
        #        vcfs=processSNP.outfile,
        #        vcfs_index=processSNP.outfile_index,
        #        out_prefix=out_prefix

        #}


        output {
            File trvcf = merge_TR_batch.outfile
            File trvcf_index = merge_TR_batch.outfile_index
            #Array[File] outvcf = run_imputation.outfile
            #Array[File] outvcf_index = run_imputation.outfile_index
            #File trvcf_index = merge_TR_batch.outfile_index
            #File snpvcf = merge_SNP_batch.outfile
            #File snpvcf_index = merge_SNP_batch.outfile_index
        }

        meta {
            description: "This workflow run imputation on batches of sample, extract TRs and merge across a single chromosome with default parameters "
    }
                   
}


#task annotaTR
#    input {
#        File vcf 
#        File vcf_index
##        File ref_panel
#        File ref_index
#        String out
#        String outtype
#    }
#    
#    command <<<
#        annotaTR --vcf chr11_2batch_merged_noheader.vcf.gz \
#                --ref-panel ensembletr_refpanel_v3_chr11.vcf.gz \
#                --out chr11_2batch_merged_annotated \
#                --vcftype hipstr \
#                --outtype vcf pgen \
#                --dosages beagleap_norm \
#                --ignore-duplicates \
#                --match-refpanel-on locid
#
#    >>>
#
#    runtime {
#        docker:"gcr.io/ucsd-medicine-cast/trtools-6.0.2:latest"
#    }
#
#    output {
#        File outvcf
#        File outvcf_index = 
#        File pgen = 
#        File psam = 
#        File pvar = 
#    }
#    """