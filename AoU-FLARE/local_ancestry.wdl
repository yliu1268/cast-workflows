version 1.0

workflow LocalAncestryWorkflow {
  ##############################
  # Input definitions
  ##############################
  File input_vcf        # e.g. "gs://.../acaf_threshold_unrelated.chr2.biallelic.phased.vcf.gz"
  File input_vcf_index  # e.g. "gs://.../acaf_threshold_unrelated.chr2.biallelic.phased.vcf.gz.tbi"
  Int chunk_size        # e.g. 12000
  String chromosome     # e.g. "2" or "X"
  String output_prefix  # e.g. "flare_output.acaf-unrelated-all.chr2"
  String output_dir     # e.g. "flare-output"

  File ref_vcf          # e.g. "gs://.../recalibrated_variants.1kgp.chr2.biallelic.phased.vcf.gz"
  File ref_vcf_index    # e.g. "gs://.../recalibrated_variants.1kgp.chr2.biallelic.phased.vcf.gz.tbi"
  File ref_panel        # e.g. "gs://.../ref.panel"
  File map_file         # e.g. "gs://.../plink.chr2.GRCh38.map"
  String flare_jar      # e.g. "gs://.../flare.jar" or local "flare.jar"
  String java_opts      # e.g. "-Xmx512g"
  Boolean probs         # e.g. true (to get posterior probabilities)
  Boolean do_em         # e.g. true (to run EM)

  ##############################
  # Step 1: Chunk the VCF
  ##############################
  call ChunkVCF {
    input:
      in_vcf      = input_vcf,
      chunkSize   = chunk_size
  }

  ##############################
  # Step 2: Run FLARE on each chunk
  ##############################
  scatter (chunkedVCF in ChunkVCF.subset_vcfs) {
    call RunFlare {
      input:
        chunk_vcf     = chunkedVCF,
        ref_vcf       = ref_vcf,
        ref_panel     = ref_panel,
        map_file      = map_file,
        flare_jar     = flare_jar,
        java_opts     = java_opts,
        probs_flag    = probs,
        em_flag       = do_em,
        out_prefix    = output_prefix + ".chunk_" + chunkedVCF.chunk_id
    }
  }

  ##############################
  # Step 3: Merge chunked outputs
  ##############################
  call IndexMergeOutputs {
    input:
      anc_vcf_files   = RunFlare.anc_vcf_gzs,
      global_anc_gzs  = RunFlare.global_anc_gzs,
      output_prefix   = output_prefix,
      output_dir      = output_dir
  }

  output {
    File merged_local_ancestry = IndexMergeOutputs.merged_local_ancestry
    File merged_global_ancestry = IndexMergeOutputs.merged_global_ancestry
  }
}

#
# Task: ChunkVCF
#
task ChunkVCF {
  input {
    File in_vcf
    Int chunkSize
  }

  command <<<
    set -euo pipefail

    # Count total samples
    TOTAL_SAMPLES=$(bcftools query -l ${in_vcf} | wc -l)

    # Create a sample list
    bcftools query -l ${in_vcf} > all_samples.txt

    # We'll create a small script that, for each chunk, extracts a sub-list
    # Because WDL does not allow direct dynamic array creation of File + chunk index,
    # we simply create a small text with chunk ranges and produce symbolic "outputs".
    echo -n "" > subset_files_manifest.txt

    start=1
    chunk_idx=0
    while [ $start -le ${TOTAL_SAMPLES} ]; do
      end=$((start + ${chunkSize} - 1))
      if [ $end -gt ${TOTAL_SAMPLES} ]; then end=${TOTAL_SAMPLES}; fi

      # chunk sample file
      chunk_sample_file="chunk_${chunk_idx}.samples.txt"
      sed -n "${start},${end}p" all_samples.txt > ${chunk_sample_file}

      # bcftools subset
      chunk_vcf="chunk_${chunk_idx}.vcf.gz"
      bcftools view -S ${chunk_sample_file} ${in_vcf} -Oz -o ${chunk_vcf}

      # Index the chunk
      bcftools index -t ${chunk_vcf}

      # Keep track in a manifest
      echo "chunk_${chunk_idx}.vcf.gz  ${chunk_idx}" >> subset_files_manifest.txt

      start=$((end + 1))
      chunk_idx=$((chunk_idx + 1))
    done
  >>>

  output {
    Array[ChunkedVCF] subset_vcfs = read_tsv("subset_files_manifest.txt") { 
      # Convert the 2 columns in the manifest to a WDL struct
      # columns: [file_name, chunk_id]
      file_name = ~{row[0]}
      chunk_id  = ~{row[1]}
    }
  }

  runtime {
    # Adjust as needed
    cpu: 2
    memory: "8G"
  }

  struct ChunkedVCF {
    String file_name
    String chunk_id
  }
}

#
# Task: RunFlare
#
task RunFlare {
  input {
    ChunkedVCF chunk_vcf
    File ref_vcf
    File ref_panel
    File map_file
    String flare_jar
    String java_opts
    Boolean probs_flag
    Boolean em_flag
    String out_prefix
  }

  command <<<
    set -euo pipefail

    # Download or local reference
    # If needed, you might do something like:
    # gsutil cp ${flare_jar} flare.jar

    # We assume chunk_vcf.file_name is local
    chunkVCF="${chunk_vcf.file_name}"

    # Build flare command
    # Use probs=${probs_flag} and em=${em_flag}
    # The Java memory is set in java_opts (e.g. -Xmx512g).
    # The "out" prefix is out_prefix, producing .anc.vcf.gz and .global.anc.gz
    # Note: We assume "em=${em_flag}" is either "em=true" or "em=false".

    java ${java_opts} -jar ${flare_jar} \
      ref=${ref_vcf} \
      ref-panel=${ref_panel} \
      gt=${chunkVCF} \
      map=${map_file} \
      probs=${probs_flag} \
      em=${em_flag} \
      out=${out_prefix}

  >>>

  output {
    File local_anc_vcf = out_prefix + ".anc.vcf.gz"
    File local_anc_vcf_index = out_prefix + ".anc.vcf.gz.tbi"
    File global_anc_gz = out_prefix + ".global.anc.gz"
  }

  runtime {
    # Adjust for how large your job is
    cpu: 16
    memory: "64G"
  }
}

#
# Task: IndexMergeOutputs
#
task IndexMergeOutputs {
  input {
    Array[File] anc_vcf_files
    Array[File] global_anc_gzs
    String output_prefix
    String output_dir
  }

  command <<<
    set -euo pipefail

    mkdir -p ${output_dir}

    echo "Indexing each local ancestry VCF (if needed) ..."
    for vcf in ${sep=' ' anc_vcf_files}; do
      if [ ! -f "${vcf}.tbi" ]; then
        echo "Indexing $vcf..."
        bcftools index -t ${vcf}
      fi
    done

    # Merge local ancestry VCFs
    echo "Merging local ancestry VCF files..."
    # We'll write a file list for bcftools merge
    ls -1 ${anc_vcf_files[@]} > anc_vcf_files.list
    bcftools merge --file-list anc_vcf_files.list -O z -o ${output_dir}/${output_prefix}.merged.anc.vcf.gz

    # Merge global ancestry files
    echo "Merging global ancestry files..."
    ls -1 ${global_anc_gzs[@]} | sort -V > global_anc_files.list
    zcat $(cat global_anc_files.list) | gzip > ${output_dir}/${output_prefix}.merged.global.anc.gz

    echo "All merging steps are complete."
  >>>

  output {
    File merged_local_ancestry = "${output_dir}/${output_prefix}.merged.anc.vcf.gz"
    File merged_global_ancestry = "${output_dir}/${output_prefix}.merged.global.anc.gz"
  }

  runtime {
    cpu: 2
    memory: "16G"
  }
}
