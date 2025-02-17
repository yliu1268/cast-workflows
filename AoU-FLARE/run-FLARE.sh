#!/bin/bash

# Prompt user for chunk size and chromosome number
echo -n "Type in the chunk size you want to use for FLARE: "
read CHUNK_SIZE
echo -n "Type in the chromosome number you want to run FLARE (1-22, X): "
read CHR_NUM

# Set variables
VCF_FILE="acaf_threshold_unrelated.chr${CHR_NUM}.biallelic.phased.vcf.gz"
OUTPUT_DIR="flare-output"

# Create output directory
mkdir -p $OUTPUT_DIR

# Get the total number of samples
TOTAL_SAMPLES=$(bcftools query -l $VCF_FILE | wc -l)

# Function to process each chunk
process_chunk() {
    local start=$1
    local end=$2
    local chunk_id=$3
    local vcf_file=$4
    local output_dir=$5

    tmp_sample_file=$(mktemp)
    bcftools query -l "$vcf_file" | sed -n "${start},${end}p" > "$tmp_sample_file"
    bcftools view -S "$tmp_sample_file" "$vcf_file" -Oz -o "${output_dir}/chunk_${chunk_id}.vcf.gz"
    rm "$tmp_sample_file"
}

export -f process_chunk

# Generate chunk ranges and parallelize processing
seq 1 $CHUNK_SIZE $TOTAL_SAMPLES | \
    parallel -j80 --colsep ' ' '
        start={};
        end=$((start + '$CHUNK_SIZE' - 1));
        if [ $end -gt '$TOTAL_SAMPLES' ]; then end='$TOTAL_SAMPLES'; fi;
        chunk_id=$(((start - 1) / '$CHUNK_SIZE'));
        process_chunk $start $end $chunk_id "'$VCF_FILE'" "'$OUTPUT_DIR'"
    '

echo "Splitting chunks complete."

# Set FLARE variables
JAVA_OPTS="-Xmx512g"
FLARE_JAR="flare.jar"
REF_VCF_FILE="recalibrated_variants.1kgp.chr${CHR_NUM}.biallelic.phased.vcf.gz"
REF_PANEL="ref.panel"
MAP_FILE="plink.chr${CHR_NUM}.GRCh38.map"
OUTPUT_PREFIX="flare-output/flare_output.acaf-unrelated-all.chr${CHR_NUM}"

# Ensure output directory exists
mkdir -p $OUTPUT_DIR

# Function to run FLARE on each chunk
process_chunk_with_flare() {
    find ${OUTPUT_DIR}/chunk_*.vcf.gz | sort -V | while read chunk_file; do
        chunk_id=$(basename $chunk_file | sed 's/chunk_\([0-9]*\).vcf.gz/\1/')
        
        time java $JAVA_OPTS -jar $FLARE_JAR \
            ref=$REF_VCF_FILE \
            ref-panel=$REF_PANEL \
            gt=$chunk_file \
            map=$MAP_FILE \
            probs=true \
            em=true \
            out=${OUTPUT_PREFIX}.chunk_${chunk_id}
    done

    echo "FLARE processing complete."
}

# Function to index and merge VCF files
index_and_merge_vcf_files() {
    echo "Indexing VCF files..."
    for vcf in ${OUTPUT_PREFIX}.chunk_*.anc.vcf.gz; do
        if [ ! -f "$vcf.tbi" ]; then
            echo "Indexing $vcf..."
            bcftools index -t $vcf
        fi
    done

    echo "Merging local ancestry VCF files..."
    bcftools merge --file-list <(ls ${OUTPUT_PREFIX}.chunk_*.anc.vcf.gz | sort -V) -O z -o ${OUTPUT_PREFIX}.merged.anc.vcf.gz

    echo "Merging global ancestry files..."
    zcat $(ls ${OUTPUT_PREFIX}.chunk_*.global.anc.gz | sort -V) | gzip > ${OUTPUT_PREFIX}.merged.global.anc.gz

    echo "Processing and merging complete."
}

process_chunk_with_flare
index_and_merge_vcf_files
