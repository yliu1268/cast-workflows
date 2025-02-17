#!/bin/bash

# Download the phased vcf files from AoU workspace bucket (as input file for FLARE)
# gsutil -u $GOOGLE_PROJECT cp gs://.../phased-vcf/acaf_threshold_unrelated.chr2.biallelic.phased.vcf.gz .
# gsutil -u $GOOGLE_PROJECT cp gs://.../phased-vcf/acaf_threshold_unrelated.chr2.biallelic.phased.vcf.gz.tbi .

##########################################################################
# Download GNU Parallel to parallel the splitting job of bcftools
# wget -qO parallel.tar.bz2 https://ftp.gnu.org/gnu/parallel/parallel-latest.tar.bz2

# Extract the tarball
# tar -xf parallel.tar.bz2

# cd parallel-*

# Use GNU Parallel without installing
# export PATH=$(pwd)/src:$PATH

# cd ../

##########################################################################
# Set variables 
VCF_FILE="acaf_threshold_unrelated.chr2.biallelic.phased.vcf.gz"
CHUNK_SIZE=12000 #Sample size for each chunk
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

    # Create a temporary file for the sample IDs of this chunk
    tmp_sample_file=$(mktemp)
    bcftools query -l "$vcf_file" | sed -n "${start},${end}p" > "$tmp_sample_file"

    # Create separate VCF files for each chunk of sample IDs
    bcftools view -S "$tmp_sample_file" "$vcf_file" -Oz -o "${output_dir}/chunk_${chunk_id}.vcf.gz"

    # Clean up temporary files
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
        process_chunk $start $end $chunk_id '"$VCF_FILE"' '"$OUTPUT_DIR"'
    '

echo "Splitting chunks complete."


##########################################################################
# Make sure to have java21 version to run FLARE, otherwise, create conda environment
# cd ../
# export CONDA_PKGS_DIRS=~/conda_pkgs/
# export CONDA_ENVS_DIRS=~/conda_envs/
# conda create -n java21_env openjdk=21 -y

# cd beagle-flare-chr2
# conda activate java21_env

##########################################################################
# Download reference: we phased 2504 unrelated-1kgp samples with BEAGLE
# gsutil -u $GOOGLE_PROJECT cp gs://.../unrelated-1kgp/recalibrated_variants.1kgp.chr2.biallelic.phased.vcf.gz.tbi .
# gsutil -u $GOOGLE_PROJECT cp gs://.../unrelated-1kgp/recalibrated_variants.1kgp.chr2.biallelic.phased.vcf.gz .


# Set variables
JAVA_OPTS="-Xmx512g"
FLARE_JAR="flare.jar"
REF_VCF_FILE="recalibrated_variants.1kgp.chr2.biallelic.phased.vcf.gz"
REF_PANEL="ref.panel"
MAP_FILE="plink.chr2.GRCh38.map"
OUTPUT_PREFIX="flare-output/flare_output.acaf-unrelated-all.chr2"
OUTPUT_DIR="flare-output"

# Ensure output directory exists
mkdir -p $OUTPUT_DIR

# Function to run FLARE on each chunk
process_chunk_with_flare() {
    # Process each chunk in numerical order
    find ${OUTPUT_DIR}/chunk_*.vcf.gz | sort -V | while read chunk_file; do
        chunk_id=$(basename $chunk_file | sed 's/chunk_\([0-9]*\).vcf.gz/\1/')
        
        # Run FLARE for the chunk
        time java $JAVA_OPTS -jar $FLARE_JAR \
            ref=$REF_VCF_FILE \
            ref-panel=$REF_PANEL \
            # Optional: use a sorted-ref-panel to ensure that ancestry names are always in alphabetical order in the output
            ### sort -k2,2 $REF_PANEL > sorted_ref_panel.txt
            gt=$chunk_file \
            map=$MAP_FILE \
            probs=true \
            em=true \
            # Optional: provide the model file here to specify the ancestry order for FLARE output
            ### model=$MODEL_FILE \ 
            # Otherwaise FLARE will use the reference panel names in the ref-panel file, in the order each name first appears in the ref-panel file, as the ancestry names.
            # In our case: EUR EAS AMR CSA AFR
            out=${OUTPUT_PREFIX}.chunk_${chunk_id}
    done

    echo "FLARE processing complete." 
}

# Function to index and merge VCF files
index_and_merge_vcf_files() {
    # Step 1: Index all VCF files if not already indexed
    echo "Indexing VCF files..."
    for vcf in ${OUTPUT_PREFIX}.chunk_*.anc.vcf.gz; do
        if [ ! -f "$vcf.tbi" ]; then
            echo "Indexing $vcf..."
            bcftools index -t $vcf
        fi
    done

    # Step 2: Merge local ancestry VCF files
    echo "Merging local ancestry VCF files..."
    bcftools merge --file-list <(ls ${OUTPUT_PREFIX}.chunk_*.anc.vcf.gz | sort -V) -O z -o ${OUTPUT_PREFIX}.merged.anc.vcf.gz

    # Step 3: Concatenate global ancestry files using zcat
    echo "Merging global ancestry files..."
    zcat $(ls ${OUTPUT_PREFIX}.chunk_*.global.anc.gz | sort -V) | gzip > ${OUTPUT_PREFIX}.merged.global.anc.gz

    echo "Processing and merging complete."
}

# Run the functions in order
process_chunk_with_flare
index_and_merge_vcf_files
