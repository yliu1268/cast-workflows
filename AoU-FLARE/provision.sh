#!/usr/bin/env bash
set -euo pipefail

echo "Installing conda dependencies and bcftools..."

# Create conda environment (if not already existing).
# Adjust conda environment name or Java version as desired.
conda create -y -n java21_env openjdk=21 bcftools
conda activate java21_env

# Optional: install GNU Parallel in the environment as well
conda install -y parallel

# Download FLARE jar
# (This example uses the official browning-lab release)
wget -O flare.jar https://faculty.washington.edu/browning/flare.jar

echo "Provisioning complete."
