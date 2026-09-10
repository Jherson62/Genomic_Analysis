#!/bin/bash
# Script to download sources to gatk performance
# Germline variant discovery

# Exit when any command fails
# e: Stop on first failure
# u: Stop if missing input
set -euo pipefail

# Conda initialize
# shellcheck disable=SC1091
source "$HOME/miniconda3/etc/profile.d/conda.sh"

# Directories
recursos="$HOME/Documentos/Fertility/Exomas/recursos"
mkdir -p "${recursos}"

# Environment 1: BWA
conda activate NGStools

# Human references fasta file
wget -nc -P "${recursos}" "https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"

# Remove innecesary information from headers
sed -i "/^>/ s/ .*$//" Homo_sapiens_assembly38.fasta

# Indexing references genome to alignment
echo "-------------"
echo "bwa index ..."
echo "-------------"

bwa index "${recursos}/Homo_sapiens_assembly38.fasta"

# Environment 2: Samtools
conda deactivate
conda activate medaka_env

# Indexing references genome to running haplotype caller
echo "------------------"
echo "Samtools index ..."
echo "------------------"
samtools faidx "${recursos}/Homo_sapiens_assembly38.fasta"

# Environment 3: GATK
conda deactivate
conda activate gatk_env

echo "-------------------------------"
echo "Creating Sequence Dictionary..."
echo "-------------------------------"

gatk CreateSequenceDictionary \
    R="${recursos}/Homo_sapiens_assembly38.fasta" \
    O="${recursos}/Homo_sapiens_assembly38.dict"

echo "Finished"