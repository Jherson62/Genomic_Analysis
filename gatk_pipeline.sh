#!/bin/bash
# exit when any command fails
set -eou pipefail

# Conda initialize
# shellcheck disable=SC1091
source "$HOME/miniconda3/etc/profile.d/conda.sh"


# Directories
sample="EX2601"
# quality="$HOME/Documentos/Fertility/Exomas/${sample}/quality"
resources="$HOME/Documentos/Fertility/Exomas/resources/" # Resources
data="$HOME/Documentos/Fertility/Exomas/${sample}/data" # Paciente data
aligned="$HOME/Documentos/Fertility/Exomas/${sample}/aligned"
results="$HOME/Documentos/Fertility/Exomas/${sample}/results"


# Enviroment 1: BWA
conda activate NGStools

echo "---------------------------------------"
echo " Map to reference using BWA-MEM"
echo "---------------------------------------"

bwa mem -t 8 -R "@RG\tID:${sample}\tPL:ILLUMINA\tSM:${sample}" \
    "${resources}"/Homo_sapiens_assembly38.fasta \
    "${data}"/${sample}_1.fastq.gz \
    "${data}"/${sample}_2.fastq.gz > "${aligned}"/${sample}.paired.sam

# Enviroment 2: GATK
conda activate gatk_env

echo "------------------------------------"
echo "Mark Duplicates and add Tags ..."
echo "------------------------------------"

gatk --java-options "-Xms8G -Xmx8G -XX:+UseStringDeduplication" MarkDuplicatesSpark \
    -I "${aligned}"/${sample}_paired.sam \
    -O "${aligned}"/${sample}_sort_dedup.bam

gatk SetNmMdAndUqTags \
    -I "${aligned}"/${sample}_sort_dedup.bam \
    -O "${aligned}"/${sample}_sort_dedup_tag.bam \
    -R "${resources}"/Homo_sapiens_assembly38.fasta

echo "------------------------------------"
echo "Base Quality Score Recalibration ..."
echo "------------------------------------"

# make a math model
 gatk --java-options "-Xms8G -Xmx8G -XX:+UseStringDeduplication" BaseRecalibrator \
    -I "${aligned}"/${sample}_sort_dedup_tag.bam \
    -R "${resources}"/Homo_sapiens_assembly38.fasta \
    --known-sites "${resources}"/Homo_sapiens_assembly38.dbsnp138.vcf \
    -O "${data}"/recal_data.table

# apply the math model
gatk --java-options "-Xms8G -Xmx8G -XX:+UseStringDeduplication" ApplyBQSR \
    -I "${aligned}"/${sample}_sort_dedup_tag.bam \
    -R "${resources}"/Homo_sapiens_assembly38.fasta \
    --bqsr-recal-file "${data}"/recal_data.table \
    -O "${aligned}"/${sample}_sort_dedup_tag_bqsr.bam

echo "-------------------------------------------"
echo "Collect Alignment & Insert Size Metrics ..."
echo "-------------------------------------------"

gatk CollectAlignmentSummaryMetrics \
    R="${resources}"/Homo_sapiens_assembly38.fasta \
    I="${aligned}"/${sample}_sort_dedup_tag_bqsr.bam \
    O="${aligned}"/alignment_metrics.txt


gatk CollectInsertSizeMetrics \
    INPUT="${aligned}"/${sample}_sort_dedup_tag_bqsr.bam \
    OUTPUT="${aligned}"/insert_size_metrics.txt \
    HISTOGRAM_FILE="${aligned}"/insert_size_histogram.pdf


echo "--------------------"
echo "HaplotypeCaller ..."
echo "--------------------"

# Call Variants . . .
gatk --java-options "-Xms8G -Xmx8G -XX:+UseStringDeduplication" HaplotypeCaller \
    -R "${resources}"/Homo_sapiens_assembly38.fasta \
    -I "${aligned}"/${sample}_sort_dedup_tag_bqsr.bam \
    -O "${results}"/raw_variants.vcf \
    -L "${resources}"/cromosomas_principales.list


echo "-----------------"
echo "Separate SNPs ..."
echo "-----------------"

gatk SelectVariants \
    -R "${resources}"/Homo_sapiens_assembly38.fasta \
    -V "${results}"/raw_variants.vcf \
    --select-type SNP \
    -O "${results}"/raw_snps.vcf

# Filter SNPs
gatk VariantFiltration \
	-R "${resources}"/Homo_sapiens_assembly38.fasta \
	-V "${results}"/raw_snps.vcf \
	-O "${results}"/filtered_snps.vcf \
	-filter-name "QD_filter" -filter "QD < 2.0" \
	-filter-name "FS_filter" -filter "FS > 60.0" \
	-filter-name "MQ_filter" -filter "MQ < 40.0" \
	-filter-name "SOR_filter" -filter "SOR > 4.0" \
	-filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
	-filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0" \
	-genotype-filter-expression "DP < 10" \
	-genotype-filter-name "DP_filter" \
	-genotype-filter-expression "GQ < 10" \
	-genotype-filter-name "GQ_filter"

# Select Variants - In case you need to annotate with gatk database
# gatk SelectVariants \
# 	--exclude-filtered \
# 	-V "${results}"/filtered_snps.vcf \
# 	-O "${results}"/analysis_ready_snps.vcf

echo "-----------------"
echo "Separate INDELs ..."
echo "-----------------"

gatk SelectVariants \
    -R "${resources}"/Homo_sapiens_assembly38.fasta \
    -V "${results}"/raw_variants.vcf \
    --select-type INDEL \
    -O "${results}"/raw_indels.vcf

# Filter INDELs
gatk VariantFiltration \
	-R "${resources}"/Homo_sapiens_assembly38.fasta \
	-V "${results}"/raw_indels.vcf \
	-O "${results}"/filtered_indels.vcf \
	-filter-name "QD_filter" -filter "QD < 2.0" \
	-filter-name "FS_filter" -filter "FS > 200.0" \
	-filter-name "SOR_filter" -filter "SOR > 10.0" \
	-genotype-filter-expression "DP < 10" \
	-genotype-filter-name "DP_filter" \
	-genotype-filter-expression "GQ < 10" \
	-genotype-filter-name "GQ_filter"

# Selecting
# gatk SelectVariants \
# 	--exclude-filtered \
# 	-V "${results}"/filtered_indels.vcf \
# 	-O "${results}"/analysis_ready_indels.vcf

echo "--------------"
echo "Merge VCFs ..."
echo "--------------"

# File for Exomiser! 
gatk MergeVcfs \
    -I "${results}"/filtered_snps.vcf \
    -I "${results}"/filtered_indels.vcf \
    -O "${results}"/merge_to_exomiser_FINAL.vcf