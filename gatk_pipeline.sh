#!/bin/bash
# exit when any command fails
set -eou pipefail

# Conda initialize
# shellcheck disable=SC1091
source "$HOME/miniconda3/etc/profile.d/conda.sh"

for i in $(seq -w 1 4); do

sample="EX260${i}"

# Directories
data="$HOME/Documentos/Fertility/Exomas/${sample}/data"
quality="$HOME/Documentos/Fertility/Exomas/${sample}/quality"
resources="$HOME/Documentos/Fertility/Exomas/resources"
aligned="$HOME/Documentos/Fertility/Exomas/${sample}/aligned"
results="$HOME/Documentos/Fertility/Exomas/${sample}/results"
stats_vcf="${results}/stats_vcf"
tmp_dir="$HOME/Documentos/Fertility/Exomas/${sample}/tmp_dir"

# create directory if necessary
mkdir -p "$quality" "$aligned" "$results" "$data" "$stats_vcf" "$tmp_dir"

# Files
ref="${resources}/Homo_sapiens_assembly38.fasta"
snpdb="${resources}/Homo_sapiens_assembly38.dbsnp138.vcf"

# Java options
JAVA_OPTS="-Xms8G -Xmx8G -XX:+UseG1GC -XX:+UseStringDeduplication -Djava.io.tmpdir=${tmp_dir}"

# Enviroment 1: BWA
conda activate NGStools

echo "----------------"
echo " Quality Control"
echo "----------------"

fastqc "${data}"/* -o "${quality}"
multiqc "${quality}" -o "${quality}" --force

echo "---------------------------------------"
echo " Map to reference using BWA-MEM"
echo "---------------------------------------"

bwa mem -t 8 -R "@RG\tID:${sample}\tPL:ILLUMINA\tSM:${sample}" \
    "${ref}" \
    "${data}/${sample}.cleaned_1.fastq.gz" \
    "${data}/${sample}.cleaned_2.fastq.gz" | \
samtools sort -@ 4 -o "${aligned}/${sample}_paired.bam"

#Enviroment 2: GATK
conda activate gatk_env

echo "------------------------------------"
echo "Mark Duplicates and add Tags ..."
echo "------------------------------------"

gatk --java-options "${JAVA_OPTS}" MarkDuplicatesSpark \
    -I "${aligned}/${sample}_paired.bam" \
    -O "${aligned}/${sample}_sort_dedup.bam" \
    --temp-dir "${tmp_dir}" \
    --spark-master "local[8]" #local threads 

gatk SetNmMdAndUqTags \
    -I "${aligned}/${sample}_sort_dedup.bam" \
    -O "${aligned}/${sample}_sort_dedup_tag.bam" \
    -R "${ref}"

echo "------------------------------------"
echo "Base Quality Score Recalibration ..."
echo "------------------------------------"

# make a math model
 gatk --java-options "${JAVA_OPTS}" BaseRecalibrator \
    -I "${aligned}/${sample}_sort_dedup_tag.bam" \
    -R "${ref}" \
    --known-sites "${snpdb}" \
    -O "${data}/recal_data.table"

# apply the math model
gatk --java-options "${JAVA_OPTS}" ApplyBQSR \
    -I "${aligned}/${sample}_sort_dedup_tag.bam" \
    -R "${ref}" \
    --bqsr-recal-file "${data}/recal_data.table" \
    -O "${aligned}/${sample}_sort_dedup_tag_bqsr.bam"

echo "-------------------------------------------"
echo "Collect Alignment & Insert Size Metrics ..."
echo "-------------------------------------------"


gatk --java-options "${JAVA_OPTS}" CollectAlignmentSummaryMetrics \
    -R "${ref}" \
    -I "${aligned}/${sample}_sort_dedup_tag_bqsr.bam" \
    -O "${aligned}/alignment_metrics.txt"


gatk --java-options "${JAVA_OPTS}" CollectInsertSizeMetrics \
    -I "${aligned}/${sample}_sort_dedup_tag_bqsr.bam" \
    -O "${aligned}/insert_size_metrics.txt" \
    -H "${aligned}/insert_size_histogram.pdf"


echo "-------------------"
echo "HaplotypeCaller ..."
echo "-------------------"

# Call Variants . . .
gatk --java-options "${JAVA_OPTS}" HaplotypeCaller \
    -R "${ref}" \
    -I "${aligned}/${sample}_sort_dedup_tag_bqsr.bam" \
    -O "${results}/raw_variants.vcf" \
    -L "${resources}/cromosomas_principales.list"


echo "-----------------"
echo "Separate SNPs ..."
echo "-----------------"

gatk SelectVariants \
    -R "${ref}" \
    -V "${results}/raw_variants.vcf" \
    --select-type SNP \
    -O "${results}/raw_snps.vcf"

# Filter SNPs
gatk VariantFiltration \
	-R "${ref}" \
	-V "${results}/raw_snps.vcf" \
	-O "${results}/filtered_snps.vcf" \
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

echo "-------------------"
echo "Separate INDELs ..."
echo "-------------------"

gatk SelectVariants \
    -R "${ref}" \
    -V "${results}/raw_variants.vcf" \
    --select-type INDEL \
    -O "${results}/raw_indels.vcf"

# Filter INDELs
gatk VariantFiltration \
	-R "${ref}" \
	-V "${results}/raw_indels.vcf" \
	-O "${results}/filtered_indels.vcf" \
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
    -I "${results}/filtered_snps.vcf" \
    -I "${results}/filtered_indels.vcf" \
    -O "${results}/merge_${sample}.vcf"


echo "------------------------"
echo "Annotation with GATK ..."
echo "------------------------"

# This part is just to get metrics of SNPS/INDELS

gatk --java-options "${JAVA_OPTS}" VariantAnnotator \
    -R "${ref}" \
    -V "${results}/merge_${sample}.vcf" \
    --dbsnp "${snpdb}" \
    -O "${results}/merge_annotated_${sample}.vcf"

bcftools view -H -f PASS -i 'ID!="."' "${results}/merge_annotated_${sample}.vcf" | \
    wc -l > "${stats_vcf}/reporte_SNPs_INDELs.txt"

bcftools stats -s - "${results}/merge_${sample}.vcf" > \
    "${stats_vcf}/all_stats${sample}.vchk"

plot-vcfstats -p "${stats_vcf}" "${stats_vcf}/all_stats${sample}.vchk"; done 
