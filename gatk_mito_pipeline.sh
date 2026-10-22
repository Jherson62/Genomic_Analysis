#!/bin/bash
set -eou pipefail

# conda
source "$HOME/anaconda3/etc/profile.d/conda.sh"
conda activate gatk_env

# Directoios
tmp_dir="$HOME/Documents/Exomas/tmp_dir"
aligned="$HOME/Documents/Exomas/${sample}/aligned"
mitocondria="$HOME/Documents/Exomas/${sample}/mitocondria"
recursos="$HOME/Documents/Exomas/recursos"
results="${mitocondria}/results"

mkdir -p "${tmp_dir}" "${mitocondria}" "${results}"

# Files
ref_canonico="${recursos}/mito/Homo_sapiens_assembly38.chrM.fasta" 
ref_shifted="${recursos}/mito/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta" 
shift_chain="${recursos}/mito/ShiftBack.chain"

sample="GEN2600"

JAVA_OPT="-Xms22G -Xmx22G -XX:+UseG1GC -XX:+UseStringDeduplication -Djava.io.tmpdir=${tmp_dir}"

# --------------------------------
# 01. Extraer reads mitocondriales
# ---------------------------------
gatk --java-options "${JAVA_OPT}" PrintReads \
    -I "${aligned}/${sample}_sort_dedup_tag.bam" \
    -L chrM \
    --read-filter MateOnSameContigOrNoMappedMateReadFilter \
    --read-filter MateUnmappedAndUnmappedReadFilter \
    -O "${sample}_chrM.bam"

# --------------------------------
# 02. mapped BAM -> unmapped BAM
# ---------------------------------
# RevertSam borra atributos ligados a alineamiento.
# ATTRIBUTE_TO_CLEAR (atributos no directamente ligado con alineamiento para eliminar)

 gatk --java-options "-Xms22G -Xmx22G -XX:+UseG1GC -XX:+UseStringDeduplication" RevertSam \
    -I "${sample}_chrM.bam" \
    -O "${sample}_chrM_unmapped.bam" \
    --OUTPUT_BY_READGROUP false \
    --VALIDATION_STRINGENCY LENIENT \
    --ATTRIBUTE_TO_CLEAR FT \
    --ATTRIBUTE_TO_CLEAR CO \
    --ATTRIBUTE_TO_CLEAR XS \
    --ATTRIBUTE_TO_CLEAR MC \
    --SORT_ORDER queryname \
    --RESTORE_ORIGINAL_QUALITIES false

# --------------------------------
# 03. Realignmet
# --------------------------------

conda activate NGStools

# ALINEAMIENTO CON CANONICO
samtools fastq -T RG "${sample}_chrM_unmapped.bam" | \
	bwa mem -t 12 -K 100000000 -p -Y -R "@RG\tID:${sample}\tSM:${sample}" \
	"${ref_canonico}" - | \
	samtools view -b -@ 4 -o "${sample}_chrM_canonico_aligned.bam" -

# ALINEAMIENTO DESPLAZADO (SHIFTED)
samtools fastq -T RG "${sample}_chrM_unmapped.bam" | \
	bwa mem -t 12 -K 100000000 -p -Y -R "@RG\tID:${sample}\tSM:${sample}" \
	"${ref_shifted}" - | \
	samtools view -b -@ 4 -o "${sample}_chrM_shifted_aligne"

# --------------------------
# 04. Uniendo informacion 
# --------------------------
# info del uBAM mitocondrial -> alineamiento canonico
conda activate gatk_env

 gatk --java-options ${JAVA_OPT} MergeBamAlignment \
    -UNMAPPED "${sample}_chrM_unmapped.bam" \
    -ALIGNED "${sample}_chrM_canonico_aligned.bam" \
    -O "${sample}_chrM_merged_canonico_aln.bam" \
    -R "${ref_canonico}" \
    --SORT_ORDER coordinate \
    --ADD_MATE_CIGAR true \
    --MAX_INSERTIONS_OR_DELETIONS -1 \
    --ALIGNER_PROPER_PAIR_FLAGS true

 gatk --java-options ${JAVA_OPT} MergeBamAlignment \
    -UNMAPPED "${sample}_chrM_unmapped.bam" \
    -ALIGNED "${sample}_chrM_shifted_aligned.bam" \
    -O "${sample}_chrM_merged_shifted_aln.bam" \
    -R ~"${ref_shifted}" \
    --SORT_ORDER coordinate \
    --ADD_MATE_CIGAR true \
    --MAX_INSERTIONS_OR_DELETIONS -1 \
    --ALIGNER_PROPER_PAIR_FLAGS true

# --------------------------
# 05. MarkDuplicates 
# --------------------------

gatk --java-options ${JAVA_OPT} MarkDuplicates \
    -I "${sample}_chrM_merged_canonico_aln.bam" \
    -O "${sample}_chrM_canonico_markdup.bam" \
    -M "${sample}_chrM_canonico_dup_metrics.txt" \
    --VALIDATION_STRINGENCY SILENT \
    --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500

gatk --java-options ${JAVA_OPT} MarkDuplicates \
    -I "${sample}_chrM_merged_shifted_aln.bam" \
    -O "${sample}_chrM_shifted_markdup.bam" \
    -M "${sample}_chrM_shifted_dup_metrics.txt" \
    --VALIDATION_STRINGENCY SILENT \
    --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500

conda activate NGStools

samtools index -@ 8 "${sample}_chrM_canonico_markdup.bam"
samtools index -@ 8 "${sample}_chrM_shifted_markdup.bam" 

# --------------------------
# 06. Mutect2 
# --------------------------

conda activate gatk_env

# Canonico: Variant Calling excluyendo la region control
gatk --java-options ${JAVA_OPT} Mutect2 \
    -R "${ref_canonico}" \
    -I "${sample}_chrM_canonico_markdup.bam" \
    -L chrM:576-16024 \
    -O "${sample}_non_control.vcf" \
    --mitochondria-mode true \
    --annotation StrandBiasBySample \
    --max-reads-per-alignment-start 75 \
    --max-mnp-distance 0

# Shifted: Variant Calling en la region control desplazada
gatk --java-options ${JAVA_OPT} Mutect2 \
    -R "${ref_shifted}" \
    -I "${sample}_chrM_shifted_markdup.bam" \
    -L chrM:8025-9144 \
    -O "${sample}_control_shifted.vcf" \
    --mitochondria-mode true \
    --annotation StrandBiasBySample \
    --max-reads-per-alignment-start 75 \
    --max-mnp-distance 0

# Corregir cordenadas de desplazamiento del VCF en la region control
gatk LiftoverVcf \
    -I "${sample}_control_shifted.vcf" \
    -O "${sample}_control.vcf" \
    -CHAIN ~/Documents/Exomas/recursos/mito/ShiftBack.chain \
    -REJECT "${sample}_control_rejected.vcf" \
    -R "${ref_canonico}" 

# --------------------------
# 07. Unir VCF 
# --------------------------
gatk MergeVcfs \
    -I "${sample}_non_control.vcf" \
    -I "${sample}_control.vcf" \
    -O "${sample}_chrM_total.vcf"

gatk MergeMutectStats \
    -stats "${sample}_non_control.vcf.stats" \
    -stats "${sample}_control_shifted.vcf.stats" \
    -O "${sample}_chrM_total.vcf.stats"

# --------------------------
# 08. Filtros 
# --------------------------
# Filtro base.
gatk FilterMutectCalls \
    -V "${sample}_chrM_total.vcf" \
    -R "${ref_canonico}" \
    -O "${sample}_chrM_filtro_base.vcf" \
    -stats "${sample}_chrM_merged.vcf.stats" \
    --mitochondria-mode true \
    --min-allele-fraction 0.01

# Filtrar probables variantes nucleares
# profundidad genomica: 30X
gatk NuMTFilterTool \
    -V "${sample}_chrM_filtro_base.vcf" \
    -R "${ref_canonico}" \
    -O "${sample}_chrM_final_filtered.vcf" \
    --autosomal-coverage 30
