#!/bin/bash
set -eou pipefail

# conda
source "$HOME/anaconda3/etc/profile.d/conda.sh"

gatk --java-options "-Xms22G -Xmx22G -XX:+UseG1GC -XX:+UseStringDeduplication" PrintReads \
    -I ../aligned/GEN2600_sort_dedup_tag.bam \
    -L chrM \
    --read-filter MateOnSameContigOrNoMappedMateReadFilter \
    --read-filter MateUnmappedAndUnmappedReadFilter \
    -O GEN2600_chrM.bam

 gatk --java-options "-Xms22G -Xmx22G -XX:+UseG1GC -XX:+UseStringDeduplication" RevertSam \
    -I GEN2600_chrM.bam \
    -O GEN2600_chrM_unmapped.bam \
    --OUTPUT_BY_READGROUP false \
    --VALIDATION_STRINGENCY LENIENT \
    --ATTRIBUTE_TO_CLEAR FT \
    --ATTRIBUTE_TO_CLEAR CO \
    --SORT_ORDER queryname \
    --RESTORE_ORIGINAL_QUALITIES false



 samtools fastq -T RG GEN2600_chrM_unmapped.bam | \
bwa mem -t 8 -K 100000000 -p -Y -R '@RG\tID:GEN2600\tSM:GEN2600' \
~/Documents/Exomas/recursos/chrM/Homo_sapiens_assembly38.chrM.fasta - | \
samtools view -b -@ 4 -o GEN2600_chrM_aligned.bam -

 samtools fastq -T RG GEN2600_chrM_unmapped.bam | \
bwa mem -t 8 -K 100000000 -p -Y -R '@RG\tID:GEN2600\tSM:GEN2600' \
~/Documents/Exomas/recursos/chrM/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta - | \
samtools view -b -@ 4 -o GEN2600_chrM_shifted_aligned.bam -

 gatk MergeBamAlignment \
    -UNMAPPED GEN2600_chrM_unmapped.bam \
    -ALIGNED GEN2600_chrM_aligned.bam \
    -O GEN2600_chrM_merged_aligned.bam \
    -R ~/Documents/Exomas/recursos/mito/Homo_sapiens_assembly38.chrM.fasta \
    --SORT_ORDER coordinate \
    --ADD_MATE_CIGAR true \
    --MAX_INSERTIONS_OR_DELETIONS -1 \
    --ALIGNER_PROPER_PAIR_FLAGS true

gatk MergeBamAlignment \
    -UNMAPPED GEN2600_chrM_unmapped.bam \
    -ALIGNED GEN2600_chrM_shifted_aligned.bam \
    -O GEN2600_chrM_shifted_merged_aligned.bam \
    -R ~/Documents/Exomas/recursos/mito/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta \
    --SORT_ORDER coordinate \
    --ADD_MATE_CIGAR true \
    --MAX_INSERTIONS_OR_DELETIONS -1 \
    --ALIGNER_PROPER_PAIR_FLAGS true

gatk MarkDuplicates \
    -I GEN2600_chrM_merged_aligned.bam \
    -O GEN2600_chrM_final.bam \
    -M GEN2600_chrM_dup_metrics.txt \
    --VALIDATION_STRINGENCY SILENT \
    --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500

gatk MarkDuplicates \
    -I GEN2600_chrM_shifted_merged_aligned.bam \
    -O GEN2600_chrM_shifted_final.bam \
    -M GEN2600_chrM_shifted_dup_metrics.txt \
    --VALIDATION_STRINGENCY SILENT \
    --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500

samtools index -@ 8 GEN2600_chrM_final.bam
samtools index -@ 8 GEN2600_chrM_shifted_final.bam 

# Llamado en la region no control
# referencia standard

gatk Mutect2 \
    -R /home/gseq/Documents/Exomas/recursos/mito/Homo_sapiens_assembly38.chrM.fasta \
    -I GEN2600_chrM_final.bam \
    -L chrM:576-16024 \
    -O GEN2600_non_control.vcf \
    --mitochondria-mode true \
    --annotation StrandBiasBySample \
    --max-reads-per-alignment-start 75 \
    --max-mnp-distance 0

gatk Mutect2 \
    -R /home/gseq/Documents/Exomas/recursos/mito/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta \
    -I GEN2600_chrM_shifted_final.bam \
    -L chrM:8025-9144 \
    -O GEN2600_control_shifted.vcf \
    --mitochondria-mode true \
    --annotation StrandBiasBySample \
    --max-reads-per-alignment-start 75 \
    --max-mnp-distance 0

gatk LiftoverVcf \
    -I GEN2600_control_shifted.vcf \
    -O GEN2600_control.vcf \
    -CHAIN ~/Documents/Exomas/recursos/mito/ShiftBack.chain \
    -REJECT GEN2600_control_rejected.vcf \
    -R ~/Documents/Exomas/recursos/mito/Homo_sapiens_assembly38.chrM.fasta

gatk MergeVcfs \
    -I GEN2600_non_control.vcf \
    -I GEN2600_control.vcf \
    -O GEN2600_chrM_merged.vcf

gatk MergeMutectStats \
    -stats GEN2600_non_control.vcf.stats \
    -stats GEN2600_control_shifted.vcf.stats \
    -O GEN2600_chrM_merged.vcf.stats

gatk FilterMutectCalls \
    -V GEN2600_chrM_merged.vcf \
    -R /home/gseq/Documents/Exomas/recursos/mito/Homo_sapiens_assembly38.chrM.fasta \
    -O GEN2600_chrM_filtered_base.vcf \
    -stats GEN2600_chrM_merged.vcf.stats \
    --mitochondria-mode true \
    --min-allele-fraction 0.01

gatk NuMTFilterTool \
    -V GEN2600_chrM_filtered_base.vcf \
    -R /home/gseq/Documents/Exomas/recursos/mito/Homo_sapiens_assembly38.chrM.fasta \
    -O GEN2600_chrM_final_filtered.vcf \
    --autosomal-coverage 30
