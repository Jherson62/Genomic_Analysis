#!/usr/bin/env  bash
# TEST: Genes Mola

set -euo pipefail

dorado basecaller ~/dorado_models/dna_r10.4.1_e8.2_400bps_sup@v5.2.0 pod5 \
	--modified-bases-models ~/dorado_models/dna_r10.4.1_e8.2_400bps_sup@v5.2.0_5mC_5hmC@v2 \
	--min-qscore 5 \
	--output-dir basecall \
	--reference ~/Documents/Bioinformatica/hg38/hg38.fa \
	--kit-name SQK-NBD114-24 \
	--barcode-both-ends \
	--trim all \
	--device cuda:all

# --barcode-both-end puede ser muy agresivo para adaptive

modkit pileup alineado_metilado.bam metilacion_cpg.bed \
	--ref ~/Documents/Bioinformatica/hg38/hg38.fa \
	--cpg \
	--modified-bases 5mC 5hmC \
	--threads 8

# Los de mayor a 15x 

awk -F'\t' '$10 >= 15' metilacion_cpg.bed > metilacion_targets.bed
