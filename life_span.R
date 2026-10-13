library(tidyverse)
library(BSgenome)
library(VariantAnnotation)
library(GenomicRanges)
library(vcfR)
library(Rsamtools)
library(stringr)

# Primer Intento
# hm: son los datos armonizados a CRCh38.
# son los que se deben utilizar para el análisis.
lf_tabla <- read.table(
  file = "../Documents/Bioinformatica/LifeSpan/PGS906/PGS000906_hmPOS_GRCh38.txt",
  skip = 19,
  sep = "\t",
  header = TRUE,
  quote = ""
)

lf_clean <- lf_tabla |> 
  dplyr::mutate(cromosoma = paste0("chr", hm_chr)) |>
  dplyr::mutate(posicion = paste0(cromosoma, "_", hm_pos)) |> 
  dplyr::select(
    cromosoma,
    hm_pos,
    posicion,
    hm_rsID,
    effect_allele,
    other_allele,
    effect_weight
  )

# Anotacion: snp-gene ------------------------------------------------------
# Solo para reporte de genes utilizados, no involucrado en el analisis
# library(biomaRt)
# 
# # Base de datos de snp
# snp_db <- useEnsembl(
#   biomart = "snps",
#   dataset = "hsapiens_snp"
# )
# 
# # Atrbutos para anotar
# # hg38
# attrs <- listAttributes(snp_db)
# attrs[grep("gene", attrs$name, ignore.case = TRUE), ]
# 
# lf_snps <- lf_tabla$hm_rsID
# 
# anotacion_snp <- getBM(
#   attributes = c("refsnp_id", "associated_gene"),
#   filters = "snp_filter",
#   values = lf_snps,
#   mart = snp_db
# )
# 
# # Genes incluidos
# anotacion_snp[anotacion_snp == ""] <- NA
# anotacion_snp[anotacion_snp == " "] <- NA
# 
# genes <- anotacion_snp[!is.na(anotacion_snp$associated_gene),
#                        "associated_gene"]
# 
# genes <- unique(trimws(unlist(strsplit(genes, "[,/]"))))
# 
# genes <- genes[!tolower(genes) %in% c("intergenic", "nrxnr", "9p21")]
# genes <- sort(genes)
# 
# cat(genes, sep = ", ")
# cat(lf_snps, sep = ", ")

# Prueba con genoma --------------------------------------------------------
# bgzip (para comprimir)
# tabix -p vcf <file> (para indexar)

vcf_file <- vcfR::read.vcfR(
  file = "../Documents/Exomas/GEN2600/results/snps_annotated.vcf.gz"
  )

vcf_fix <- as.data.frame(vcf_file@fix, stringAsFactor = FALSE)
vcf_gt <- as.data.frame(
  vcfR::extract.gt(vcf_file, element = "GT")
)

vcf_df <- cbind(
  vcf_fix, vcf_gt
)

vcf_df <- vcf_df |> 
  dplyr::mutate(posicion = paste0(CHROM, "_", POS)) |> 
  dplyr::select(-INFO)

# Merge
data_longevity <- lf_clean |> 
  left_join(vcf_df, by = "posicion")

head(data_longevity, 6)

# Completamos vacíos
filas_na <- is.na(data_longevity$CHROM)

ranges_faltas <- GRanges(
  seqnames = data_longevity$cromosoma[filas_na],
  ranges = IRanges(start = data_longevity$hm_pos[filas_na],
                   end = data_longevity$hm_pos[filas_na])
)

fasta <- Rsamtools::FaFile(
  "../Documents/Exomas/recursos/Homo_sapiens_assembly38.fasta"
)

nucleotides <- getSeq(
  fasta,
  ranges_faltas
)
nucleotides_str <- as.character(nucleotides)

# Rellenamos
data_longevity$CHROM[filas_na] <- names(nucleotides_str)
data_longevity$REF[filas_na] <- nucleotides_str
data_longevity$GEN2600[filas_na] <- "0/0"

# LONGEVIDAD SCORE --------------------------------------------------------
data_score <- data_longevity |> 
  dplyr::select(
    effect_allele,
    effect_weight,
    REF,
    ALT,
    Genotype = GEN2600
  )

data_score <- data_score |> 
  dplyr::mutate(
    alelos = case_when(
      Genotype == "0/0" ~ paste0(REF, REF),
      Genotype == "0/1" ~ paste0(REF, ALT),
      Genotype == "1/0" ~ paste0(ALT, REF),
      Genotype == "1/1" ~ paste0(ALT, ALT)
    )
  ) |> 
  dplyr::mutate(
    numero = stringr::str_count(alelos, pattern = effect_allele)
  ) |> 
  dplyr::mutate(
    score = effect_weight * numero
  )

# Puntaje de mi paciente
sum(data_score$score)

# Puntaje ideal
data_max <- data_score |> 
  dplyr::mutate(dosis_max = case_when(
    effect_weight >= 0 ~ effect_weight*2,
    effect_weight < 0 ~ effect_weight*0
  ),
  dosis_min = case_when(
    effect_weight <= 0 ~ effect_weight*2,
    effect_weight > 0 ~ effect_weight*0
  )) |> 
  dplyr::select(effect_weight, dosis_max, dosis_min)

sum(data_max$dosis_max)
sum(data_max$dosis_min)

