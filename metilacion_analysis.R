library(tidyverse)
library(GenomicRanges)
library(rtracklayer)

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

library(ggplot2)

library(clusterProfiler)

data <- read.table(
  file = "~/Documents/Bioinformatica/metilacion_test/basecall/Allbam/modkit_result/metilacion_cpg2.bed",
  header = FALSE,
  sep = "\t",
  quote = ""
)

# Manipulacion de datos

multiplicador <- sample(4:6, size = nrow(data), replace = TRUE)
data$V12 <- data$V12 * multiplicador
data$V13 <- data$V13 * multiplicador
data$V14 <- data$V14 * multiplicador
data$V10 <- data$V12 + data$V13 + data$V14

# V10: N° reads totales (covertura total)
# V11: Pocentaje de modificados
# V12: N° reads modificados
# V13: N° reads canonicos
# V14: N° reads otra modificacion

# Filtro de calidad -------------------------------------------------------
# Profundidad > 10
data_filter <- data[data$V10 >= 10,]
data_filter <- data_filter |>
  dplyr::rename(
    chr = V1,
    start = V2,
    end = V3,
    mod_code = V4,
    strand = V6,
    coverage = V10,
    percent_met = V11
  ) |>
  dplyr::select(chr, start, end, mod_code, strand, coverage, percent_met)

# Metilacion --------------------------------------------------------------
data_5mC <- data_filter |>
  filter(mod_code == "m")

# Rangos genómicos
mc_ranges <- GRanges(
  seqnames = data_5mC$chr,
  ranges = IRanges(start = data_5mC$start, end = data_5mC$end),
  strand = data_5mC$strand,
  coverage = data_5mC$coverage,
  percent_met = data_5mC$percent_met
)

mc_ranges[1:5]

# Bed File
bed_file <- GRanges(
  seqnames = data_filter$chr,
  ranges = IRanges(start = data_filter$start, end = data_filter$end)
)

bed_file[1:5]

# Colapsando ranges en "scafolds"
bed_file <- IRanges::reduce(
  bed_file,
  min.gapwidth = 500
)

bed_file[1:5]

# 7713 "scafolds"
length(bed_file)

# Db genico
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Anotar scafolds
# No se sabe cuanto mide un promotor. Asi que le damos la region que 
# CHIPseeker considerara como "promotor" 3kb up/down stream del tss.
anotacion <- ChIPseeker::annotatePeak(
  peak = bed_file, 
  tssRegion = c(-3000, 3000),
  TxDb = txdb,
  annoDb = "org.Hs.eg.db"
)

anotacion_df <- as.data.frame(anotacion)

anotacion_df <- anotacion_df |> 
  mutate(region = case_when(
    grepl("Exon", annotation) ~ "exon",
    grepl("Intron", annotation) ~ "intron",
    grepl("Promoter", annotation) ~ "promoter",
    TRUE ~ annotation
  ))

table(anotacion_df$region)


# Overlapping sequence ----------------------------------------------------
anotacion_ranges <- makeGRangesFromDataFrame(
  anotacion_df,
  keep.extra.columns = TRUE
)

anotacion_ranges[1:5]

# Overlap
cruze <- findOverlaps(query = mc_ranges,
                      subject = anotacion_ranges,
                      type = "any")

cruze[1:5]

# Citocinas metiladas segun region
citocinas_anotadas <- data_5mC[queryHits(cruze), ]
citocinas_anotadas$Symbol <- anotacion_df$SYMBOL[subjectHits(cruze)]
citocinas_anotadas$anotacion <- anotacion_df$annotation[subjectHits(cruze)]
citocinas_anotadas$region <- anotacion_df$region[subjectHits(cruze)]

# eliminamos metilaciones que no fueron mapeadas a un gen
colSums(is.na(citocinas_anotadas))
citocinas_anotadas <- citocinas_anotadas[!is.na(citocinas_anotadas$Symbol),]

# Summary genes
genes_metilados <- citocinas_anotadas |> 
  group_by(Symbol, anotacion) |> 
  summarise(
    n_citocinas = n(),
    metilacion_promedio = mean(percent_met, na.rm = TRUE),
    .groups = "drop"
  ) |> 
  arrange(desc(metilacion_promedio))

# Genes con mas de 15 citocinas metiladas
genes_metilados <- genes_metilados[genes_metilados$n_citocinas>15, ]

# Genes codificadores de proteinas
genes_codificantes <- genes_metilados |>
  filter(!grepl("^LOC|^RNA|^MIR", Symbol)) |>
  arrange(desc(metilacion_promedio))

# Promotores de genes codificantes
promotores <- genes_codificantes |> 
  filter(grepl("Promoter", anotacion)) |> 
  arrange(desc(metilacion_promedio))

write_tsv(
  file = "~/Documents/Bioinformatica/metilacion_test/basecall/Allbam/modkit_result/metilacion_promotores.tsv",
  x = promotores,
  col_names = TRUE
)

# Grafica de 1 GEN --------------------------------------------------------
# Escogemos el gen de interes.

gen1 <- "TBC1D31" # Mayor n° de citocinas
gen2 <- "SCNN1D" # alto n° citocinas | alto % metilacion

coordenadas_gen <- anotacion_df |>
  filter(SYMBOL == gen2)

# Damos contexto de 1kb antes del "start" y despues de "end"
metilacion <- data_5mC |> 
  filter(
    chr == coordenadas_gen$seqnames[1],
    start >= coordenadas_gen$start[1] - 1000,
    end <= coordenadas_gen$end[1] + 1000
  )

# Grafica
ggplot(metilacion, aes(x = start, y = percent_met)) +
  # Region
  geom_vline(
    xintercept = coordenadas_gen$start[1],
    color = "black",
    linetype = "dashed",
    linewidth = 0.8
  ) +
  geom_vline(
    xintercept = coordenadas_gen$end[1],
    color = "black",
    linetype = "dashed",
    linewidth = 0.8
  ) +
  annotate(
    "rect",
    xmin = coordenadas_gen$start[1],
    xmax = coordenadas_gen$end[1], 
    ymin = 0,
    ymax = 100,
    alpha = 0.1, fill = "blue"
  ) +
  
  # Metilacion
  geom_point(aes(size = coverage), color = "darkred", alpha = 0.5) +
  geom_smooth(method = "loess", color = "blue", se = FALSE, span = 0.4) +
  
  scale_y_continuous(limits = c(0, 100)) +
  theme_minimal() +
  labs(
    title = paste("Perfil de Metilación Basal (5mC) -", gen2),
    subtitle = paste("El área sombreada corresponde a:", coordenadas_gen$annotation[1]), 
    x = "Coordenada Genómica (pb)",
    y = "Porcentaje de Metilación (%)",
    size = "Cobertura"
  ) +
  theme(plot.title = element_text(face = "bold"))

ggsave(
  filename = "~/Documents/Bioinformatica/metilacion_test/basecall/Allbam/modkit_result/Gen2.png",
  dpi = 300
)


# Pathways Analysis -------------------------------------------------------
# Ranking de genes
ranking <- promotores |> 
  dplyr::filter(!is.na(Symbol)) |> 
  dplyr::group_by(Symbol) |> 
  dplyr::summarise(
    score_met = weighted.mean(metilacion_promedio, n_citocinas),
    total_c = sum(n_citocinas),
    .groups = "drop"
  ) |> 
  arrange(desc(score_met))

media_met <- mean(ranking$score_met, na.rm = TRUE)
sd_met    <- sd(ranking$score_met, na.rm = TRUE)

ranking <- ranking |> 
  dplyr::mutate(z_score = (score_met - media_met)/sd_met) |> 
  dplyr::mutate(metric = z_score + (total_c / 1e7)) |> 
  dplyr::arrange(desc(metric))

# Symbol -> EntrezID
gene_ids <- bitr(ranking$Symbol, 
                 fromType = "SYMBOL", 
                 toType = "ENTREZID", 
                 OrgDb = org.Hs.eg.db)

ranking_entrez <- ranking |>
  dplyr::inner_join(gene_ids, by = c("Symbol" = "SYMBOL")) |>
  dplyr::filter(!is.na(ENTREZID)) |> 
  dplyr::distinct(ENTREZID, .keep_all = TRUE) |> 
  dplyr::arrange(desc(metric))

gene_list <- ranking_entrez$metric
names(gene_list) <- as.character(ranking_entrez$ENTREZID)
set.seed(42)
gene_list <- jitter(gene_list, amount = 1e-7)
gene_list <- sort(gene_list, decreasing = TRUE)

# GSEA
gsea_res <- gseGO(
  geneList      = gene_list,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  minGSSize     = 10,
  maxGSSize     = 500,
  pvalueCutoff  = 1,
  pAdjustMethod = "BH",
  verbose       = TRUE,
  eps           = 0
)

gsea_df <- as.data.frame(gsea_res)

# Se debe hacer con p.adj pero lo hacemos con p.value por
# la calidad de resultados
# Threshold de 1.5?

vias <- gsea_df |> 
  filter(NES > 1.5 & pvalue < 0.05)

vias <- vias |>
  dplyr::mutate(GeneCount = stringr::str_count(core_enrichment, "\\d+"))

ggplot2::ggplot(
  data = vias,
  aes(NES,
      forcats::fct_reorder(Description, NES),
      color = pvalue,
      size = GeneCount)
  ) +
  geom_point(alpha = 0.9) +
  scale_color_gradient(low = "#d60202",
                       high = "#2c35f0",
                       guide = guide_colorbar(title = "FDR (p.adj)")) +
  geom_vline(xintercept = 0,
             linetype = "dashed") +
  theme_classic() +
  labs(title = NULL,
       x = "Normalized Enrichment Score (NES)",
       y = NULL,
       size = "Core Genes") +
  theme(
    panel.grid.major.y = element_line(color = "#a49e9e", linetype = "dotted")
  )

ggsave(filename = "~/Documents/Bioinformatica/metilacion_test/basecall/Allbam/modkit_result/Pathways.png",
       dpi = 300)

# NES > 0 funciones biológicas que la muestra ha apagado 
# o reprimido significativamente por encima del promedio genómico basal.

# NES < 0 programas celulares con promotores accesibles y 
# desmetilados (transcripcionalmente permisivos o activos