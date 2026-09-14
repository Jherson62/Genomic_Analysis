############################################################
# QDNAseq Multi-sample Pipeline
# Optimizado para WGA y detección de aneuploidías
# Genoma: hg19
############################################################

# BiocManager::install("QDNAseq") #nolint
# BiocManager::install("Biobase") #nolint
# BiocManager::install("QDNAseq.hg19") #nolint

rm(list = ls())

library(QDNAseq)
library(Biobase)
library(QDNAseq.hg19)

setwd("/home/jherson/Documentos/Fertility/")
bam_folder <- "RUN41/"
output_dir <- "RUN41/karios/"
bin_size <- 1000   # 1 Mb bins

bamfiles <- list.files(
  path = bam_folder,
  pattern = "\\.bam$",
  full.names = TRUE
)

sample_names <- gsub(".bam", "", basename(bamfiles))
cat(sample_names, sep = "\n")

bins <- getBinAnnotations(
  binSize = bin_size,
  genome = "hg19"
)

read_counts <- binReadCounts(
  bins,
  bamfiles,
  cache = FALSE # every new session, new bam files
)

sampleNames(read_counts) <- sample_names

# Raw copy number profile
#plot(
#  read_counts,
#  main = "1 kbp bins",
#  logTransform = FALSE,
#  ylim = c(-50, 200)
#)

#highlightFilters(
#  read_counts,
#  logTransform = FALSE,
#  residual = TRUE,
#  blacklist = TRUE
#)

############################################################
# 5️⃣ Filtrado y correcciones
############################################################

read_counts <- applyFilters(
  read_counts
)

# Correction for GC content and mappability
# Estimations
read_counts <- estimateCorrection(read_counts)

read_counts <- applyFilters(
  read_counts,
  chromosomes = NA
)

# Apply the GC/mappeability correction
read_counts <- correctBins(read_counts)

# Normalize and smooth outliers
read_counts <- normalizeBins(read_counts)
copy_numbers_smooth <- smoothOutlierBins(read_counts)

############################################################
# 6️⃣ Segmentación
############################################################

copy_numbers_segmented <- segmentBins(
  copy_numbers_smooth,
  transformFun = "sqrt"
)

copy_numbers_segmented <- normalizeSegmentedBins(
  copy_numbers_segmented
)

############################################################
# 7️⃣ Llamado de CNVs
############################################################

copy_numbers_called <- callBins(
  copy_numbers_segmented,
  organism = "human",
  method = "cutoff",
  cutoffs = log2(
    c(
      deletion = 0.5,
      loss = 1.5,
      gain = 2.5,
      amplification = 10
    ) / 2
  )
)

############################################################
# 8️⃣ Análisis por muestra (sexo, MAPD, aneuploidías)
############################################################

summary_table <- data.frame()

for (i in seq_along(sample_names)) {
  sample <- sample_names[i]

  cn <- assayDataElement(copy_numbers_called, "copynumber")[, i]
  chr <- fData(copy_numbers_called)$chromosome

  ########## DETECTAR SEXO ##########
  y_signal <- median(cn[chr == "Y"], na.rm = TRUE)

  if (!is.na(y_signal) && y_signal > 0.1) {
    sex  <- "XY"
  } else {
    sex <- "XX"
  }

  ########## DETECTAR ANEUPLOIDIAS ##########
  chr_median <- tapply(cn, chr, median, na.rm = TRUE)

  status <- rep("Normal", length(chr_median))
  status[chr_median < 1.5] <- "Monosomy"
  status[chr_median > 2.5] <- "Trisomy"

  aneuploid_chr <- names(chr_median)[status != "Normal"]

  if (length(aneuploid_chr) == 0) {
    aneuploidy <- "Euploid"
  } else {
    aneuploidy <- paste(
      aneuploid_chr,
      status[status != "Normal"],
      collapse = "; "
    )
  }

  ########## MAPD (Calidad) ##########
  log2ratios <- assayDataElement(copy_numbers_smooth, "copynumber")[, i]
  mapd <- median(abs(diff(log2ratios)), na.rm = TRUE)

  ########## RESULTADOS ##########
  summary_table <- rbind(
    summary_table,
    data.frame(
      Sample = sample,
      Sex = sex,
      MAPD = round(mapd, 3),
      Aneuploidy = aneuploidy
    )
  )
}

############################################################
# 9️⃣ Guardar tabla resumen
############################################################

write.table(
  summary_table, sep = "\t",
  file = paste0(output_dir, "PGTA_summary_results.tsv"),
  row.names = FALSE
)

############################################################
# 🔟 Exportar tabla completa de bins
############################################################

cn_matrix <- assayDataElement(copy_numbers_called, "copynumber")
bin_info <- fData(copy_numbers_called)

results_table <- cbind(
  bin_info[, c("chromosome", "start", "end", "gc", "mappability")],
  cn_matrix
)

write.table(
  results_table, sep = "\t",
  file = paste0(output_dir, "QDNAseq_bins_results.tsv"),
  row.names = FALSE
)

############################################################
# 1️⃣1️⃣ Reporte PDF (un gráfico por muestra)
############################################################

pdf(
  paste0(output_dir, "QDNAseq_PGTA_Report.pdf"),
  width = 14,
  height = 8
)

for (i in seq_along(sample_names)) {

  sample_obj <- copy_numbers_called[, i]   # extrae solo una muestra

  plot(
    sample_obj,
    main = paste(
      sample_names[i],
      "| Sexo:", summary_table$Sex[i],
      "| MAPD:", summary_table$MAPD[i]
    )
  )
}

dev.off()

############################################################
# 1️⃣2️⃣ Plot de sesgo GC
############################################################

pdf(paste0(output_dir, "GC_bias_plot.pdf"))

plot(read_counts, type = "gc")

dev.off()

############################################################

cat("\n---------------------------------\n")
cat("ANÁLISIS FINALIZADO\n")
cat("Resultados guardados en:\n")
cat(output_dir,"\n")
cat("---------------------------------\n")
