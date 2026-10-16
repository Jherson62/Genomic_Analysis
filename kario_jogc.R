# BiocManager::install("QDNAseq")
# BiocManager::install("Biobase")
# BiocManager::install("QDNAseq.hg19")

library(QDNAseq)
library(Biobase)
library(QDNAseq.hg19)

# Directories
setwd("../Documents/Bioinformatica/PGTA_RUNS/")
bam_folder <- "RUN38/bamFiles/"
output_dir <- "RUN38/karyos/"
bin_size <- 1000 # 1 Mb

# Crear si no existe
if(!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

bamfiles <- list.files(
  path = bam_folder,
  pattern = "\\.bam$",
  full.names = TRUE
)

sample_names <- gsub(
  ".bam", "",
  basename(bamfiles)
)

# Bins
bins <- getBinAnnotations(
  binSize = bin_size,
  genome = "hg19"
)

read_counts <- binReadCounts(
  bins = bins,
  bamfiles = bamfiles,
  bamnames = sample_names,
  cache = FALSE
)

# Filtrado y correcciones -------------------------------------------------
# Por default filtra sexuales y MT
read_counts <- applyFilters(
  read_counts
)

# estima la correccion
read_counts <- estimateCorrection(
  read_counts
)

# Aplica a todos los cromosomas (incluido sexual)
read_counts <- applyFilters(
  read_counts,
  chromosomes = NA
)

# Aplicar correccion GC/mapeabilidad
read_counts <- correctBins(
  read_counts
)

# Outliers
read_counts <- normalizeBins(read_counts)
copy_numbers_smooth <- smoothOutlierBins(read_counts)

# Segmentacion ------------------------------------------------------------
copy_numbers_segmented <- segmentBins(
  copy_numbers_smooth,
  transformFun = "sqrt"
)

copy_numbers_segmented <- normalizeSegmentedBins(
  copy_numbers_segmented
)

# CNVs --------------------------------------------------------------------
# Se asume 2 copias. 
# Escala log2. Log2(#copias/2)
# Log2(2/2) = 0 (normal)

# < 0.5 (log2(0.5/2) = -2) Delecion homocigota: 0 copias
# 0.5 a 1.5 (log2(1.5/2) = -0.415) Delecion heterocigotica: 1 copia
# 1.5 a 2.5 (log2 entre -0.415 y 0.322) Normal (ronda el 0)
# 2.5 a 10 (log2(10/2) = 2.32) Ganancia (Trisomia)
# Amplificacion mas de 10 (log2 = 2.32)

# Luego de callBins, se etiqueta:
# Delecion: -2
# Del. heterocigotica: -1
# Normal: 0
# Ganancia: 1
# Amplificacion:2 

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

# Analisis por muestra ----------------------------------------------------

summary_table <- data.frame()

for (i in seq_along(sample_names)) {
  sample <- sample_names[i]
  
  # promedio en cada segmento/normalizado
  # 2n = 1 (mediana)
  seg <- assayDataElement(
    copy_numbers_called,
    "segmented"
  )[, sample]
  
  calls <- assayDataElement(
    copy_numbers_called,
    "calls"
  )[, sample]
  
  chr <- fData(
    copy_numbers_called
  )$chromosome
  
  # Deteccion de SEXO / ANEUPLOIDIAS SEXUALES
  x_signal <- median(
    seg[chr == "X"],
    na.rm = TRUE
  )
  y_signal <- median(
    seg[chr == "Y"],
    na.rm = TRUE
  )
  
  # Masculino XY: X ~ 0.5 & Y > 0.2
  # Femenino XX: X ~ 1.0 & Y ~ 0.0
  get_copy <- function(signal) {
    if (is.na(signal) || signal < 0.2) return(0)
    if (signal < 0.75) return(1)
    if (signal <= 1.25) return(2)
    return(3)
  }
  
  x_copies <- get_copy(x_signal)
  y_copies <- get_copy(y_signal)
  
  sex_karyotype <- switch(
    paste0("X", x_copies, "Y", y_copies),
    "X2Y0" = list(sex = "XX", status = "Female"),
    "X1Y1" = list(sex = "XY", status = "Male"),
    "X1Y0" = list(sex = "X0", status = "Turner syndrome"),
    "X2Y1" = list(sex = "XXY", status = "Klinefelter syndrome"),
    "X1Y2" = list(sex = "XYY", status = "Jacobs syndrome"),
    "X3Y0" = list(sex = "XXX", status = "Triple X syndrome"),
    list(sex = paste0(
      "X=", round(x_signal, 2),
      ", Y=", round(y_signal, 2)
    ), status = "Indeterminado")
  )
  
  # Deteccion de ANEUPLOIDIAS
  autosomes <- as.character(1:22)
  chr_autosomes <- chr[chr %in% autosomes]
  seg_autosomes <- seg[chr %in% autosomes]
  
  # Mediana normal = 1.0 (2n)
  chr_median <- tapply(seg_autosomes, chr_autosomes, median, na.rm = TRUE)
  
  status <- rep("Normal", length(chr_median))
  names(status) <- names(chr_median)
  status[chr_median < 0.75] <- "Monosomy"
  status[chr_median > 1.25] <- "Trisomy"
  
  # Chr con aneuploidias
  trisomy_chr <- status[status == "Trisomy"]
  monosomy_chr <- status[status == "Monosomy"]
  
  # Numero de cromosomas
  n_chr <- 46 + (length(trisomy_chr) - length(monosomy_chr))
  
  if (nchar(sex_karyotype$sex) == 3) {
    n_chr <- n_chr + 1
  } else if (sex_karyotype$sex == "X0") {
    n_chr <- n_chr -1
  } else {
    n_chr <- n_chr
  }
  
  # Reporte total
  aneuploidy <- c()
  if (length(trisomy_chr) == 0 & (sex_karyotype$sex %in% c("XX", "XY") & length(monosomy_chr) == 0 )
  ) {
    aneuploidy <- c(aneuploidy, "Euploid")
  } else if (!(sex_karyotype$sex %in% c("XX", "XY"))) {
    aneuploidy <- c(aneuploidy, paste0("TDS:", sex_karyotype$status))
  }
  
  if (length(trisomy_chr) != 0) {
    aneuploidy <- c(aneuploidy, paste("Trisomy:", paste(names(trisomy_chr), collapse = ", ")))
  }
  
  if (length(monosomy_chr) != 0) {
    aneuploidy <- c(aneuploidy, paste("Monosomy:", paste(names(monosomy_chr), collapse = ", ")))
  }
  
  aneuploidy <- paste(aneuploidy, collapse = " | ")
  
  # MAPD (Control de Calida)
  ratios <- assayDataElement(copy_numbers_smooth, "copynumber")[, sample]
  ratios[ratios <= 0] <- NA 
  log2_signal <- log2(ratios)
  
  mapd <- median(
    abs(diff(log2_signal)), na.rm = TRUE
  )
  
  # RESULTADOS
  summary_table <- rbind(
    summary_table,
    data.frame(
      Sample = sample,
      Chromosomas = n_chr,
      Sex = sex_karyotype$sex,
      MAPD = round(mapd, 2),
      Aneuploidies = aneuploidy
    )
  )
}

write.table(
  summary_table, sep = "\t",
  file = paste0(output_dir, "PGTA_summary_results.tsv"),
  row.names = FALSE
)

# PDF Reporte -------------------------------------------------------------

pdf(
  paste0(output_dir, "QDNAseq_PGTA_Report.pdf"),
  width = 14,
  height = 8
)

for (i in seq_along(sample_names)) {
  sample_obj <- copy_numbers_called[, i]
  
  plot(
    sample_obj,
    main = paste(
      sample_names[i],
      "| Chr:", summary_table$Chromosomas[i],
      "| Sexo:", summary_table$Sex[i],
      "| Dx:", summary_table$Aneuploidies[i],
      "| MAPD:", summary_table$MAPD[i]
    )
  )
}

dev.off()






