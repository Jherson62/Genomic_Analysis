# BiocManager::install("QDNAseq")
# BiocManager::install("Biobase")
# BiocManager::install("QDNAseq.hg19")

library(QDNAseq)
library(Biobase)
library(QDNAseq.hg19)
library(tidyverse)

# Directories
setwd("~/Documents/Bioinformatica/PGTA_RUNS/")
bam_folder <- "RUN42ab/bamFiles/"
output_dir <- "RUN42ab/karyos/"
bin_size <- 1000 # 1 Mb

# Crear directorio
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
cat(sample_names, sep = "\n")

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
set.seed(124)
read_counts <- applyFilters(
  read_counts
)

# estima la correccion
set.seed(124)
read_counts <- estimateCorrection(
  read_counts
)

# Aplica a todos los cromosomas (incluido sexual)
set.seed(124)
read_counts <- applyFilters(
  read_counts,
  chromosomes = NA
)

# Aplicar correccion GC/mapeabilidad
set.seed(124)
read_counts <- correctBins(
  read_counts
)

# Outliers
set.seed(124)
read_counts <- normalizeBins(read_counts)
copy_numbers_smooth <- smoothOutlierBins(read_counts)

# Segmentacion ------------------------------------------------------------
set.seed(124)
copy_numbers_segmented <- segmentBins(
  copy_numbers_smooth,
  transformFun = "sqrt"
)

set.seed(124)
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


# Tabla resultados --------------------------------------------------------

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

# Mosaicismo ------------------------------------------------------

# 2n = 1
seg_data <- assayDataElement(copy_numbers_called, "segmented")
bin_info <- fData(copy_numbers_called)

mosaico_data <- cbind(
  data.frame("chromosome" = bin_info$chromosome),
  seg_data
)

data <- mosaico_data |> 
  pivot_longer(cols = 2:dim(mosaico_data)[2],
               names_to = "Barcode",
               values_to = "Bins") |> 
  dplyr::mutate(
    twobins = 2 * Bins
  ) |> 
  dplyr::group_by(Barcode, chromosome) |> 
  dplyr::summarise(
    mean_bins = mean(twobins, na.rm = TRUE),
    .groups = "drop"
  )

sex_info <- data |> 
  dplyr::filter(chromosome == "Y") |> 
  dplyr::mutate(sex = ifelse(mean_bins > 1, "XY", "XX")) |> 
  select(Barcode, sex)

data_classified <- data |> 
  left_join(sex_info, by = "Barcode") |> 
  dplyr::mutate(
    expected_copies = case_when(
      chromosome %in% as.character(1:22) ~ 2,
      sex == "XX" & chromosome == "X" ~ 2,
      sex == "XX" & chromosome == "Y" ~ 0,
      sex == "XY" & chromosome %in% c("X", "Y") ~ 1
    ),
    porcentaje = (mean_bins - expected_copies) * 100,
    clasificacion = case_when(
      sex == "XX" & chromosome == "Y" ~ "Ruido",
      
      porcentaje >= 80 ~ "Ganancia completa",
      porcentaje <= -80 ~ "Perdida completa",
      
      porcentaje >= 50 & porcentaje < 80 ~ "Mosaico alto (Ganancia)",
      porcentaje <= -50 & porcentaje > -80 ~ "Mosaico alto (Perdida)",
      
      porcentaje >= 20 & porcentaje < 50 ~ "Mosaico bajo (Ganancia)",
      porcentaje <= -20 & porcentaje > -50 ~ "Mosaico bajo (Perdida)",
      
      TRUE ~ "Euploide"
    )
  )

# Solo valores de media
data_wide <- data |> 
  pivot_wider(
    names_from = Barcode,
    values_from = mean_bins
  )

data_wide <- data_wide |> 
  dplyr::mutate(
    chromosome = factor(chromosome,
                        levels = c(as.character(1:22), "X", "Y"))
  ) |> 
  arrange(chromosome)

# Solo porcentajes
data_wide <- column_to_rownames(data_wide, var = "chromosome")
porcentaje <- (data_wide - 2)*100
porcentaje <- rownames_to_column(porcentaje, var = "chromosome")

write.table(
  porcentaje,
  sep = "\t",
  file = paste0(output_dir, "QDNAseq_mosaic_percentaje.tsv"),
  row.names = FALSE
)

# Heat Map Categorico
library(ggplot2)

data_classified <- data_classified |>
  mutate(chromosome = factor(
    chromosome, levels = c(as.character(1:22), "X", "Y")
  ))

colores_pgta <- c(
  "Ganancia completa"       = "red", 
  "Mosaico alto (Ganancia)" = "#D6604D",
  "Mosaico bajo (Ganancia)" = "#F4A582",
  "Euploide"                = "#A6DBA0",
  "Mosaico bajo (Perdida)"  = "#92C5DE", 
  "Mosaico alto (Perdida)"  = "#4393C3",
  "Perdida completa"        = "blue",
  "Ruido"                   = "white" 
)

# 3. Generar el Heatmap
grafico_pgta <- ggplot(data_classified, 
                       aes(x = chromosome,
                           y = Barcode,
                           fill = clasificacion)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_manual(
    values = colores_pgta, 
    name = NULL,
    drop = FALSE
  ) +
  theme_minimal(base_size = 15) +
  labs(
    title = paste("PGT-A" , substr(bam_folder, 1, 5)),
    x = "Cromosoma",
    y = "Muestra"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_blank()
  ) +
  coord_fixed(ratio = 1) 

print(grafico_pgta)

ggsave(
  filename = paste0(output_dir, "Heatmap_PGTA.png"),
  plot = grafico_pgta,
  width = 16, 
  height = 8, 
  dpi = 300,
  bg = "white"
)

cat("\n---------------------------------\n")
cat("ANÁLISIS FINALIZADO\n")
cat("Resultados guardados en:\n")
cat(output_dir,"\n")
cat("---------------------------------\n")




