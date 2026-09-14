library(readxl)
library(tidyverse)
library(karyoploteR)

setwd("Documentos/Fertility/RUN41/")

data_bins_chr <- read.table(
  "karios/QDNAseq_bins_results.tsv",
  sep = "\t", header = TRUE
)

dim(data_bins_chr)

view(karyoPlot)

karyoPlot <- data_bins_chr %>%
  pivot_longer(
    cols = 6:24,
    names_to = "Barcode",
    values_to = "Bins"
  ) %>%
  mutate(
    bins = 2 * Bins
  ) %>%
  group_by(Barcode, chromosome) %>%
  summarise(
    mean_bins = mean(bins, na.rm = TRUE),
    .groups = "drop"
  )

karyo_wide <- karyoPlot %>%
  pivot_wider(
    names_from = Barcode,
    values_from = mean_bins
  )

karyo_wide <- karyo_wide %>%
  mutate(
    chromosome = factor(
      chromosome,
      levels = c(as.character(1:22), "X", "Y")
    )
  ) %>%
  arrange(chromosome)

openxlsx::write.xlsx(karyo_wide, "RUN41_mosaico.xlsx")

mydata <- toGRanges(data.frame(
  chr = paste0("chr", karyoPlot$chromosome),
  start = karyoPlot$start,
  end = karyoPlot$end,
  y = 2 * karyoPlot$Bins,
  mean_bins = karyoPlot$mean_bins
))

ymin <- floor(min(mydata$y, na.rm = T))

mydata

ymax <- ceiling(max(mydata$y, na.rm = T))

kp <- plotKaryotype(plot.type=4, ideogram.plotter = NULL,
                    labels.plotter = NULL)
kpAddCytobandsAsLine(kp)
kpAddChromosomeNames(kp, srt=45)


cols <- sapply(split(karyoPlot$color, paste0("chr",karyoPlot$chromosome)), unique)

kpPoints(kp,
        data = mydata,
         ymin = ymin,
         ymax = ymax,
        col=colByChr(mydata, colors = cols))

kpAxis(kp, ymin = ymin, ymax = ymax)
kpAbline(kp, h=2, ymin=ymin, ymax=ymax, lty=2, col="#666666")



library(dplyr)
library(tidyr)
library(writexl)

sexo <- "M"

summary_mosaic_all <- data_bins_chr %>%
  pivot_longer(
    cols = 6:24,
    names_to = "Barcode",
    values_to = "Bins"
  ) %>%

  mutate(
    Barcode = gsub(
      ".*(barcode[0-9]+).*",
      "\\1",
      Barcode
    )
  ) %>%

  filter(Bins > 0.1) %>%

  mutate(copy_number = 2 * Bins) %>%

  group_by(Barcode, chromosome) %>%
  summarise(
    CN = median(copy_number),
    .groups = "drop"
  ) %>%

  mutate(

    esperado = case_when(
      chromosome == "X" & sexo == "M" ~ 1,
      chromosome == "Y" & sexo == "M" ~ 1,
      chromosome == "X" & sexo == "F" ~ 2,
      chromosome == "Y" & sexo == "F" ~ 0,
      TRUE ~ 2
    ),

    mosaicismo_pct = (CN - esperado) * 100,

    CN = round(CN, 3),
    mosaicismo_pct = round(mosaicismo_pct, 1),
  ) %>%
  arrange(as.numeric(chromosome))

# Separar por barcode
lista_excel <- split(summary_mosaic_all, summary_mosaic_all$Barcode)

# Exportar múltiples hojas
write_xlsx(
  lista_excel,
  path = "barcode_mosaicismo.xlsx"
)
