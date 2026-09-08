# Packages ----------------------------------------------------------------
library(tidyverse)

setwd("/home/jherson/Documentos/Fertility/Exomas/exomiser_result")

###################
id <- "EX2603"
###################

genes <- read.table(
  file = sprintf("%s_out/%s_exomiser.genes.tsv", id, id),
  sep = "\t",
  comment.char = "",
  quote = "",
  header = TRUE,
  fill = TRUE,
  na.strings = "NA",
  check.names = FALSE,
  stringsAsFactors = FALSE
)
colnames(genes)[1] <- gsub("^#", "", colnames(genes)[1])

variants <- read.table(
  file = sprintf("%s_out/%s_exomiser.variants.tsv", id, id),
  sep = "\t",
  comment.char = "",
  quote = "",
  header = TRUE,
  fill = TRUE,
  na.strings = "NA",
  check.names = FALSE,
  stringsAsFactors = FALSE
)
colnames(variants)[1] <- gsub("^#", "", colnames(variants[1]))

#   RANK = Ranking de genes. El #1 es el más asociado al fenotipo (si le damos algun fenotipo en el .yml)
#   MOI = mode of inheritance
#   p-value = mide el azar de tomar un gen del genoma y tenga ese ranking fenotipico.

# Filtering most important gene variants
top_variant <- variants |> 
  dplyr::filter(
    EXOMISER_ACMG_CLASSIFICATION %in% c("PATHOGENIC", "LIKELY_PATHOGENIC") |
    CLINVAR_PRIMARY_INTERPRETATION %in% c("PATHOGENIC", "LIKELY_PATHOGENIC")
  ) |> 
  dplyr::arrange(desc(EXOMISER_GENE_VARIANT_SCORE))

# Readable table
top_variant_clean <- top_variant |> 
  select(
    Gen = GENE_SYMBOL,
    Info = ID,
    Herencia = MOI,
    Genotipo = GENOTYPE,
    Efecto = FUNCTIONAL_CLASS,
    Nomenclatura_HGVS = HGVS,
    Frecuencia_Max = MAX_FREQ,
    Lugar_Frecuencia = MAX_FREQ_SOURCE,
    Clasificacion_ACMG = EXOMISER_ACMG_CLASSIFICATION,
    Evidencia_ACMG = EXOMISER_ACMG_EVIDENCE,
    Enfermedad_ACMG =EXOMISER_ACMG_DISEASE_NAME,
    ClinVar_Clasificacion = CLINVAR_PRIMARY_INTERPRETATION,
    ClinVar_ID = CLINVAR_VARIATION_ID,
    ClinVar_Star = CLINVAR_STAR_RATING
  )

# FUNCTIONAL CLASS --------------------------------------------------------

mut_class <- variants |> 
  dplyr::group_by(FUNCTIONAL_CLASS) |> 
  summarise(Count = n()) |> 
  mutate(Percent = (Count/length(variants$FUNCTIONAL_CLASS))*100) |> 
  mutate(Text = case_when(
    Percent < 1                                        ~ "Otros",
    FUNCTIONAL_CLASS == "disruptive_inframe_deletion"  ~ "Deleción en marco disruptiva",
    FUNCTIONAL_CLASS == "disruptive_inframe_insertion" ~ "Inserción en marco disruptiva",
    FUNCTIONAL_CLASS == "frameshift_elongation"        ~ "Cambio de marco (Elongación)",
    FUNCTIONAL_CLASS == "frameshift_truncation"        ~ "Cambio de marco (Truncamiento)",
    FUNCTIONAL_CLASS == "frameshift_variant"           ~ "Cambio de marco (Variante)",
    FUNCTIONAL_CLASS == "inframe_deletion"             ~ "Deleción en marco",
    FUNCTIONAL_CLASS == "inframe_insertion"            ~ "Inserción en marco",
    FUNCTIONAL_CLASS == "missense_variant"             ~ "Variante Missense",
    FUNCTIONAL_CLASS == "splice_acceptor_variant"      ~ "Sitio aceptor de splicing",
    FUNCTIONAL_CLASS == "splice_donor_variant"         ~ "Sitio donador de splicing",
    FUNCTIONAL_CLASS == "splice_region_variant"        ~ "Región de splicing",
    FUNCTIONAL_CLASS == "stop_gained"                  ~ "Variante Nonsense (Stop gained)",
    FUNCTIONAL_CLASS == "stop_lost"                    ~ "Pérdida de codón de parada",
    FUNCTIONAL_CLASS == "synonymous_variant"           ~ "Variante Sinónima",
    TRUE ~ FUNCTIONAL_CLASS
  )) |> 
  dplyr::group_by(Text) |> 
  summarise(
    Count = sum(Count),
    Percent = sum(Percent)
  )
  
mut_class$Text <- relevel(as.factor(mut_class$Text),
                          ref = "Otros")

colors_class <- c(
  "Variante Sinónima"               = "#4A90E2",
  "Variante Missense"               = "#F39C12",
  
  "Variante Nonsense (Stop gained)" = "#D73027",
  "Cambio de marco (Truncamiento)"  = "#800026",
  "Cambio de marco (Variante)"      = "#C51B7D",
  "Cambio de marco (Elongación)"    = "#E7298A",
  "Pérdida de codón de parada"      = "#E6550D",
  
  "Sitio aceptor de splicing"       = "#1B9E77",
  "Sitio donador de splicing"       = "#31A354",
  "Región de splicing"              = "#A1D99B", #Menor impacto

  "Deleción en marco"               = "#756BB1",
  "Inserción en marco"              = "#BCBDDC",
  "Deleción en marco disruptiva"    = "#54278F",
  "Inserción en marco disruptiva"   = "#9E9AC8",
  
  "Otros"                           = "grey20"
)

mut_class_graph <- ggplot(
  data = mut_class,
  aes(x = 2, y = Percent, fill = Text)
) +
  xlim(c(0.5, 2.5)) +
  geom_col(color = "white") +
  coord_polar(theta = "y", start = 0) +
  theme_void(base_size = 17) +
  scale_fill_manual(values = colors_class) +
  labs(
    title = "EFECTOS",
    fill = NULL
  ) +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5),
    legend.position = "right"
  )

ggsave(
  file = sprintf("%s_out/Graphs_efectos-donna.png", id),
  dpi = 300,
  width = 6,
  height = 4,
  plot = mut_class_graph
)

# ACMG - CLASS ------------------------------------------------------------
acmg_class <- variants |> 
  dplyr::group_by(EXOMISER_ACMG_CLASSIFICATION) |> 
  summarise(Count = n()) |> 
  mutate(Percent = (Count/sum(Count))*100) |> 
  mutate(Text = case_when(
    EXOMISER_ACMG_CLASSIFICATION == "BENIGN" ~ "Benigno",
    EXOMISER_ACMG_CLASSIFICATION == "LIKELY_BENIGN" ~ "Probablemente Benigno",
    EXOMISER_ACMG_CLASSIFICATION == "LIKELY_PATHOGENIC" ~ "Probablemente Patogenico",
    EXOMISER_ACMG_CLASSIFICATION == "NOT_AVAILABLE" ~ "Sin evidencia",
    EXOMISER_ACMG_CLASSIFICATION == "PATHOGENIC" ~ "Patogenico",
    EXOMISER_ACMG_CLASSIFICATION == "UNCERTAIN_SIGNIFICANCE" ~ "Significado incierto",
    TRUE ~ EXOMISER_ACMG_CLASSIFICATION
  ))

acmg_class_graph <- ggplot(data = acmg_class, 
                           mapping = aes(x = 2, y = Percent, fill = Text)) +
  geom_col(color = "white") +
  coord_polar(theta = "y", start = 0) +
  xlim(0.5, 2.5) +
  theme_void(base_size = 17) +
  labs(
    title = "EVIDENCIA ACMG",
    fill = NULL
  ) +
  theme (
    plot.title = element_text(size = 20, hjust = 0.5)
  )
  
ggsave(
  file = sprintf("%s_out/Graphs_acmg-donna.png", id),
  dpi = 300,
  plot = acmg_class_graph
)


# PharmCAT Report ---------------------------------------------------------
library(tidyverse)
library(jsonlite)
library(purrr)

setwd("/home/jherson/Documentos/Fertility/Exomas/")
json_path <- "EX2604/results/pharmcat_result/reporte_final_pgx/ex2604_pgx_ready.vcf.preprocessed.report.json"

data <- fromJSON(
  txt = json_path,
  flatten = TRUE
)

str(data)





