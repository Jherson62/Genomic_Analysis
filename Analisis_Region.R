library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)

# Barcode14 RUN42.
# Region del chrX
# PGTA: hg19

region <- GRanges(
  seqnames = "chrX",
  ranges = IRanges(start = 1, end = 6000000)
)

# Genes del hg19
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
genes_hg19 <- genes(txdb)

# Filtrar por region
genes_region <- subsetByOverlaps(genes_hg19,
                                 region)


columns(org.Hs.eg.db)
keytypes(org.Hs.eg.db)

# EntreID -> Symbol
genes_region$symbol <- mapIds(
  org.Hs.eg.db,
  keys = genes_region$gene_id,
  column = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first"
)

# Ver la lista de genes encontrados
df_genes <- as.data.frame(genes_region)

writexl::write_xlsx(
  df_genes, 
  path = "~/Downloads/genes_region.xlsx"
)

# Informacion OMIM
info_genes <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = genes_region$gene_id,
  column = c("SYMBOL", "GENENAME", "OMIM"),
  keytype = "ENTREZID",
  multiVals = "first"
)

writexl::write_xlsx(
  info_genes, 
  path = "~/Downloads/genes_informacion.xlsx"
)
