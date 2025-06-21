# Gene datasets
# https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE264108&format=file&file=GSE264108%5Freadcounts%2Etxt%2Egz
id <- c("GSM8211616", "GSM8211617", "GSM8211618", "GSM8211619", "GSM8211620",
        "GSM8211621", "GSM8211622", "GSM8211623", "GSM8211624", "GSM8211625",
        "GSM8211626", "GSM8211627", "GSM8211628", "GSM8211629")
metadata <- data.frame(id, condition = rep(c("mTNBC", "HD"), each = 7)) %>%
    mutate(across(everything(), as.factor)) %>%
    column_to_rownames("id")
metadata$condition <- relevel(metadata_samples$condition, ref = "HD")
gene_counts <- select(gene_counts, seq(15))
colnames(gene_counts)[2:15] <- id

tnbc_neutrophils <- list(gene_counts = gene_counts, metadata = metadata)
use_data(tnbc_neutrophils, overwrite = TRUE)

---

# Biomart
genome_species <- c("hsapiens", "mmusculus")
genome_version <- 113

biomart_annotation <- map(
    genome_species,
    ~ {
    db_annot <- useEnsembl(
        biomart = "genes",
        dataset = paste0(., "_gene_ensembl"),
        version = genome_version
    )

    getBM(
        attributes = c("ensembl_gene_id", "external_gene_name", "description", "chromosome_name","start_position", "end_position", "gene_biotype", "entrezgene_id"),
        mart = db_annot
    ) %>%
        rename(gene_name = "external_gene_name") %>%
        as_tibble()
}) %>%
    set_names(paste(genome_species, genome_version, sep = "_"))

use_data(biomart_annotation, overwrite = TRUE)

---
  
surface_genes <- map(
    genome_species,
    ~ {
      db_annot <- useEnsembl(
        biomart = "genes",
        dataset = paste0(., "_gene_ensembl"),
        version = genome_version
      )
      
      getBM(
        attributes = c("hgnc_symbol", "ensembl_gene_id", "go_id", "name_1006", "definition_1006"),
        filters = "go",
        values = c("GO:0005887", "GO:0009986"),
        mart = ensembl
      ) %>%
        pull(ensembl_gene_id) %>%
        unique()
}) %>%
  set_names(paste(genome_species, genome_version, sep = "_"))

use_data(surface_genes, overwrite = TRUE)
