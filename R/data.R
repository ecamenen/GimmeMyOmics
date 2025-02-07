#' Bulk RNA-seq Data from Neutrophils in Triple-Negative Breast Cancer (TNBC) Patients and Healthy Donors
#'
#' A dataset containing bulk RNA sequencing (RNA-seq) data from neutrophils isolated from patients with metastatic triple-negative breast cancer (mTNBC) and healthy donors (HDs). The dataset highlights systemic immune alterations and changes in neutrophil functionality in TNBC.
#'
#' @format A `list` containing two elements:
#' \describe{
#'   \item{gene_counts}{A `tibble` with 23,567 rows (genes) and 15 columns:
#'     \describe{
#'       \item{ensembl_gene_id}{Character, the Ensembl gene identifier (e.g., ENSG00000000457, ENSG00000000460).}
#'       \item{GSM8211616, GSM8211617, ..., GSM8211629}{Numeric, the normalized gene expression counts for each sample.}
#'     }
#'   }
#'   \item{metadata}{A `tibble` with 14 rows (samples) and 2 columns:
#'     \describe{
#'       \item{id}{Character, the GEO sample identifier (e.g., GSM8211616, GSM8211617).}
#'       \item{condition}{Character, the group identifier (e.g., mTNBC, HD).}
#'     }
#'   }
#' }
#'
#' @source The dataset was created from RNA-seq data downloaded from the NCBI GEO database under the accession numbers:
#'   - GSM8211616, GSM8211617, GSM8211618, GSM8211619, GSM8211620, GSM8211621, GSM8211622 (mTNBC patients)
#'   - GSM8211623, GSM8211624, GSM8211625, GSM8211626, GSM8211627, GSM8211628, GSM8211629 (healthy donors, HDs).
#'   The data was generated as part of a study on the systemic immune landscape and neutrophil functionality in triple-negative breast cancer.
#'
#' @references
#' Bakker NAM, Garner H, van Dyk E, Champanhet E et al. "Triple-Negative Breast Cancer Shapes the Systemic Immune Landscape and Alters Neutrophil Functionality."
#'
#' @examples
#' # Load the dataset
#' data("tnbc_neutrophils")
#'
#' # Access gene counts
#' head(tnbc_neutrophils$gene_counts)
#'
#' # Access metadata
#' head(tnbc_neutrophils$metadata)
"tnbc_neutrophils"

#' Gene Annotation Data from Ensembl using biomaRt
#'
#' A dataset containing gene annotations for two species (*Homo sapiens* and *Mus musculus*) retrieved from Ensembl using the `biomaRt` package. The annotations include gene IDs, gene names, descriptions, chromosomal locations, and biotypes.
#'
#' @format A named `list` of length 2, where each element is a `tibble` containing gene annotations for a specific species and Ensembl version:
#' \describe{
#'   \item{hsapiens_113}{A `tibble` with 86,402 rows and 7 columns:
#'     \describe{
#'       \item{ensembl_gene_id}{Character, the Ensembl gene identifier (e.g., ENSG00000210049).}
#'       \item{gene_name}{Character, the official gene symbol (e.g., MT-TF).}
#'       \item{description}{Character, a brief description of the gene (e.g., "mitochondrially encoded tRNA-Phe").}
#'       \item{chromosome_name}{Character, the chromosome where the gene is located (e.g., MT, 1, 2).}
#'       \item{start_position}{Integer, the start position of the gene on the chromosome.}
#'       \item{end_position}{Integer, the end position of the gene on the chromosome.}
#'       \item{gene_biotype}{Character, the biotype of the gene (e.g., "protein_coding", "tRNA").}
#'     }
#'   }
#'   \item{mmusculus_113}{A `tibble` with 78,298 rows and 7 columns:
#'     \describe{
#'       \item{ensembl_gene_id}{Character, the Ensembl gene identifier (e.g., ENSMUSG00000064336).}
#'       \item{gene_name}{Character, the official gene symbol (e.g., mt-Tf).}
#'       \item{description}{Character, a brief description of the gene (e.g., "mitochondrially encoded tRNA-Phe").}
#'       \item{chromosome_name}{Character, the chromosome where the gene is located (e.g., MT, 1, 2).}
#'       \item{start_position}{Integer, the start position of the gene on the chromosome.}
#'       \item{end_position}{Integer, the end position of the gene on the chromosome.}
#'       \item{gene_biotype}{Character, the biotype of the gene (e.g., "protein_coding", "tRNA").}
#'     }
#'   }
#' }
#'
#' @source The data was retrieved from Ensembl (version 113).
#'
#' @examples
#' # Load the dataset
#' data("biomart_annotation")
#'
#' # Access human gene annotations
#' head(biomart_annotation$hsapiens_113)
#'
#' # Access mouse gene annotations
#' head(biomart_annotation$mmusculus_113)
"biomart_annotation"
