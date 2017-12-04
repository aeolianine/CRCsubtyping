#' iNMF (Schlicker) subtypes signatures.
#'
#' The genes defining each of the 5 subtypes in the Schlicker (iNMF) signature.
#'
#' @format List of 7 character vectors. Each vector contains entrez gene identifiers
#' \describe{
#'   \item{1}{entrez identifiers of signature genes for subtype 1}
#'   \item{2}{entrez identifiers of signature genes for subtype 2}
#'   \item{1.1}{entrez identifiers of signature genes for subtype 1.1}
#'   \item{1.2}{entrez identifiers of signature genes for subtype 1.2}
#'   \item{1.3}{entrez identifiers of signature genes for subtype 1.3}
#'   \item{2.1}{entrez identifiers of signature genes for subtype 2.1}
#'   \item{2.2}{entrez identifiers of signature genes for subtype 2.2}
#'   ...
#' }
#' @source \url{https://bmcmedgenomics.biomedcentral.com/articles/10.1186/1755-8794-5-66}
"inmf_signatures_entrez"

#' iNMF (Schlicker) subtypes signatures.
#'
#' The genes defining each of the 5 subtypes in the Schlicker (iNMF) signature.
#'
#' @format List of 7 character vectors. Each vector contains gene symbol identifiers
#' \describe{
#'   \item{1}{symbol identifiers of signature genes for subtype 1}
#'   \item{2}{symbol identifiers of signature genes for subtype 2}
#'   \item{1.1}{symbol identifiers of signature genes for subtype 1.1}
#'   \item{1.2}{symbol identifiers of signature genes for subtype 1.2}
#'   \item{1.3}{symbol identifiers of signature genes for subtype 1.3}
#'   \item{2.1}{symbol identifiers of signature genes for subtype 2.1}
#'   \item{2.2}{symbol identifiers of signature genes for subtype 2.2}
#'   ...
#' }
#' @source \url{internal; to be edited}
"inmf_signatures_sym"

#' The 786 genes identified in the PAM analysis by Sadandandam et al, including the 5 subtype centroids.
#'
#' Source: Sadanandam et al (2013). A colorectal cancer classification system that associates cellular
#'          phenotype and responses to therapy. Nature Medicine, 19(5), 619–25. doi:10.1038/nm.3175.
#'          Suppl Table 1, "PAM" sheet
#'
#' @format Data frame of 786 genes by 5 subtypes. Each entry represents the relative expression of that
#'         gene in the centroid of the corresponding subtype;
#' @source \url{https://images.nature.com/original/nature-assets/nm/journal/v19/n5/extref/nm.3175-S2.xls}
"sadanandam_786_genes"
