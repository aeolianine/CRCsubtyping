#' Loading the Schlicker (iNMF) subtype signatures from the paper supplement.
#'
#' Reference: Schlicker et al (2012). Subtypes of primary colorectal tumors correlate with response to targeted
#'            treatment in colorectal cell lines. BMC Medical Genomics, 5(1), 66. doi:10.1186/1755-8794-5-66.
#'
#' @source \url{goo.gl/fghKBY},\url{goo.gl/vCzAwf}, \url{goo.gl/C5FpoE}, \url{goo.gl/uwABi1},
#'                             \url{goo.gl/pvzt33}, \url{goo.gl/b1eGqo}, \url{goo.gl/S6UhM3}
#'
#' @param annotation (character), can be gene 'symbol' or 'entrez' identifier; default is 'symbol'
#' @seealso \code{\link{loadSadanandamSignature}} and \code{\link{loadCMSgenes}}
#' @return a list of signatures
#' @export
#' @examples
#' sig = loadSchlickerSignature()
#' print(sig$'1.1')

loadSchlickerSignature = function(annotation = 'symbol'){

  stopifnot(annotation %in% c('symbol', 'entrez'))

  if (grepl('entrez',annotation)){
    data('inmf_signatures_entrez')
    return(inmf_signatures_entrez)
  }
  data("inmf_signatures_sym")
  return(inmf_signatures_sym)

}


#' Loading Sadanandam's subtype signatures.
#'
#' Reference: Sadanandam et al (2013). A colorectal cancer classification system that associates cellular phenotype and responses to therapy. Nature Medicine, 19(5), 619–25. doi:10.1038/nm.3175.
#' Note: In the paper these genes are refered to as "CRCassigner-786"
#' @param NULL, no inputs
#' @seealso \code{\link{loadSchlickerSignature}} and \code{\link{loadCMSgenes}}
#' @return data frame of 786 genes by 5 subtypes; entries are relative expression values in the 5 subtype centroids
#' @export
#' @examples
#' sigs = loadSadanandamSignature()
#' print(sigs$Inflammatory)

loadSadanandamSignature = function(){
  data('sadanandam_786_genes')
  return(sadanandam_786_genes)
}


#' Convert Entrez ID to gene symbols.
#' Note: uses the "org.Hs.eg.db" package.
#'
#' @param character vector of entrez ids
#' @return character vector of gene symbols, with names as entrez ids
#' @export
#' @examples
#' print(entrez2symbol('1281') == 'COL3A1')

entrez2symbol = function(geneids){

  if (!require(org.Hs.eg.db)){
    source('https://bioconductor.org/biocLite.R')
    biocLite('org.Hs.eg.db')
  }
  if (!require(annotate)){
    source('https://bioconductor.org/biocLite.R')
    biocLite('annotate')
  }
  return( annotate::getSYMBOL(geneids, data='org.Hs.eg') )

}

#' Converts gene symbols to entrez ids using org.Hs.egSYMBOL2EG
#'          from the org.Hs.eg.db package.
#'
#' @param geneNames - character vector
#' @return character vector, the entrez identifiers corresponding to geneName, or the geneName itself if it is not mapped
#' @export
#' @examples
#' print(symbol2entrez('KRAS'))
#' print(symbol2entrez('shouldnotbemapped'))
#' print(symbol2entrez(c('KRAS', 'AAAAAAAA', 'TP53', 'KRAS')))

symbol2entrez = function(geneNames){

  if (!require(org.Hs.eg.db)){
    source('https://bioconductor.org/biocLite.R')
    biocLite('org.Hs.eg.db')
  }

  # the entrez to symbol map can also be accessed like this:
  # links(org.Hs.egSYMBOL2EG[mappedkeys(org.Hs.egSYMBOL2EG)])[1:10,c('gene_id', 'symbol')]

  L = as.list(org.Hs.eg.db::org.Hs.egSYMBOL2EG[AnnotationDbi::mappedkeys(org.Hs.eg.db::org.Hs.egSYMBOL2EG)])
  L = L[geneNames]
  L = L[!is.na(names(L))]

  unmappedGenes = setdiff(geneNames, names(L))
  for (gene in unmappedGenes){ L[[gene]] = gene}
  names = names(L)
  L = as.character(L)
  names(L) = names

  # clean up
  rm(geneNames, unmappedGenes, gene, names)
  return(L)
}


#' Load random-forest genes from CMS classifier,
#' stored in CMSclassifier::finalModel$importance, with their entrez identifiers
#'
#' @param geneSymbol, boolean: should entrez identifiers be converted to gene symbols; default is FALSE
#' @seealso \code{\link{loadSchlickerSignature}} and \code{\link{loadSadanandamSignature}}
#' @return a matrix with CMS genes as rownames, and random forest importance scores per subtype;
#'         this is basically "CMSclassifier::finalModel$importance", with row names converted to
#'         gene symbols if geneSymbol=TRUE
#' @export
#' @examples
#' mat = loadCMSgenes(geneSymbol=FALSE)
#' print(head(mat))

loadCMSgenes = function(geneSymbol = FALSE){

  if (!require(CMSclassifier)){
    library(devtools)
    install_github("Sage-Bionetworks/CMSclassifier")
  }

  mat = CMSclassifier::finalModel$importance

  if (geneSymbol){
    rownames(mat) = entrez2symbol(rownames(mat))
  }

  remove(geneSymbol)
  return( mat )
}

#' Intersect a signature with the data's available genes
#'
#' Produces a signature with the same subtypes, but smaller sets of genes, corresponding to what is available in the data.
#'
#' @param genes - the genes available in the data, usually rownames(data)
#' @param sigs - a list of the genes in every subtype (according to this signature);
#'               or a character vector of genes, if that signature is not split into subtypes
#' @param verbose - logical: print how many genes are in the intersection or not, default is TRUE
#' @return intsig - the "intersection" signature, comprising of only signature genes available in "genes"
#' @export
#' @examples
#' intersectSignature(c('ARID1B', 'PTEN'),c('KRAS', 'TGFB1', 'ARID1B'), verbose=TRUE)

intersectSignature = function(genes, sig, verbose = TRUE){

  if (grepl('list', class(sig))){
    intsig = lapply(sig, function(x){ intersect(x, genes) })

  } else if (grepl('character', class(sig))){
    intsig = intersect(genes, sig)
  }

  if (verbose & grepl('list', class(sig))){
    library(knitr,quietly = TRUE)
    print('intersectSignature(): Percentage of genes available for each subtype:')
    tab = rbind(NULL, 100*as.numeric(lapply(intsig, length))/as.numeric(lapply(sig, length)))
    colnames(tab) = names(sig)
    print(kable(tab))

    remove(tab)
  } else if (verbose & grepl('character', class(sig))){
    print('intersectSignature(): Percentage of genes available for this signature:')
    print(100*length(intsig)/length(sig))
  }

  remove(genes, sig, verbose)
  return(intsig)
}

