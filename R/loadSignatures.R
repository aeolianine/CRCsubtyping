#' Loading the Schlicker (iNMF) subtype signatures as a list.
#' 
#' Original data source: "/srv/nfs4/medoid-bulk/NKI/a.schlicker/AZTS/RDATA_PRE_RENAME/inmf_signatures_sym.rdata".
#' Reference: Schlicker et al (2012). Subtypes of primary colorectal tumors correlate with response to targeted treatment in colorectal cell lines. BMC Medical Genomics, 5(1), 66. doi:10.1186/1755-8794-5-66.
#'
#' @param no inputs needed
#' @seealso \code{\link{loadSadanandamSignature}} which is another function loading signatures
#' @return a list of signatures
#' @export
#' @examples
#' sigs = loadSchlickerSignature()
#' print(sigs$'1.1')

loadSchlickerSignature = function(){

    load("~/projects/CRCorganoids/subtypingSchlicker/inmf_signatures_sym.rdata")    
    sigs = list('1'=ts.az.clust1.sym,'2'=ts.az.clust2.sym,'1.1'=ts.az.clust1.1.sym, '1.2'=ts.az.clust1.2.sym, '1.3'=ts.az.clust1.3.sym, '2.1'=ts.az.clust2.1.sym, '2.2'=ts.az.clust2.2.sym)

    # some tests to make sure data loaded properly
    stopifnot(length(sigs)==7)
    stopifnot(names(sigs)==c("1","2","1.1","1.2","1.3","2.1","2.2"))
    stopifnot(length(sigs$'1.3')==216)
    
    return(sigs)
}



#' Loading Sadanandam's subtype signatures as a list.
#' 
#' Reference: Sadanandam et al (2013). A colorectal cancer classification system that associates cellular phenotype and responses to therapy. Nature Medicine, 19(5), 619–25. doi:10.1038/nm.3175.
#'
#' @param no inputs needed
#' @seealso \code{\link{loadAndreasSignature}} which is another function loading signatures
#' @return a list of signatures
#' @export
#' @examples
#' sigs = loadSadanandamSignature()
#' print(sigs$Inflammatory)

loadSadanandamSignature = function(){
    
    sigs = list('Enterocyte' = c(), 'TA' = c(), 'Stem.like' = c(), 'Inflammatory'=c(), 'Goblet.like' = c())
    t <- read.table('~/projects/CRCorganoids/subtypingSadanandam/signature.csv', skip = 6, header = T, sep=';', stringsAsFactors = F)

    for (i in 1:dim(t)[1]){                             # for all columns
        ind = which(t[i,]==max(as.numeric(t[i,2:6])) )    # find max score
        subT = colnames(t)[ind]                           # subtype corresponding to max score
        sigs[[subT]] = append( sigs[[subT]], t[i,1] )   # assign sample to max-score subtype above
    }

    # some tests to make sure data loaded properly
    stopifnot(length(sigs)==5)
    stopifnot(names(sigs)==c('Enterocyte', 'TA', 'Stem.like', 'Inflammatory', 'Goblet.like'))
    return(sigs)
}



#' Plotting the scores from the signature table (CRCassigner-786) as a heatmap
#' 
#' Reference: Sadanandam et al (2013). A colorectal cancer classification system that associates cellular phenotype and responses to therapy. Nature Medicine, 19(5), 619–25. doi:10.1038/nm.3175
#'
#' @param no inputs needed
#' @seealso \code{\link{loadSadanandamSignature}} - this is a function that actually loads the signature as a list
#' @export
#' @examples
#' plotSadanandamSignatureScores()

plotSadanandamSignatureScores = function(){
    
    t <- read.table('~/projects/CRCorganoids/subtypingSadanandam/signature.csv', skip = 6, header = T, sep=';', stringsAsFactors = F)
    rownames(t) <- t[,1]
    t <- t[,-1]
    require(gplots)
    require(RColorBrewer)
    palette = colorRampPalette(c('blue','white'))(n=1000)
    heatmap.2(as.matrix(t), dendrogram = c('none'), scale=c('row'), col = palette, cexRow = 0.1, las=1, cexCol = 1)

    # attempt to plot column labels on top, didn't work
    #axis(3, 1:ncol(t), labels = colnames(t), las = 2, tick = 0, cex.axis = 1) 
}



#' Comparing signatures: overlapping genes per subtype
#' 
#' This function computes the overlapping genes per signature subtype for two signatures.
#' This is an indirect comparison of signature similarity since genes expression scores can be correlated.
#'
#' @param sig1 the first signature
#' @param sig2 the second signature
#' @return matrix of number of common genes per subtype
#' @export
#' @examples
#' signaturesGeneOverlap(loadAndreasSignature(), loadSadanandamSignature())

signaturesGeneOverlap = function(sig1,sig2){
   
    m = matrix(NA, nrow=length(sig1), ncol=length(sig2))
    i=1
    j=1
    for (s1 in names(sig1)){
        for (s2 in names(sig2)){
            m[i,j] = length( intersect(sig1[[s1]], sig2[[s2]]) )
            j = j + 1
        }
        i = i + 1
        j = 1
    }
    rownames(m) = names(sig1)
    colnames(m) = names(sig2)

    require(knitr)
    kable(m)
    return(m)
}


#' Intersect a signature with the data's available genes
#' 
#' Produces a signature with a same subtypes, but smaller sets of genes, corresponding to what is available in the data.
#'
#' @param genes - the genes available in the data, usually rownames(data)
#' @param sigs - a list of the genes in every subtype (according to this signature)
#' @param print - logical: print how many genes are in the intersection or not, default is TRUE
#' @return intsig - the "intersection" signature, comprising of only genes available in "genes"
#' @export
#' @examples
#' intsig = intersectSignature(rownames(data),sig,print=TRUE)
#' stopifnot(names(intsig)==names(sig))
#' print(length(intsig))

intersectSignature = function(genes, sig, print = TRUE){

    intsig = sig                             # copy for the subtype names
    for (name in names(sig)){
        intsig[[name]] = intersect(genes, sig[[name]]) 
    }

    if (print){
        for (name in names(intsig)){
            print(paste(length(intsig[[name]]),'genes in subtype',name,'out of',length(genes),'samples genes and',length(sig[[name]]),'signature genes, (',as.integer(100*length(intsig[[name]])/length(sig[[name]])),'% of signature)',sep=' '))
        }
    }

    # tests
    stopifnot(length(intsig)==length(sig))
    stopifnot(names(intsig)==names(sig))
    
    return(intsig)
}

