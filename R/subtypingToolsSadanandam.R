#' Subtype samples into Sadanandam's subtypes using Andreas's subtype() function.
#'
#' @param whichSamples - a data frame (gene by samples) to be subtyped
#' @param types - samples annotation (a character or numeric vector, size 1xncol(whichSamples))
#' @param plotHeatmaps - TRUE/FALSE, whether to create heatmap plots per subtype, default is FALSE
#' @param plotSilWidths - TRUE/FALSE, whether to create silhouette width plots per subtype, default is FALSE
#' @param directory - string, naming the directory in which to save the heatmap/silhouette width plots, default is '.'
#' @param filePrefix - string, indicating a prefix to be added to the heatmap/sil width plot names, default is ''
#' @param distance - what is the distance metric for the clustering
#' @return res - clustering and silhouette widths of the samples
#' @export
#' @examples
#' ...

subtypeSamplesSadanandam = function(whichSamples, types = NULL, plotHeatmaps = FALSE, plotSilWidths = FALSE, directory='.', filePrefix='', distance = 'spearman', silhouette=TRUE){
    # first intersect signatures
    intsig = intersectSignature(rownames(whichSamples), loadSadanandamSignature())
    res = subtype(whichSamples, intsig, distance = distance, silhouette=silhouette)
            
    if (plotSilWidths==TRUE){
        png(paste(directory, "/", filePrefix, "sw_sad.png", sep=""), width=4000, height=3000, res=300)
        plotSilhouetteWidths( res$silhouette )
        dev.off()
    }
    if (plotHeatmaps==TRUE){
        png(paste(directory, "/", filePrefix, ".png", sep=""), width=4000, height=3000, res=300)

        subtype_colors = c('red', 'magenta', 'blue', 'brown', 'green')
        names(subtype_colors) = c('Enterocyte', 'Goblet.like', 'Inflammatory', 'Stem.like', 'TA')
        dataType_colors = c('orange','skyblue')
        names(dataType_colors) = c('organoids','TCGA-RNAseq')

        anno_colors = list(subtypes = subtype_colors, geneSig = subtype_colors, datatypes = dataType_colors)
        
        createHeatmap(whichSamples, res$clustering, intsig, types = types, cellwidth=1, cellheight = 0.5, anno_colors = anno_colors)
        dev.off()
    }
    return(res)
}



#' Convert a simple classes character vector to a clustering object. 
#'
#' This produces an output similar to the iNMF() function output. This applies to the Sadanandam clusters only.
#'
#' @param whichSamples - data (columns/samples) to which the clustering applies
#' @param clust - character vector of classes, names(clust) are the samples
#' @param signature - default is loadAndreasSignature()
#' @param distance - what is the distance metric for the clustering
#' @return clustering object with silhouette widths
#' @export
#' @examples
#' ...

createClusteringDataStructureSadanandam = function(whichSamples, clust, distance = 'spearman'){
    
    # recreate the clustering object for all samples (same as the output of iNMF())
    res = list()
    res$clustering = list("Enterocyte"=c(), "Goblet.like"=c(), "Stem.like"=c(), "Inflammatory"=c(), "TA" = c())
    
    for (sample in names(clust)){
        
        if (clust[[sample]] == "Enterocyte"){ 
            res$clustering[["Enterocyte"]] = append(res$clustering[["Enterocyte"]], sample)
      
        } else if (clust[[sample]] == "Goblet.like"){ 
            res$clustering[["Goblet.like"]] = append(res$clustering[["Goblet.like"]], sample)
      
        } else if (clust[[sample]] == "Stem.like"){ 
            res$clustering[["Stem.like"]] = append(res$clustering[["Stem.like"]], sample)
      
        } else if (clust[[sample]] == "Inflammatory"){ 
            res$clustering[["Inflammatory"]] = append(res$clustering[["Inflammatory"]], sample)
      
        } else if (clust[[sample]] == "TA"){
            res$clustering[["TA"]] = append(res$clustering[["TA"]], sample)
        }
    }
  
    # compute silhouette widths .....................
    # first for step1
    samples = res$clustering[['Enterocyte']]
    samples = append( samples, res$clustering[['Goblet.like']])
    samples = append( samples, res$clustering[['Stem.like']])
    samples = append( samples, res$clustering[['Inflammatory']])
    samples = append( samples, res$clustering[['TA']] )

    clusters = c(rep(1, length(res$clustering[['Enterocyte']])), rep(2, length(res$clustering[['Goblet.like']])), rep(3, length(res$clustering[['Stem.like']])), rep(4, length(res$clustering[['Inflammatory']])), rep(5, length(res$clustering[['TA']]))  )
    names(clusters) = samples  # this is the clustering object

    clustersDict = list('1'= 'Enterocyte', '2'= 'Goblet.like', '3'= 'Stem.like', '4'= 'Inflammatory', '5'= 'TA')
    
    # intersect the signature ....
    signature = intersectSignature(rownames(whichSamples), loadSadanandamSignature() )
    
    dist.mat = 1 - cor(whichSamples[unlist(signature), samples], use="pairwise.complete.obs", method = distance)
    sil = silhouette(clusters, dist.mat)

    if ( sum(is.na(sil))==0 ){
        rownames(sil) = samples
        for (i in 1:dim(sil)[1]){
            sil[i,1] = clustersDict[[ as.character(sil[i,1]) ]]
            sil[i,2] = clustersDict[[ as.character(sil[i,2]) ]]
        }
    }
    res$silhouette = sil

    return(res)
}


#' Subtyping samples one at a time, using Sadanandam's subtypes.
#'
#' Right now the implementation uses the subtype() function to assign 1 sample to a subtype. The subtypes are Inflammatory, TA, Stem.like, Goblet.like and Enterocyte.
#' 
#' @param whichSamples - a data frame (gene by samples) to be subtyped
#' @param plotHeatmaps - TRUE/FALSE, whether to create heatmap plots per subtype, default is FALSE
#' @param plotSilWidths - TRUE/FALSE, whether to create silhouette width plots per subtype, default is FALSE
#' @param method - default 'maxmean': assigning cluster expression means to signature subtype expression means
#' @param directory - string, naming the directory in which to save the heatmap/silhouette width plots, default is '.'
#' @param filePrefix - string, indicating a prefix to be added to the heatmap/sil width plot names, default is ''
#' @param distance - what is the distance metric for the clustering
#' @return res - clustering and silhouette widths of the samples
#' @export
#' @examples
#' res = subtypeOneAtATimeSadanandam(organoidsT)

subtypeOneAtATimeSadanandam = function(whichSamples, types = NULL, plotHeatmaps = FALSE, plotSilWidths = FALSE, method = 'maxmean', directory='.', filePrefix = '', distance = 'spearman', silhouette=TRUE){
    
    intsig = intersectSignature(rownames(whichSamples), loadSadanandamSignature())
  
    # Perform max-mean subtyping
    clust = list()
  
    for (column in 1:ncol(whichSamples)){

        if (method=='maxmean'){
            # this is the only method implemented so far
            res = subtype(whichSamples[,column,drop=FALSE], signatures=intsig, silhouette = silhouette)
        }
        clust[[colnames(whichSamples[,column,drop=FALSE])]] = names( res$clustering )
    
    }
    remove(res)
    res = list(clustering=NULL, silhouette=NA)
    res$clustering = list("Stem.like"=c(), "Goblet.like"=c(), "Inflammatory"=c(), "TA"=c(), "Enterocyte"=c())
    
    for (sample in names(clust)){
        subT = clust[[sample]]
        res$clustering[[subT]] = append( res$clustering[[subT]], sample )
    }

    # compute silhouette widths
    samples = append(res$clustering$Inflammatory, res$clustering$Goblet.like)
    samples = append(samples, res$clustering$TA)
    samples = append(samples, res$clustering$Enterocyte)
    samples = append(samples, res$clustering$Stem.like)
    clusters = c(rep(1, length(res$clustering$Inflammatory)), rep(2, length(res$clustering$Goblet.like)), rep(3, length(res$clustering$TA)), rep(4, length(res$clustering$Enterocyte)), rep(5, length(res$clustering$Stem.like)) )
    names(clusters) = samples  # this is the clustering object
  
    dict = list('1'= 'Inflammatory', '2'= 'Goblet.like', '3'= 'TA', '4'= 'Enterocyte', '5'= 'Stem.like')
  
    dist.mat = 1 - cor(whichSamples[unlist(intsig), samples], use="pairwise.complete.obs", method = distance)
    sil = silhouette(clusters, dist.mat)

    if (sum(is.na(sil))==0){
        rownames(sil) = samples
        for (i in 1:dim(sil)[1]){
            sil[i,1] = dict[[ as.character(sil[i,1]) ]]
            sil[i,2] = dict[[ as.character(sil[i,2]) ]]
        }
    }
    res$silhouette = sil

    
    if (plotSilWidths==TRUE && sum(is.na(sil))==0){
        png(paste(directory, "/", filePrefix, "sw_sad.png", sep=""), width=4000, height=3000, res=300)
        plotSilhouetteWidths( res$silhouette )
        dev.off()
    }
    if (plotHeatmaps==TRUE){
        png(paste(directory, "/", filePrefix, ".png", sep=""), width=4000, height=3000, res=300)
        
        subtype_colors = c('red', 'magenta', 'blue', 'brown', 'green')
        names(subtype_colors) = c('Enterocyte', 'Goblet.like', 'Inflammatory', 'Stem.like', 'TA')
        dataType_colors = c('red','gray90')
        names(dataType_colors) = c('organoids','TCGA-RNAseq')

        if (is.null(types)){ dataType_colors = NULL}
        anno_colors = list(subtypes = subtype_colors, geneSig = subtype_colors, datatypes = dataType_colors)
        
        createHeatmap(whichSamples, res$clustering, intsig, types = types, cellwidth=1, cellheight = 0.6, anno_colors = anno_colors)

        dev.off()
    }

    
    return(res)
}


#' Apply the shrunken centroids method to assign a sample to a subtype. Use the Sadanandam signature here.
#'
#' Ref: Tibshirani, R. (2002). Diagnosis of multiple cancer types by shrunken centroids of gene expression. PNAS, 99(10), 6567–72.
#'
#' @param df - data frame/matrix of expression values (genes by samples)
#' @param types: samples annotation (optional, character vector)
#' @param distance - what is the distance metric for the clustering
#' @return final predicted classes for the samples in df
#' @export
#' @examples
#' ...

shrunkenCentroidsSadanandam = function(df, types = NULL, plotHeatmap = F, directory = '.', filePrefix = '', distance = 'spearman'){

    res = subtypeOneAtATimeSadanandam(df, distance = distance)

    # pick positive silhouette widths for training data
    classes = c()
    trainSamples = c()
    for (i in 1:nrow(res$silhouette)){

        sw = as.numeric( res$silhouette[i,"sil_width"] )
        if (sw < 0.1){ next }
        classes = append(classes, res$silhouette[i,"cluster"]  )
        trainSamples = append(trainSamples, rownames(res$silhouette)[i])
        
    }
    print(paste('shrunkenCentroidsSadanandam(): Training on',length(trainSamples),'samples out of',ncol(df),sep=' '))
    
    require(pamr)
    
    dataTrain = list(x = as.matrix(df[,trainSamples]), y=classes)
    trainFit <- pamr.train(dataTrain, n.threshold = 30)
    
    #cvfit <- pamr.cv(trainFit, dataTrain, nfold=10, n.threshold = 30)
    ind = max( which( trainFit$errors==min(trainFit$errors) ) )

    threshold = trainFit$threshold[ind]
    print(paste('shrunkenCentroidsSadanandam(): Min error threshold',threshold,sep=' '))

    p = pamr.predict(trainFit, df, threshold, type=c("class"))
    
    # export the resulting clustering
    clust = list()
    for (i in 1:ncol(df)){ clust[[ as.character(p[i]) ]] = append(clust[[ as.character(p[i]) ]], colnames(df)[i]) }

    if (plotHeatmap){
        png(paste(directory, "/", filePrefix, ".png", sep=""), width=4000, height=3000, res=300)

        subtype_colors = c('red', 'magenta', 'blue', 'brown', 'green')
        names(subtype_colors) = c('Enterocyte', 'Goblet.like', 'Inflammatory', 'Stem.like', 'TA')
        #dataType_colors = c('orange','lightgreen','skyblue')
        #names(dataType_colors) = c('organoids','TCGA-Agilent','TCGA-RNAseq')
        dataType_colors = c('orange','skyblue')
        names(dataType_colors) = c('organoids','TCGA-RNAseq')
        
        anno_colors = list(subtypes = subtype_colors, geneSig = subtype_colors, datatypes = dataType_colors)

        intsig = intersectSignature(rownames(df), loadSadanandamSignature())
        createHeatmap(df, clust, intsig, types = types, cellwidth=1, cellheight = 0.5, anno_colors = anno_colors)
        dev.off()
    }

    names(p) = colnames(df)
    res = createClusteringDataStructureSadanandam(df, p, distance = distance)
    return(res)
}


#' Use the centroid matrix from the Sadanandam paper to assign samples to the nearest centroid.
#'
#' Using Euclidean distance for sample to centroid distance.
#'
#' @param df - a data frame representing the data matrix
#' @param types: samples annotation (optional, character vector)
#' @param distance - what is the distance metric for the clustering
#' @return final predicted classes for the samples in df
#' @export
#' @examples
#' ...

nearestCentroidsSadanandam = function(df, types = NULL, plotHeatmaps = F, plotSilWidths=F, directory = '.', filePrefix = '', distance = 'spearman'){
    
    # read-in the centroid matrix
    t <- read.table('~/projects/CRCorganoids/subtypingSadanandam/signature.csv', skip = 6, header = T, sep=';', stringsAsFactors = FALSE)
    rownames(t) <- t[,1]  # make the gene names the rownames
    t = t[,-1]            # remove the gene names column
    t <- as.matrix(t)

    # intersect the genes in "t" with the genes in "df"
    commonGenes = intersect(rownames(df), rownames(t))
    t = t[commonGenes,]
    df = df[commonGenes,]

    #euc = function(x,y) { return(  sqrt( sum((x-y)^2) )  ) }
    dcor = function(x,y) { return( 1-cor(x,y,method=distance)  ) }

    res = list(clustering=NULL, silhouette=NA)
    res$clustering = list("Stem.like"=c(), "Goblet.like"=c(), "Inflammatory"=c(), "TA"=c(), "Enterocyte"=c())
    
    for (sample in 1:ncol(df)){
        #euc_dist = apply(t,2,euc,df[,sample])

        dcor_dist = apply(t,2,dcor,df[,sample,drop=FALSE])
        subtypeName = colnames(t)[ which.min(dcor_dist) ]
        res$clustering[[ subtypeName ]] = append(res$clustering[[ subtypeName ]], colnames(df)[sample])
    }


    # compute silhouette widths
    samples = append(res$clustering$Inflammatory, res$clustering$Goblet.like)
    samples = append(samples, res$clustering$TA)
    samples = append(samples, res$clustering$Enterocyte)
    samples = append(samples, res$clustering$Stem.like)
    clusters = c(rep(1, length(res$clustering$Inflammatory)), rep(2, length(res$clustering$Goblet.like)), rep(3, length(res$clustering$TA)), rep(4, length(res$clustering$Enterocyte)), rep(5, length(res$clustering$Stem.like)) )
    names(clusters) = samples  # this is the clustering object
  
    dict = list('1'= 'Inflammatory', '2'= 'Goblet.like', '3'= 'TA', '4'= 'Enterocyte', '5'= 'Stem.like')
  
    dist.mat = 1 - cor(df[rownames(t), samples], use="pairwise.complete.obs", method = distance)
    sil = silhouette(clusters, dist.mat)

    if ( sum(is.na(sil))==0 ){
        rownames(sil) = samples
        for (i in 1:dim(sil)[1]){
            sil[i,1] = dict[[ as.character(sil[i,1]) ]]
            sil[i,2] = dict[[ as.character(sil[i,2]) ]]
        }
    }
    res$silhouette = sil

    if (plotSilWidths==TRUE){
        png(paste(directory, "/", filePrefix, "sw_sad.png", sep=""), width=4000, height=3000, res=300)
        plotSilhouetteWidths( res$silhouette )
        dev.off()
    }
    if (plotHeatmaps==TRUE){
        png(paste(directory, "/", filePrefix, ".png", sep=""), width=4000, height=3000, res=300)
        
        subtype_colors = c('red', 'magenta', 'blue', 'brown', 'green')
        names(subtype_colors) = c('Enterocyte', 'Goblet.like', 'Inflammatory', 'Stem.like', 'TA')
        dataType_colors = c('orange','skyblue')
        names(dataType_colors) = c('organoids','TCGA-RNAseq')

        if (is.null(types)){ dataType_colors = NULL}
        anno_colors = list(subtypes = subtype_colors, geneSig = subtype_colors, datatypes = dataType_colors)

        intsig = intersectSignature(rownames(df), loadSadanandamSignature())
        createHeatmap(df, res$clustering, intsig, types = types, cellwidth=1.2, cellheight = 0.6, anno_colors = anno_colors)

        dev.off()
    }

    
    return(res)
}

#' Computes the number of mis-classified samples between two clusterings.
#'
#' This is implemented for results of the Sadanandam algorithm.
#'
#' @param resInmf1 - first clustering object to be compared
#' @param resInmf2 - second clustering object to be compared
#' @param samplesOfInterest - out of the samples in the data, the subset of samples for which the comparison is performed
#' @return list -
#' @export
#' @examples
#' ...

compareTwoClusteringsSadanandam = function(res1, res2, samplesOfInterest){

    c1 = res1$silhouette
  
    clusterings = c1[,1]
    silWidths = c1[,3]

    # make sure samples are ordered the same way
    stopifnot( sum( sort(names(clusterings)) != sort(names(silWidths)) )==0 )
    firstSet = cbind(clusterings[samplesOfInterest], silWidths[samplesOfInterest])
    colnames(firstSet) = c('cluster','sil_width')
  
    c1 = res2$silhouette
    clusterings = c1[,1]
    silWidths = c1[,3]

    stopifnot( sum( sort(names(clusterings)) != sort(names(silWidths)) )==0 )
    secondSet = cbind(clusterings[samplesOfInterest], silWidths[samplesOfInterest])
    colnames(secondSet) = c('cluster','sil_width')

    compare = cbind( firstSet[samplesOfInterest,], secondSet[samplesOfInterest,]  )
    isitmatching = compare[,1] == compare[,3]
    print(paste('compareTwoClusteringsSadanandam(): Mis-classified samples',sum(isitmatching==FALSE),'out of',length(isitmatching),sep=' '))

    compare = cbind( compare[samplesOfInterest,], isitmatching[samplesOfInterest]  )  # add clustering matches/mismatches
  
    return(list(compare, sum(isitmatching==FALSE)/length(isitmatching)))
}
