#' Run the iNMF (Schlicker) subtyping on a set of samples.
#'
#' @param whichSamples - a data frame (gene by samples) to be subtyped
#' @param plotHeatmaps - TRUE/FALSE, whether to create heatmap plots per subtype, default is FALSE
#' @param plotSilWidths - TRUE/FALSE, whether to create silhouette width plots per subtype, default is FALSE
#' @param directory - string, naming the directory in which to save the heatmap/silhouette width plots, default is '.'
#' @param filePrefix - string, indicating a prefix to be added to the heatmap/sil width plot names, default is ''
#' @param distance - what is the distance metric for the clustering
#' @param silhouette - Boolean; compute silhouette widths or not; default is TRUE
#' @return resInmf - clustering and silhouette widths of the samples
#' @export
#' @examples
#' print('Examples are not written yet.')

subtypeSamplesINMF = function(whichSamples, types = NULL, plotHeatmaps = FALSE, plotSilWidths = FALSE, directory='.', filePrefix='', distance = 'spearman', silhouette=TRUE){

  # first intersect signatures
  inmfSigs = intersectSignature(rownames(whichSamples), loadSchlickerSignature())

  # Perform iNMF subtyping
  # exprs, signatures, bootstrap=FALSE, randomizeFeatures=FALSE, silhouette=TRUE, plotHeatmaps=TRUE, directory=".", filePrefix="")
  resInmf = iNMF(whichSamples, types = types, signatures=inmfSigs, bootstrap = FALSE, randomizeFeatures=FALSE, silhouette=silhouette, plotHeatmaps = plotHeatmaps, directory=directory, filePrefix = filePrefix, distance = distance)

  if (plotSilWidths){
      #png(paste(directory, "/", filePrefix, "sw_step1.png", sep=""), width=4000, height=3000, res=300)
      png(paste(directory, "/", filePrefix, "sw_step1.png", sep=""), width=4000, height=3000, res=300 )
      sw = resInmf$step1$silhouette
      plotSilhouetteWidths(sw)
      dev.off()

      #png(paste(directory, "/", filePrefix, "sw_step2_c1.png", sep=""), width=4000, height=3000, res=300)
      png(paste(directory, "/", filePrefix, "sw_step2_c1.png", sep=""), width=4000, height=3000, res=300 )
      sw = resInmf$step2.c1$silhouette
      plotSilhouetteWidths(sw)
      dev.off()

      #png(paste(directory, "/", filePrefix, "sw_step2_c2.png", sep=""), width=4000, height=3000, res=300)
      png(paste(directory, "/", filePrefix, "sw_step2_c2.png", sep=""), width=4000, height=3000, res=300 )
      sw = resInmf$step2.c2$silhouette
      plotSilhouetteWidths(sw)
      dev.off()
  }

  return(resInmf)
}


#' Subtyping samples one at a time, using the Schlicker subtypes (1.1, 1.2, 1.3, 2.1, 2.2)
#'
#' Right now the implementation uses the iNMF function to assign 1 sample to a subtype.
#'
#' @param whichSamples - a data frame (gene by samples) to be subtyped
#' @param types - samples annotation vector, default is NULL
#' @param plotHeatmaps - TRUE/FALSE, whether to create heatmap plots per subtype, default is FALSE
#' @param plotSilWidths - TRUE/FALSE, whether to create silhouette width plots per subtype, default is FALSE
#' @param method - default 'maxmean': assigning cluster expression means to signature subtype expression means
#' @param directory - string, naming the directory in which to save the heatmap/silhouette width plots, default is '.'
#' @param filePrefix - string, indicating a prefix to be added to the heatmap/sil width plot names, default is ''
#' @param distance - what is the distance metric for the clustering
#' @param silhouette - Boolean, compute silhouette width or not; default is TRUE
#' @return resInmf - clustering and silhouette widths of the samples
#' @export
#' @examples
#' print('Examples not written yet')

subtypeOneAtATimeINMF = function(whichSamples, types = NULL, plotHeatmaps = FALSE, plotSilWidths = FALSE, method = 'maxmean', directory='.', filePrefix = '', distance = 'spearman', silhouette=TRUE){

  # this runs Andreas's "subtype" function but with length(samples)=1 (method = 'maxmean')
  # and combines the results
  signatures = intersectSignature(rownames(whichSamples),loadSchlickerSignature())

  # Perform iNMF subtyping
  # exprs, signatures, bootstrap=FALSE, randomizeFeatures=FALSE, silhouette=TRUE, plotHeatmaps=TRUE, directory=".", filePrefix="")
  clust = list()

  # for every sample, subtype
  for (column in 1:ncol(whichSamples)){

    if (method=='maxmean'){   # only the 'maxmean' method is implemented so far
      resInmf = iNMF(whichSamples[,column,drop=FALSE], signatures=signatures, plotHeatmaps = FALSE, bootstrap = FALSE, randomizeFeatures=FALSE, distance = distance, silhouette=silhouette)
    }

    # here, in "clust" save the step2 labels only
    if (names(resInmf$step1$clustering)=="1"){
      # take: resInmf$step2.c1$clustering
      clust[[colnames(whichSamples[,column,drop=FALSE])]] = names( resInmf$step2.c1$clustering )

    } else if (names(resInmf$step1$clustering)=="2"){
      # take: resInmf$step2.c2$clustering
      clust[[colnames(whichSamples[,column,drop=FALSE])]] = names( resInmf$step2.c2$clustering )
    }
  }

  remove(resInmf)

  resInmf = createClusteringDataStructureINMF(whichSamples, clust, distance = distance)

  if (plotSilWidths){
      png(paste(directory, "/", filePrefix, "sw_step1.png", sep=""), width=4000, height=3000, res=300)
      sw = resInmf$step1$silhouette
      plotSilhouetteWidths(sw)
      dev.off()

      png(paste(directory, "/", filePrefix, "sw_step2_c1.png", sep=""), width=4000, height=3000, res=300)
      sw = resInmf$step2.c1$silhouette
      plotSilhouetteWidths(sw)
      dev.off()

      png(paste(directory, "/", filePrefix, "sw_step2_c2.png", sep=""), width=4000, height=3000, res=300)
      sw = resInmf$step2.c2$silhouette
      plotSilhouetteWidths(sw)
      dev.off()
  }

  if (plotHeatmaps) {

      sigs = list(step1 = list(), step2.c1 = list(), step2.c2 = list())
      sigs$step1 = list('1'=signatures[['1']], '2'=signatures[['2']])
      sigs$step2.c1 = list('1.1'=signatures[['1.1']], '1.2'=signatures[['1.2']], '1.3'=signatures[['1.3']])
      sigs$step2.c2 = list('2.1'=signatures[['2.1']], '2.2'=signatures[['2.2']])


      subtype_colors = c("gray","green","gray10", "gray40", "gray70", "springgreen4","springgreen2")
    names(subtype_colors) = c('1','2','1.1','1.2','1.3','2.1','2.2')

      geneSig_colors = c("gray","green","gray10", "gray40", "gray70", "springgreen4","springgreen2")
    names(geneSig_colors) = c('1','2','1.1','1.2','1.3','2.1','2.2')

      # THIS IS HARD-CODED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      dataType_colors = c('red','gray90')
      names(dataType_colors) = c('organoids','TCGA-RNAseq')

      cellwidth = 1.2
      cellheight = 0.3
      png(paste(directory, "/", filePrefix, "step1.png", sep=""), width=4000, height=3000, res=300)
      anno_colors = list(subtypes=subtype_colors[1:2], geneSig=geneSig_colors[1:2], datatypes=dataType_colors)
      createHeatmap(whichSamples, resInmf$step1$clustering, sigs$step1, types = types, anno_colors = anno_colors, cellwidth=cellwidth, cellheight=cellheight)
      dev.off()

      cellwidth = 2.3
      cellheight = 0.7
      png(paste(directory, "/", filePrefix, "step2_c1.png", sep=""), width=4000, height=3000, res=300)
      anno_colors = list(subtypes=subtype_colors[3:5], geneSig=geneSig_colors[3:5], datatypes=dataType_colors)
      createHeatmap(whichSamples, resInmf$step2.c1$clustering, sigs$step2.c1, types = types, anno_colors = anno_colors, cellwidth=cellwidth, cellheight=cellheight)
      dev.off()

      cellwidth = 1.9
      cellheight = 0.9
      png(paste(directory, "/", filePrefix, "step2_c2.png", sep=""), width=4000, height=3000, res=300)
      anno_colors = list(subtypes=subtype_colors[6:7], geneSig=geneSig_colors[6:7], datatypes=dataType_colors)
      createHeatmap(whichSamples, resInmf$step2.c2$clustering, sigs$step2.c2, types = types, anno_colors = anno_colors, cellwidth=cellwidth, cellheight=cellheight)
      dev.off()
    }

  return(resInmf)
}


#' Convert a simple classes character vector to a clustering object.
#'
#' This produces an output similar to the iNMF() function output.
#'
#' @param whichSamples - data (columns/samples) to which the clustering applies
#' @param clust - character vector of classes, names(clust) are the samples
#' @param distance - what is the distance metric for the clustering
#' @return clustering object with silhouette widths
#' @export
#' @examples
#' print('Examples not written yet')

createClusteringDataStructureINMF = function(whichSamples, clust, distance = 'spearman'){

    # recreate the resInmf structure but for all samples (same as the output of iNMF())
    resInmf = list(step1=list(), step2.c1=list(), step2.c2=list())
    resInmf$step1$clustering = list("1"=c(),"2"=c())
    resInmf$step2.c1$clustering = list("1.1"=c(),"1.2"=c(),"1.3"=c())
    resInmf$step2.c2$clustering = list("2.1"=c(),"2.2"=c())

    for (sample in names(clust)){

        if (clust[[sample]] == "1.1"){
            resInmf$step1$clustering[["1"]] = append(resInmf$step1$clustering[["1"]], sample)
            resInmf$step2.c1$clustering[["1.1"]] = append(resInmf$step2.c1$clustering[["1.1"]], sample)

        } else if (clust[[sample]] == "1.2"){
            resInmf$step1$clustering[["1"]] = append(resInmf$step1$clustering[["1"]], sample)
            resInmf$step2.c1$clustering[["1.2"]] = append(resInmf$step2.c1$clustering[["1.2"]], sample)

        } else if (clust[[sample]] == "1.3"){
            resInmf$step1$clustering[["1"]] = append(resInmf$step1$clustering[["1"]], sample)
            resInmf$step2.c1$clustering[["1.3"]] = append(resInmf$step2.c1$clustering[["1.3"]], sample)

        } else if (clust[[sample]] == "2.1"){
            resInmf$step1$clustering[["2"]] = append(resInmf$step1$clustering[["2"]], sample)
            resInmf$step2.c2$clustering[["2.1"]] = append(resInmf$step2.c2$clustering[["2.1"]], sample)

        } else if (clust[[sample]] == "2.2"){
            resInmf$step1$clustering[["2"]] = append(resInmf$step1$clustering[["2"]], sample)
            resInmf$step2.c2$clustering[["2.2"]] = append(resInmf$step2.c2$clustering[["2.2"]], sample)
        }

    }

    # compute silhouette widths .....................
    # first for step1
    samples = append( resInmf$step1$clustering[['1']], resInmf$step1$clustering[['2']] )

    clusters = c(rep(1, length(resInmf$step1$clustering[['1']])), rep(2, length(resInmf$step1$clustering[['2']])))
    names(clusters) = samples  # this is the clustering object

    # intersect the signature ....
    signature = intersectSignature(rownames(whichSamples), loadSchlickerSignature())

    # re-do the signatures in the step1-step2.c1-step2-c2 form
    sigs = list(step1 = list('1'=signature[['1']],'2'=signature[['2']] ),
                step2.c1=list('1.1'=signature[['1.1']], '1.2'=signature[['1.2']], '1.3'=signature[['1.3']]),
                step2.c2=list('2.1'=signature[['2.1']], '2.2'=signature[['2.2']]))


    if (length(samples)<=1){
      sil = NA
    } else {
      dist.mat = 1 - cor(whichSamples[unlist(sigs$step1), samples], use="pairwise.complete.obs", method = distance)
      sil = silhouette(clusters, dist.mat)
    }

    if (sum(is.na(sil))>0){
      resInmf$step1$silhouette = NA
      resInmf$step2.c1$silhouette = NA
      resInmf$step2.c2$silhouette = NA
      return(resInmf)
    }

    rownames(sil) = samples

    for (i in 1:nrow(sil)){
        sil[i,1] = as.character(sil[i,1])
        sil[i,2] = as.character(sil[i,2])
    }
    resInmf$step1$silhouette = sil


    # then for step2.c1
    samples = append( resInmf$step2.c1$clustering[['1.1']], resInmf$step2.c1$clustering[['1.2']] )
    samples = append( samples, resInmf$step2.c1$clustering[['1.3']] )
    clusters = c(rep(11, length(resInmf$step2.c1$clustering[['1.1']])), rep(12, length(resInmf$step2.c1$clustering[['1.2']])), rep(13, length(resInmf$step2.c1$clustering[['1.3']])))

    names(clusters) = samples  # this is the clustering object

    if (length(samples)<=1){
      sil = NA
    } else {
      dist.mat = 1 - cor(whichSamples[unlist(sigs$step2.c1), samples], use="pairwise.complete.obs", method = distance)
      sil = silhouette(clusters, dist.mat)
    }

    # if there is only one cluster, sil=NA
    if (sum(is.na(sil))==0){

        # this changes the '11', '12', '13' to '1.1', '1.2', '1.3'

        rownames(sil) = samples
        for (i in 1:nrow(sil)){
            s1 = as.character(sil[i,1])
            s1 = unlist(strsplit(s1,''))
            sil[i,1] = paste(c(s1[1],'.',s1[2]),collapse='')

            s2 = as.character(sil[i,2])
            s2 = unlist(strsplit(s2,''))
            sil[i,2] = paste(c(s2[1],'.',s2[2]),collapse='')
        }

    }
    resInmf$step2.c1$silhouette = sil


    # then for step2.c2
    samples = append( resInmf$step2.c2$clustering[['2.1']], resInmf$step2.c2$clustering[['2.2']] )
    clusters = c(rep(21, length(resInmf$step2.c2$clustering[['2.1']])), rep(22, length(resInmf$step2.c2$clustering[['2.2']])))
    names(clusters) = samples  # this is the clustering object

    if (length(samples)<=1){
      sil = NA
    } else {
      dist.mat = 1 - cor(whichSamples[unlist(sigs$step2.c2), samples], use="pairwise.complete.obs", method = distance)
      sil = silhouette(clusters, dist.mat)
    }


    # if there is only one cluster, sil=NA
    if (sum(is.na(sil))==0){

        rownames(sil) = samples
        for (i in 1:nrow(sil)){
            s1 = as.character(sil[i,1])
            s1 = unlist(strsplit(s1,''))
            sil[i,1] = paste(c(s1[1],'.',s1[2]),collapse='')

            s2 = as.character(sil[i,2])
            s2 = unlist(strsplit(s2,''))
            sil[i,2] = paste(c(s2[1],'.',s2[2]),collapse='')
        }
    }
    resInmf$step2.c2$silhouette = sil

    return(resInmf)
}



#' Apply the shrunken centroids method to assign a sample to a subtype. Use the iNMF signature here.
#'
#' Ref: Tibshirani, R. (2002). Diagnosis of multiple cancer types by shrunken centroids of gene expression. PNAS, 99(10), 6567–72.
#'
#' @param df - data frame/matrix of expression values (genes by samples)
#' @param types: samples annotation (optional, character vector)
#' @param distance - what is the distance metric for the clustering
#' @return final predicted classes for the samples in df
#' @export
#' @examples
#' print('Examples not written yet')

shrunkenCentroidsINMF = function(df, types = NULL, plotHeatmap = F, directory = '.', filePrefix = '', distance = 'spearman'){

    # these are the initial cluster estimates given to PAM
    res = subtypeOneAtATimeINMF(df, distance = distance)

    # pick positive silhouette widths for training data
    classes = c()
    trainSamples = c()
    for (i in 1:nrow(res$step1$silhouette)){

        sw = as.numeric( res$step1$silhouette[i,"sil_width"] )
        if (sw < 0.1){ next }  # keep high-silhouette width samples only
        classes = append(classes, res$step1$silhouette[i,"cluster"]  )
        trainSamples = append(trainSamples, rownames(res$step1$silhouette)[i])

    }
    print(paste('shrunkenCentroidsINMF(): Step1 - Training on',length(trainSamples),'samples out of',ncol(df),sep=' '))

    require(pamr)

    dataTrain = list(x = as.matrix(df[,trainSamples]), y=classes)
    trainFit <- pamr.train(dataTrain, n.threshold = 30)

    #cvfit <- pamr.cv(trainFit, dataTrain, nfold=10, n.threshold = 30)
    ind = max( which( trainFit$errors==min(trainFit$errors) ) )

    threshold = trainFit$threshold[ind]
    print(paste('shrunkenCentroidsINMF(): Step1 - Min error threshold',threshold,sep=' '))

    p = pamr.predict(trainFit, df, threshold, type=c("class"))
    p = as.character(p)
    names(p) = colnames(df)

    print('shrunkenCentroidsINMF(): Errors in step 1:')
    pamr.confusion(trainFit, threshold)


    # do step2.c1
    classes1 = c()
    trainSamples = c()

    for (i in 1:nrow(res$step2.c1$silhouette)){

        thisSample = rownames(res$step2.c1$silhouette)[i]
        indThisSample = which(colnames(df)==thisSample)
        if (p[indThisSample]==2){ next }   # it was misclassified, will not train on it

        sw = as.numeric( res$step2.c1$silhouette[i,"sil_width"] )
        if (sw < 0.1){ next }
        classes1 = append(classes1, res$step2.c1$silhouette[i,"cluster"]  )
        trainSamples = append(trainSamples, thisSample)

    }
    if (length(trainSamples)<=1){

      trainSamples = names(p)[p==1]
      print(paste('shrunkenCentroidsINMF(): Step2.c1 - Training on',length(trainSamples),'samples out of',sum(p==1),'in cluster 1 [all samples, but very low sil.widths !!!!!!!!!!!]',sep=' '))
      classes1 = res$step2.c1$silhouette[trainSamples,"cluster"]

    } else {
      print(paste('shrunkenCentroidsINMF(): Step2.c1 - Training on',length(trainSamples),'samples out of',sum(p==1),'in cluster 1',sep=' '))
    }

    if (length(levels(as.factor(classes1)))==1){
      print('shrunkenCentroidsINMF(): Cannot train on a single class variable, returning
                                      the initialization class membership.')
      return(res)
    }

    dataTrain = list(x = as.matrix(df[,trainSamples]), y=classes1)
    trainFit <- pamr.train(dataTrain, n.threshold = 30)

    #cvfit <- pamr.cv(trainFit, dataTrain, nfold=10, n.threshold = 30)
    ind = max( which( trainFit$errors==min(trainFit$errors) ) )

    threshold = trainFit$threshold[ind]
    print(paste('shrunkenCentroidsINMF(): Step2.c1 - Min error threshold',threshold,sep=' '))

    p1 = pamr.predict(trainFit, df[,p=="1"], threshold, type=c("class"))
    p1 = as.character(p1)
    names(p1) = colnames(df[,p=="1"])

    print('shrunkenCentroidsINMF(): Errors in step2.c1:')
    pamr.confusion(trainFit, threshold)


    # do step2.c2
    classes2 = c()
    trainSamples = c()
    for (i in 1:nrow(res$step2.c2$silhouette)){

        thisSample = rownames(res$step2.c2$silhouette)[i]
        indThisSample = which(colnames(df)==thisSample)
        if (p[indThisSample]==1){ next }   # it was misclassified, will not train on it

        sw = as.numeric( res$step2.c2$silhouette[i,"sil_width"] )
        if (sw < 0.1){ next }
        classes2 = append(classes2, res$step2.c2$silhouette[i,"cluster"]  )
        trainSamples = append(trainSamples, thisSample)

    }
    if (length(trainSamples)<=1){
      trainSamples = names(p)[p==2]
      print(paste('shrunkenCentroidsINMF(): Step2.c2 - Training on',length(trainSamples),'samples out of',sum(p==2),'in cluster 2 [all samples, but very low sil.widths !!!!!!!!!!!]',sep=' '))
      classes2 = res$step2.c2$silhouette[trainSamples,"cluster"]
    } else {
      print(paste('shrunkenCentroidsINMF(): Step2.c2 - Training on',length(trainSamples),'samples out of',sum(p==2),'in cluster 2',sep=' '))

    }

    if (length(levels(as.factor(classes2)))==1){
      print('shrunkenCentroidsINMF(): Cannot train on a single class variable, returning
                                      the initialization class membership.')
      return(res)
    }

    dataTrain = list(x = as.matrix(df[,trainSamples]), y=classes2)
    trainFit <- pamr.train(dataTrain, n.threshold = 30)

    #cvfit <- pamr.cv(trainFit, dataTrain, nfold=10, n.threshold = 30)
    ind = max( which( trainFit$errors==min(trainFit$errors) ) )

    threshold = trainFit$threshold[ind]
    print(paste('shrunkenCentroidsINMF(): Step2.c2 - Min error threshold',threshold,sep=' '))

    p2 = pamr.predict(trainFit, df[,p=="2"], threshold, type=c("class"))
    p2 = as.character(p2)
    names(p2) = colnames(df[,p=="2"])

    print('shrunkenCentroidsINMF(): Errors in step2.c2:')
    pamr.confusion(trainFit, threshold)


    # export the resulting clustering
    clust = list()
    for (i in 1:ncol(df)){
        # find in p:
        ind = which(colnames(df)[i] == names(p))
        #print( colnames(df)[i] )
        cluster = p[ind]
        #print(cluster)
        clust[[ cluster ]] = append(clust[[ cluster ]], colnames(df)[i])
    }
    clust1 = list()
    for (i in 1:ncol(df)){
        # find in p1:
        ind = which(colnames(df)[i] == names(p1))
        if (length(ind) == 0){ next }
        cluster = p1[ind]
        clust1[[ cluster ]] = append(clust1[[ cluster ]], colnames(df)[i])
    }
    clust2 = list()
    for (i in 1:ncol(df)){
        # find in p2:
        ind = which(colnames(df)[i] == names(p2))
        if (length(ind) == 0){ next }
        cluster = p2[ind]
        clust2[[ cluster ]] = append(clust2[[ cluster ]], colnames(df)[i])
    }


    if (plotHeatmap){

        # defaults for the Cell paper:  cellwidth = 1, cellheight = 0.5
        cellwidth = 1
        cellheight = 0.5

        subtype_colors = c("gray","green","gray10", "gray40", "gray70", "springgreen4","springgreen2")
        names(subtype_colors) = c('1','2','1.1','1.2','1.3','2.1','2.2')

        geneSig_colors = c("gray","green","gray10", "gray40", "gray70", "springgreen4","springgreen2")
        names(geneSig_colors) = c('1','2','1.1','1.2','1.3','2.1','2.2')

        #dataType_colors = c('orange','lightgreen','skyblue')
        #names(dataType_colors) = c('organoids','TCGA-Agilent','TCGA-RNAseq')
        dataType_colors = c('orange','skyblue')
        names(dataType_colors) = c('organoids','TCGA-RNAseq')

        # step1
        png(paste(directory, "/", filePrefix, "step1.png", sep=""), width=4000, height=3000, res=300)
        anno_colors = list(subtypes = subtype_colors[1:2], geneSig = subtype_colors[1:2], datatypes = dataType_colors)

        intsig = intersectSignature(rownames(df), loadSchlickerSignature())
        intsig = list('1'=intsig[['1']], '2'=intsig[['2']])
        createHeatmap(df, clust, intsig, types = types, cellwidth=cellwidth, cellheight = cellheight, anno_colors = anno_colors)
        dev.off()

        # step2.c1
        png(paste(directory, "/", filePrefix, "step2_c1.png", sep=""), width=4000, height=3000, res=300)
        anno_colors = list(subtypes = subtype_colors[3:5], geneSig = subtype_colors[3:5], datatypes = dataType_colors)

        intsig = intersectSignature(rownames(df), loadSchlickerSignature())
        intsig = list('1.1'=intsig[['1.1']], '1.2'=intsig[['1.2']], '1.3'=intsig[['1.3']])
        createHeatmap(df, clust1, intsig, types = types, cellwidth=cellwidth, cellheight = cellheight, anno_colors = anno_colors)
        dev.off()

        # step2.c2
        png(paste(directory, "/", filePrefix, "step2_c2.png", sep=""), width=4000, height=3000, res=300)
        anno_colors = list(subtypes = subtype_colors[6:7], geneSig = subtype_colors[6:7], datatypes = dataType_colors)

        intsig = intersectSignature(rownames(df), loadSchlickerSignature())
        intsig = list('2.1'=intsig[['2.1']], '2.2'=intsig[['2.2']])
        createHeatmap(df, clust2, intsig, types = types, cellwidth=cellwidth, cellheight = cellheight, anno_colors = anno_colors)
        dev.off()


    }

    p = append(p1, p2)
    res = createClusteringDataStructureINMF(df, p, distance = distance)
    return(res)

}


# ................................................
# functions from Andreas's github repository below
# https://github.com/andreas-schlicker/crcsc/tree/master/groups/F/pipeline
# ................................................


#' Map sample clusters to the correct subtype.
#' Computes the average expression value for each sample cluster and each signature.
#' Then iteratively assigns subtype IDs to the clusters starting with the cluster
#' with the highest certainty for a given subtype.
#' The function assumes that all samples and features are contained in the expression
#' matrix. Missing values will be ignored.
#'
#' @param samples named list with sample clusters
#' @param named list of signatures
#' @param exprs expresssion matrix with samples in columns and features in rows
#' @return the mapping between cluster ID and subtype label
#' @export
#' @author Andreas Schlicker
assignClusterId = function(samples, signatures, exprs) {

    # Calculated mean expression for each sample cluster and each signature
    means = sapply(samples, function(samp) { sapply(signatures, function(x) { mean(exprs[x, samp], na.rm=TRUE) }) }, simplify=TRUE)

    # Calculate the margin for each cluster assignment
    margin = apply(means, 1, function(x) { sort(x, decreasing=TRUE)[1] - sort(x, decreasing=TRUE)[2] })

    # Mapping between subtype label and cluster ID
    mapping = rep("", times=length(names(signatures)))
    names(mapping) = names(signatures)
    # All cluster IDs
    clusts = colnames(means)

    # Cycle through all subtypes in order of decreasing margin
    for (n in names(sort(margin, decreasing=TRUE))) {
        if (length(clusts) > 1) {
            # There are more than one cluster left
            # Take one of the subtypes with the largest means
            mapping[n] = names(which(means[n, clusts] == max(means[n, clusts])))[1]
            clusts = setdiff(clusts, mapping[n])
        } else {
            # Assign the remaining subtype
            mapping[n] = clusts
        }
    }
    # Reverse the mapping between subtype and cluster ID
    res = names(mapping)
    names(res) = mapping

    res
}

#' Performs one step of the iNMF clustering. Essentially, this is a
#' hierarchical clustering of the given samples using the given signatures.
#' Cluster IDs are assigned using function "assignClusterId" but should be
#' manually verified using a heatmap or similar. The number of clusters is
#' equal to the number of signatures. Also computes silhouette values for the
#' clustering.
#'
#' @param exprs expression matrix with samples in columns and features in rows
#' @param signatures list of signatures to use.
#' @param samples vector with sample names. If NULL (default), all samples will be
#' clustered.
#' @param silhouette boolean value indicating whether the silhouette is to be
#' computeted
#' @return Named list with two elements. The first element "clustering" contains a
#' named list with the assignment of samples to clusters. The names correspond to the
#' names used for the signatures. The second element "silhouette" contains an object
#' of class silhouette, if silhouette==TRUE.
#' @export
#' @author Andreas Schlicker

subtype = function(exprs, signatures, samples=NULL, silhouette=TRUE, distance = 'spearman') {
    if (silhouette) {
        require(cluster) || stop("I need package \"cluster\" to calculate the silhouette!")
    }

    if (is.null(samples)) {
        # if "samples" not specified, use all samples in "exprs"
        samples = colnames(exprs)
    }

    if (length(samples) == 1) {
        # single-sample subtyping:
        # assign to the highest expression-mean subset of genes from the signature
        means = unlist(lapply(signatures, function(x) { mean(exprs[x, samples]) }))
        res = list()
        res[names(means)[which(means == max(means))][1]] = c(samples)
        return(list(clustering=res))
    }

    # compute 1-correlation distance between the samples, using genes from the signature only
    dist.mat = 1 - cor(exprs[unlist(signatures), samples], use="pairwise.complete.obs", method = distance)

    # cut the hierarchical clustering tree at either the number of clusters from the signatures variable,
    #                   or the number of samples if length(samples)<length(signatures)

    # "clustering" is a sample:cluster type of list
    clustering = cutree(hclust(as.dist(dist.mat), method="complete"), k=min(length(signatures), length(samples)))


    # flip the list assignment: in "clust.list" the names are the values in "clustering", and vice versa
    # "clust.list" is a cluster:samples type of list
    # these clusters are simply hierarchical clustering "clusters" - they are not mapped to the clusters from the signatures yet
    clust.list = list()
    for (i in unique(clustering)) {
        clust.list[[as.character(i)]] = names(clustering[clustering == i])
    }


    # in assignClusterId() the clusters from the hierarchical clustering above get mapped to the signature clusters
    # this is done by sequential maximum [gene expression] mean matching
    idMap = assignClusterId(clust.list, signatures, exprs)

    # here the "clust.list" is re-written as "res", where the hierarchical clustering id's are named with the assigned "idMap" id's from above
    res = list()
    for (n in names(clust.list)) {
        res[[idMap[n]]] = clust.list[[n]]
    }

    # compute silhouette widths for the hierarchical clustering above
    # rename clusters with the signature cluster labels from idMap
    sil = NULL
    if (silhouette) {
        sil = silhouette(clustering, dist.mat)
    }
    if (silhouette && !is.na(sil[1])){
        for (i in names(idMap)) {
            sil[which(sil[, "cluster"] == i), "cluster"] = idMap[[i]]
            sil[which(sil[, "neighbor"] == i), "neighbor"] = idMap[[i]]
        }
        rownames(sil) = samples
    }

    if (is.null(sil)) {
        list(clustering=res)
    } else {
        list(clustering=res, silhouette=sil)
    }
}


#' Run iNMF subtyping of the specified data set. This functions rests primarily on the subtype() function.
#'
#' @param exprs expression matrix with samples in columns and features in rows
#' @param signatures named list with subtyping signatures mapped to feature ids
#'        used in the expression matrix. The following names have to be used "step1",
#'        "step2.c1" and "step2.c2". The function "signaturesList" can be used to
#'        obtain this list.
#' @param types: an annotation vector of size 1xncol(exprs)
#' @param bootstrap boolean indicating whether running in bootstrap mode. If TRUE, samples will be randomly selected
#' @param randomizeFeatures boolean indicating if features should be randomized before running the subtyping; default: FALSE
#' @param silhouette boolean indicating whether the silhouette should be calculated
#' @param plotHeatmaps boolean indicating whether clustering heatmaps should be saved in files
#' @param directory path to the directory where heatmaps will be saved
#' @param filePrefix file name prefix for the heatmap files
#' @param distance - what is the distance metric for the clustering
#' @return a named list with the clustering results
#' @export
#' @author Andreas Schlicker, modified by Gergana Bounova

iNMF = function(exprs, signatures, types = NULL, bootstrap=FALSE, randomizeFeatures=FALSE, silhouette=TRUE, plotHeatmaps=TRUE, directory=".", filePrefix="", distance = 'spearman') {

    if (bootstrap) {
        samples = sample(colnames(exprs), ncol(exprs), TRUE)
    } else {
        samples = colnames(exprs)
    }

    if (randomizeFeatures) {
        rownames(exprs) = sample(rownames(exprs))
    }

    sigs = list(step1 = list('1'=signatures[['1']],'2'=signatures[['2']] ),
                step2.c1=list('1.1'=signatures[['1.1']], '1.2'=signatures[['1.2']], '1.3'=signatures[['1.3']]),
                step2.c2=list('2.1'=signatures[['2.1']], '2.2'=signatures[['2.2']]))


    clust1 = subtype(exprs, sigs$step1, samples, silhouette, distance = distance)
    clust2.c1 = subtype(exprs, sigs$step2.c1, clust1$clustering$'1', silhouette, distance = distance)
    clust2.c2 = subtype(exprs, sigs$step2.c2, clust1$clustering$'2', silhouette, distance = distance)

    subtype_colors = c("gray","green","gray10", "gray40", "gray70", "springgreen4","springgreen2")
    names(subtype_colors) = c('1','2','1.1','1.2','1.3','2.1','2.2')

    geneSig_colors = c("gray","green","gray10", "gray40", "gray70", "springgreen4","springgreen2")
    names(geneSig_colors) = c('1','2','1.1','1.2','1.3','2.1','2.2')

    # THIS IS HARD-CODED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    dataType_colors = c('orange','skyblue')
    names(dataType_colors) = c('organoids','TCGA-RNAseq')

    if (is.null(types)){ dataType_colors = NULL }

    if (plotHeatmaps){
        cellwidth = 1.1
        cellheight = 0.2
        png(paste(directory, "/", filePrefix, "step1.png", sep=""), width=4000, height=3000, res=300)
        anno_colors = list(subtypes=subtype_colors[1:2], geneSig=geneSig_colors[1:2], datatypes=dataType_colors)
        createHeatmap(exprs, clust1$clustering, sigs$step1, types = types, anno_colors = anno_colors, cellwidth=cellwidth, cellheight=cellheight)
        dev.off()

        cellwidth = 2.3
        cellheight = 0.7
        png(paste(directory, "/", filePrefix, "step2_c1.png", sep=""), width=4000, height=3000, res=300)
        anno_colors = list(subtypes=subtype_colors[3:5], geneSig=geneSig_colors[3:5], datatypes=dataType_colors)
        createHeatmap(exprs, clust2.c1$clustering, sigs$step2.c1, types = types, anno_colors = anno_colors, cellwidth=cellwidth, cellheight=cellheight)
        dev.off()

        cellwidth = 1.9
        cellheight = 0.9
        png(paste(directory, "/", filePrefix, "step2_c2.png", sep=""), width=4000, height=3000, res=300)
        anno_colors = list(subtypes=subtype_colors[6:7], geneSig=geneSig_colors[6:7], datatypes=dataType_colors)
        createHeatmap(exprs, clust2.c2$clustering, sigs$step2.c2, types = types, anno_colors = anno_colors, cellwidth=cellwidth, cellheight=cellheight)
        dev.off()
    }

    list(step1=clust1, step2.c1=clust2.c1, step2.c2=clust2.c2)
}


#' Computes the number of mis-classified samples between two clusterings.
#'
#' This is implemented for results of the iNMF algorithm, and could be ran for step1 or step2.
#'
#' @param resInmf1 - clustering object for the first clustering
#' @param resInmf2 - clustering object for the second clustering
#' @param samplesOfInterest - out of the samples in the data, the subset of samples for which the comparison is performed
#' @param whichStep - which step in the clustering is being compared, 'step1' or 'step2', default is 'step2'
#' @return list -
#' @export
#' @examples
#' print('Examples not written yet')

compareTwoClusteringsINMF = function(resInmf1, resInmf2, samplesOfInterest, whichStep = 'step2'){

  c1 = resInmf1$step1$silhouette
  if (whichStep == 'step2'){
    c1 = resInmf1$step2.c1$silhouette
    c2 = resInmf1$step2.c2$silhouette
  }

  clusterings = c1[,1]
  if (whichStep == 'step2'){ clusterings = append(clusterings, c2[,1]) }

  silWidths = c1[,3]
  if (whichStep == 'step2'){ silWidths = append(silWidths, c2[,3]) }

  # make sure samples are ordered the same way
  stopifnot( sum( sort(names(clusterings)) != sort(names(silWidths)) )==0 )
  firstSet = cbind(clusterings[samplesOfInterest], silWidths[samplesOfInterest])
  colnames(firstSet) = c('cluster','sil_width')


  # now prepare the second set of samples

  c1 = resInmf2$step1$silhouette
  if (whichStep == 'step2'){
    c1 = resInmf2$step2.c1$silhouette
    c2 = resInmf2$step2.c2$silhouette
  }

  clusterings = c1[,1]
  if (whichStep == 'step2'){ clusterings = append(clusterings, c2[,1]) }

  silWidths = c1[,3]
  if (whichStep == 'step2'){ silWidths = append(silWidths, c2[,3]) }

  stopifnot( sum( sort(names(clusterings)) != sort(names(silWidths)) )==0 )
  secondSet = cbind(clusterings[samplesOfInterest], silWidths[samplesOfInterest])
  colnames(secondSet) = c('cluster','sil_width')


  compare = cbind( firstSet[samplesOfInterest,], secondSet[samplesOfInterest,]  )
  isitmatching = compare[,1] == compare[,3]
  print(paste('compareTwoClusteringsINMF(): Mis-classified samples',sum(isitmatching==FALSE),'out of',length(isitmatching),sep=' '))

  compare = cbind( compare[samplesOfInterest,], isitmatching[samplesOfInterest]  )  # add clustering matches/mismatches

  return(list(compare, sum(isitmatching==FALSE)/length(isitmatching)))
}


#' Convert the INMF step1, step2.c1, step2.c2 structure to lists of clusters.
#'
#' @param resInmf - [list] =, result of iNMF()
#' @return simple list of clusters
#' @export
#' @examples
#' print('Examples not written yet')

unpackINMFclustering = function(resInmf){

    clust = list('1.1'=c(), '1.2'=c(), '1.3'=c(), '2.1'=c(), '2.2'=c())

    clust[['1.1']] = resInmf$step2.c1$clustering[['1.1']]
    clust[['1.2']] = resInmf$step2.c1$clustering[['1.2']]
    clust[['1.3']] = resInmf$step2.c1$clustering[['1.3']]
    clust[['2.1']] = resInmf$step2.c2$clustering[['2.1']]
    clust[['2.2']] = resInmf$step2.c2$clustering[['2.2']]

    return(clust)
}


#' Convert list of clusters to a character vector
#'
#' @param clust - list of clusters
#' @return a character vector representation of the clusters (names are cluster members)
#' @export
#' @examples
#' print('Examples not written yet')

clust2list = function(clust){

  L = c()
  names = c()
  for ( c in names(clust) ){
    for ( member in clust[[c]] ){
      L = append(L, c)
      names = append(names, member)
    }
  }

  names(L) = names
  return(L)

}


#' Unpack the INMF step1, step2.c1, step2.c2 structure into a dataframe with assignment, neighbor and sil width
#'
#' @param resInmf - [list] =, result of iNMF()
#' @return data frame with cluster assignment, neighboring cluster and sil width
#' @export
#' @examples
#' print('Examples not written yet')

unpackINMF = function(resInmf){

    df = data.frame(stringsAsFactors = FALSE)

    if ( length(resInmf$step2.c1$silhouette)==1 & is.na(resInmf$step2.c1$silhouette) ){

      cl = resInmf$step2.c1$clustering
      l = clust2list(cl)

      mat = l
      mat = cbind(mat, rep(NA, length(l)))
      mat = cbind(mat, rep(NA, length(l)))

      rownames(mat) = names(l)
      colnames(mat) = c('cluster', 'neighbor', 'sil_width')

      df = rbind(df, as.matrix(mat))

    } else {

      df = rbind(df, as.matrix( resInmf$step2.c1$silhouette[, c('cluster', 'neighbor', 'sil_width')] ) )

    }

    if ( is.na(resInmf$step2.c2$silhouette) ){

      cl = resInmf$step2.c2$clustering
      l = clust2list(cl)

      mat = l
      mat = cbind(mat, rep(NA, length(l)))
      mat = cbind(mat, rep(NA, length(l)))

      rownames(mat) = names(l)
      colnames(mat) = c('cluster', 'neighbor', 'sil_width')


      df = rbind(df, as.matrix(mat))


    } else {

      df = rbind(df, as.matrix( resInmf$step2.c2$silhouette[, c('cluster', 'neighbor', 'sil_width')] ) )
    }


    return(df)
}
