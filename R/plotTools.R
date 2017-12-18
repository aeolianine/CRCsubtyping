source('~/tools/generalPlottingTools.R')

#' Plotting the scores from the signature table (CRCassigner-786) as a heatmap
#'
#' Reference: Sadanandam et al (2013). A colorectal cancer classification system that associates cellular phenotype and responses to therapy. Nature Medicine, 19(5), 619–25. doi:10.1038/nm.3175
#'
#' @param no inputs needed
#' @seealso \code{\link{loadSadanandamSignature}} - this is a function that actually loads the signature as a list
#' @export
#' @examples
#' plotSadanandamSignatureScores()

plotSadanandamSubtypeCentroids = function(){

  data(loadSadanandamSignature)

  require(gplots)
  require(RColorBrewer)
  palette = colorRampPalette(c('blue','white'))(n=1000)
  heatmap.2(as.matrix(t), dendrogram = c('none'), scale=c('row'), col = palette, cexRow = 0.1, las=1, cexCol = 1)

  # attempt to plot column labels on top, didn't work
  #axis(3, 1:ncol(t), labels = colnames(t), las = 2, tick = 0, cex.axis = 1)
}



#' Plot two mds components for results summarized in a clustering object.
#'
#' Uses the iNMF type results and cmdscale, as well as "type" annotation; and the Schlicker signature.
#'
#' @param exprs - expression matrix
#' @param resInmf - a clustering object, the output of iNMF()
#' @param types - indicating different sample types (factor of level 2)
#' @param method - can be 'cmdscale', 'princomp', 'isoMDS', or 'tsne'; default is 'cmdscale'
#' @param plotName, if the plot is being saved; default is NULL
#' @param title - character string for the title of the plot, if any; default is ''
#' @export
#' @examples
#' org = loadPilotOrganoidExpression(meanCenter=TRUE)
#' res = subtypeSampleINMF(org)
#' types = rep('t',ncol(org))
#' types[grepl('n',colnames(org))] = 'n'
#' plotClusterMDSSchlicker(org, res, types)

plotClusterMDSSchlicker = function(exprs, resInmf, types, method = 'cmdscale', plotName = NULL, title = ''){

    if (grepl('cmdscale',method)){
      k = cmdscale(1-cor(exprs),k=2)
    } else if (grepl('princomp',method)){
      k = princomp(1-cor(exprs), scores=T, scale=c('none'))
      k = k$loadings[,1:2]
    } else if (grepl('isoMDS',method)){
      k = isoMDS(1-cor(exprs), y = cmdscale(1-cor(exprs), 2), k = 2, maxit = 50, trace = TRUE)
      k = k$points
    } else if (grepl('tsne',method)){
      require(tsne)
      d = 1-cor(exprs)   # 1-correlation distance
      k = tsne(d, k=2, perplexity=50, max_iter = 2000)
    }

    rownames(k) = colnames(exprs)

    levels = levels(as.factor(types))

    if (!is.null(plotName)){ png(plotName, res=200, width=2000, height=2000) }
    plot(k[,1],k[,2], col='black', pch=19, cex=0.2, main = title)

    points(k[resInmf$step2.c2$clustering$'2.2',1],k[resInmf$step2.c2$clustering$'2.2',2], col='blue', cex=0.5, pch=19)
    points(k[resInmf$step2.c2$clustering$'2.1',1],k[resInmf$step2.c2$clustering$'2.1',2], col='skyblue', cex=0.5, pch=19)
    points(k[resInmf$step2.c1$clustering$'1.2',1],k[resInmf$step2.c1$clustering$'1.2',2], col='darkred', cex=0.5, pch=19)
    points(k[resInmf$step2.c1$clustering$'1.1',1],k[resInmf$step2.c1$clustering$'1.1',2], col='pink', cex=0.5, pch=19)
    points(k[resInmf$step2.c1$clustering$'1.3',1],k[resInmf$step2.c1$clustering$'1.3',2], col='orange', cex=0.5, pch=19)

    #for (i in 1:length(levels)){
    for (i in 1:2){
      points(k[types==levels[i],1],k[types==levels[i],2], col='black', pch=20+i, cex=1.1)
    }

    legend('bottomleft',c('2.2','2.1','1.2','1.1','1.3'), col = c('blue','skyblue','darkred','pink','orange'), lty=c(1,1,1,1,1), bty = 'n' )

    if (!is.null(plotName)){ dev.off() }

}



#' Plot two mds components for results summarized in a clustering object.
#'
#' Uses the Sadanandam type results and cmdscale, as well as "type" annotation; and the Sadanandam signature.
#'
#' @param exprs - expression matrix
#' @param resInmf - a clustering object, the output of iNMF()
#' @param types - indicating different sample types (factor of level 2)
#' @param method - can be 'cmdscale', 'princomp', 'isoMDS', or 'tsne'; default is 'cmdscale'
#' @param plotName, if the plot is being saved; default is NULL
#' @param title - character string for the title of the plot, if any; default is ''
#' @export
#' @examples
#' ...

plotClusterMDSSadanandam = function(exprs, clust, types, method = 'cmdscale', plotName = NULL, title = ''){

    if (grepl('cmdscale',method)){
      k = cmdscale(1-cor(exprs),k=2)
    } else if (grepl('princomp',method)){
      k = princomp(1-cor(exprs), scores=T, scale=c('none'))
      k = k$loadings[,1:2]
    } else if (grepl('isoMDS',method)){
      k = isoMDS(1-cor(exprs), y = cmdscale(1-cor(exprs), 2), k = 2, maxit = 50, trace = TRUE)
      k = k$points
    } else if (grepl('tsne',method)){
      require(tsne)
      d = 1-cor(exprs)   # 1-correlation distance
      k = tsne(d, k=2, perplexity=50, max_iter = 2000)
    }

    rownames(k) = colnames(exprs)

    levels = levels(as.factor(types))

    if (!is.null(plotName)){ png(plotName, res=200, width=2000, height=2000) }
    plot(k[,1],k[,2], col='black', pch=19, cex=0.2, main = title)


    points(k[clust$clustering$Enterocyte,1],k[clust$clustering$Enterocyte,2], col='blue', cex=0.5, pch=19)
    points(k[clust$clustering$TA,1],k[clust$clustering$TA,2], col='skyblue', cex=0.5, pch=19)
    points(k[clust$clustering$Goblet.like,1],k[clust$clustering$Goblet.like,2], col='darkred', cex=0.5, pch=19)
    points(k[clust$clustering$Stem.like,1],k[clust$clustering$Stem.like,2], col='pink', cex=0.5, pch=19)
    points(k[clust$clustering$Inflammatory,1],k[clust$clustering$Inflammatory,2], col='orange', cex=0.5, pch=19)

    for (i in 1:2){
      points(k[types==levels[i],1],k[types==levels[i],2], col='black', pch=20+i, cex=1.1)
    }

    legend('topright',c('Enterocyte','TA','Goblet.like','Stem.like','Inflammatory'), col = c('blue','skyblue','darkred','pink','orange'), lty=c(1,1,1,1,1), bty='n')
    if (!is.null(plotName)){ dev.off() }

}


#' Plot silhouette widths.
#'
#' @param sw, silhouette width object
#' @export
#' @examples
#' ...

plotSilhouetteWidths = function(sw){


  clusters = sort(unique(sw[,'cluster']))
  numClusters = length(clusters)
  absoluteMax = max(as.numeric(sw[,'sil_width']))   # maximum silhouette width
  par(mfrow = c(1,numClusters))

  for (c in 1:numClusters){
    swc = sw[sw[,'cluster']==clusters[c],]
    if (length(dim(swc))==0){   # only one member in this cluster

      par(lab=c(1,1,1))
      label = names(sw[,'cluster'])[sw[,'cluster']==clusters[c]]
      value = as.numeric( sw[sw[,'cluster']==clusters[c],'sil_width'] )
      barplot(value, width = 1, axes = FALSE, ylim = c(0,absoluteMax+0.01*value), main = clusters[c])
      axis(side=1,at=c(0.5),labels = c(label), cex.lab = 0.8)
      axis(side=2, at=c(0,value), las=1)
      next

      }

    swc = swc[order(swc[,'sil_width'],decreasing=T),]

    labels = names(swc[,'sil_width'])
    par(lab=c(length(labels),1,length(labels)))

    topSilW = max(as.numeric(swc[,'sil_width']))
    barplot(as.numeric(swc[,'sil_width']), width = 0.8, axes = FALSE, ylim = c(0,absoluteMax+0.01*topSilW), main = clusters[c])
    ticksWhere = c(0, topSilW/2, topSilW)
    ticksWhere = round(ticksWhere, digits = 3)
    axis(side=1,at=1:length(labels)-0.5,labels = labels, cex.lab = 0.05, las = 2)
    axis(side=2, at = ticksWhere, las = 1)

  }
  par(mfrow=c(1,1))
}


#' Create subtypes heatmap.
#'
#' Adapted from Andreas's createHeatmap() function. Need to add examples and explain inputs. Uses the function aheatmap() from the NMF package.
#'
#' @param exprs - expression matrix (genes by samples)
#' @param clustering - list of sample clusters
#' @param signatures - list of genes per subtype
#' @param types - factor for samples types, for example c('organoids', 'TCGA-RNAseq'), used as column annotation in aheatmap()
#' @param anno_colors - used for annotation colors in aheatmap(), default is NULL
#' @param filename - to be used within aheatmap(), default is NA
#' @param cellwidth/cellheight - inputs for aheatmap(), to control the size of cells plotted
#' @seealso \code{\link[NMF]{aheatmap}}
#' @export
#' @examples
#' createHeatmap(mat,clust,intsig, types = src, filename = '~/projects/CRCorganoids/figures/heatmap_organoids_rnaSeq_inmf_complete.png')

createHeatmap = function(exprs, clustering, signatures, types = NULL, anno_colors = NULL, filename = NA, cellwidth = 1, cellheight = 1) {

  clustering = clustering[order(names(clustering))]
  signatures = signatures[order(names(signatures))]

  require(gplots) || stop("I need package \"gplots\" for doing this.")
  require(NMF) || stop("I need package \"NMF\" for doing this.")

  bounds = quantile(exprs[unlist(signatures), unlist(clustering)], probs=c(0.01, 0.99), na.rm=TRUE)

  subtypes = c()

  for (i in 1:length(clustering)) {
    subtypes = c(subtypes, rep(names(clustering)[i], times=length(clustering[[i]])) )
  }

  features = c()
  geneSig = c()
  for (i in 1:length(signatures)) {
      geneSig = c(geneSig, rep(names(signatures)[i], times=length(signatures[[i]] ) ) )
      features = c(features, signatures[[i]])
  }


  #heatmap.2(exprs[features, samples],
  #          trace="none",
  #          scale=c("none"),
  #          col=colorpanel(49, low="blue", high="yellow"),
  #          breaks=seq(bounds[1], bounds[2], length.out=50),
  #          ColSideColors=csd,
  #          RowSideColors=rsd,
  #          Colv=NA,
  #          Rowv=NA,
  #          dendrogram="none")


  # calculate a reordering, by clustering each subtype separately
  ordering = c()
  for (i in 1:length(clustering)){

      clusterName = names(clustering)[i]
      sig = signatures[[clusterName]]
      score = c()
      for (sample in clustering[[i]]){
          score = append(score, mean(exprs[sig,sample]))
      }

      names(score) = clustering[[i]]
      ordering = append(ordering, names(sort(score)))

      #if (length(clustering[[i]])==0){ next }
      #if (length(clustering[[i]])==1) { ordering = append(ordering, clustering[[i]]); next }
      ## if euclidean
      ##dmat = dist(t(exprs[features, clustering[[i]]]), method = 'euclidean')
      ##d = hclust(dmat, method = 'complete')
      #d = hclust(as.dist(1-cor(exprs[features, clustering[[i]]] ) ), method = 'complete')
      #ordering = append(ordering, d$labels[d$order])
  }
  datatypes = c()
  for (sample in ordering){
      datatypes = append(datatypes, types[ colnames(exprs)==sample ] )
  }

  aheatmap(exprs[features, ordering, drop=FALSE],
           color = colorpanel(49, low="blue", high="yellow"),
           scale='none',
           breaks=seq(bounds[1], bounds[2], length.out=50),
           labRow = NULL,
           labCol = NULL,
           Rowv=NA,
           Colv=NA,
           cexRow = 0,
           cexCol = 0,
           annCol = cbind(datatypes, subtypes),
           annRow = cbind(NULL,geneSig),
           annColors=anno_colors,
           #width = 5,
           #height = 4,
           cellwidth = cellwidth,
           cellheight = cellheight,  # used these to send figure to Marc van de Wetering
           fontsize = 15,
           filename = NA
           )
  #dev.off()
}
