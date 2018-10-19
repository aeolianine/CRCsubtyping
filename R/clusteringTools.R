#' Implement W. Rand's metric for comparing two clusters.
#'
#' Source: Rand, W. M. (1971). Objective Criteria for the Evaluation of Clustering Methods. Journal of the American Statistical Association, 66(336), 846–850.
#'
#' @param numItems - total number of items being clustered
#' @param clustering1 - list of clusters in clustering 1; each cluster is a vector of items
#' @param clustering2 - list of clusters in clustering 2; each cluster is a vector of items
#' @return Rand's metric as defined in his paper.
#' @export
#' @examples
#' c = randClusteringMetric(5, list('1'=c(1,2,3), '2'=c(4,5)), list('1'=c(1), '2'=c(2,3,4), '3'=c(5)) )
#' print(c)

randClusteringMetric = function(numItems, clustering1, clustering2){

    nij = c()
    for (c1 in clustering1){
        for (c2 in clustering2){
            nij = append( nij, length(intersect(c1,c2)) )
        }
    }
    nij = matrix(nij, byrow = TRUE, nrow = length(clustering1), ncol = length(clustering2))

    sumnijsq = 0
    sumisumjnijsq = 0
    for (i in 1:length(clustering1)){
        sumisumjnijsq = sumisumjnijsq + sum(nij[i,])^2
        sumnijsq = sumnijsq + sum(nij[i,]^2)
    }

    sumjsuminijsq = 0
    for (j in 1:length(clustering2)){
        sumjsuminijsq = sumjsuminijsq + sum(nij[,j])^2
    }

    c12 = numItems*(numItems-1)/2 - 0.5*(sumisumjnijsq + sumjsuminijsq) + sumnijsq
    c12 = c12/( numItems*(numItems-1)/2 )

    return(c12)
}


#' Calculates the average silhouette width of a clustering.
#'
#' @param clustering - a clustering object, such as the output of iNMF()
#' @return avesw - a number: the average silhouette width
#' @export
#' @examples
#' res = subtypeOneAtATimeINMF(data)
#' avesw = averageSilhouetteWidth(res)

averageSilhouetteWidth = function(sw){

    avesw = mean( as.numeric( sw[,'sil_width'] ) )
    return(avesw)

}


#' Distance function for two vectors. Auxiliary function.
#'
#' @param x - first vector
#' @param y - second vector
#' @param metric - distance metric, could be 'cor' for 1-pearson correlation or 'eucl' for Euclidean distance
#' @return the distance between the two vector [numeric, length=1]
#' @export
#' @examples
#' ....
distanceFunction = function(x,y,metric='spearman'){
       stopifnot( metric %in% c('spearman','euclidean','pearson') )
       stopifnot( length(x) == length(y) )
       stopifnot(class(x)=='numeric')
       stopifnot(class(y)=='numeric')

       if (metric!='euclidean'){
           return( 1-cor(x,y,method=metric) )
       } else if (metric=='euclidean'){
           return( dist(rbind(x,y), method=metric) )
       }
   }


#' Implement a simple version of nearest-centroid classification.
#' Centroids can be computed from labels, or given directly.
#'
#' @param df - a genes by samples matrix
#' @param samples - a single [new] sample to be labeled; if NULL, all columns of df are (re-)labeled; default is NULL
#' @param distance - vector distance, could be 'cor' (for 1-pearson correlation) or 'eucl' (for Euclidean)
#' @param labels - sample labels, could be NULL
#' @param centroids - pre-computed centroids, could be NULL
#' @return corresponding labels for all samples
#' @export
#' @examples
#' ...

nearestCentroid = function(df, samples=NULL, distance = 'spearman', labels=NULL, centroids = NULL){


    if (is.null(labels) & is.null(centroids)){
        print('clusteringTools.R::nearestCentroid(): need either class labels or pre-computed centroids, both cannot be null.')
        return(NA)
    }

   stopifnot(distance %in% c('spearman','pearson','euclidean'))   # distance can be 1-correlation or Euclidean

   if (!is.null(labels)){

       stopifnot( length(labels) == ncol(df) )  # compute the centroids using the assigned labels

       # compute centroids
       centroids = matrix(NA, nrow = nrow(df), ncol = length(unique(labels)) )
       cnt = 0
       for (label in unique(labels)){
           cnt = cnt + 1
           centroids[,cnt] = apply( df[, label==labels, drop=FALSE], 1, mean )
       }
       colnames(centroids) = unique(labels)
       rownames(centroids) = rownames(df)

   } else if (!is.null(centroids)){

       stopifnot( nrow(centroids) <= nrow(df) )  # cannot have more genes in the centroids
       # nooo should check that the rownames of centroids are a subset of the rownames of df
       if ( sum(is.na(match( rownames(centroids),rownames(df) )))>0 ){
           print('clusteringTools.R::nearestCentroid(): there are genes in the centroids that missing in the data.')
           return(NA)
       }
   }

   if (is.null(samples)){
       # compute for all samples in df
       samples = df
   }

    relabeled = c()
    for (sample in 1:ncol(samples)){
        sampleDist = apply( centroids, 2, distanceFunction, samples[,sample])
        relabeled = append( relabeled, colnames(centroids)[ which.min(sampleDist) ] )
    }
    names(relabeled) = colnames(samples)

    return(relabeled)

}

