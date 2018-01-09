#' Normalize a given dataset to a reference dataset, using M-ComBat.
#' Source: https://github.com/aeolianine/M-ComBat/blob/master/MComBatRScript.R
#' Reference: Stein et al, Removing batch effects from purified plasma cell gene expression microarrays with modified ComBat,
#'            BMC Bioinformatics2015 16:63, https://doi.org/10.1186/s12859-015-0478-3
#'
#' @param myDataset - the given dataset to be normalized
#' @param refDataset - the reference dataset to "be normalized to"
#' @return a matrix corresponding to the normalized "myDataset"
#'         together with "refDataset" (pasted together column-wise)
#' @export
#' @examples
#' print('combatToRef(): The examples are not written yet!')

combatToRef = function(myDataset, refDataset, plot=FALSE){

  # source the M-combat script from github
  require(RCurl)
  eval( parse(text = getURL('https://raw.githubusercontent.com/aeolianine/M-ComBat/master/MComBatRScript.R')) )
  # this loads the M.COMBAT function

  # combine two datasets
  myDataset = myDataset[rowSums(myDataset)!=0, ]  # remove zero rows (ComBat will not converge)
  refDataset = refDataset[rowSums(refDataset)!=0, ]  # remove zero rows (ComBat will not converge)

  commonGenes = intersect(rownames(myDataset), rownames(refDataset))
  if (isTRUE(all.equal(myDataset[commonGenes,], refDataset[commonGenes,]))){
    print('combatToRef(): own and reference datasets are the same dataset. Nothing to do.')
    return(myDataset)
  }

  mat = cbind(myDataset[commonGenes,], refDataset[commonGenes,])
  batch = c(rep('myDataset', ncol(myDataset)), rep('refDataset', ncol(refDataset)))

	# remove zero rows from mat
	mat = mat[rowSums(mat)!=0,]

	matNorm = M.COMBAT(as.matrix(mat), batch = batch, center = 'refDataset', mod = c( rep(1,ncol(mat)) ) )

	# clean up
  rm(myDataset, refDataset, commonGenes, mat, batch)
	return( matNorm )

}

#' M-combat with the GSE35896 dataset as reference.
#' @param mat - the matrix to be normalized to GSE35896
#' @return a list of two variables: the normalized dataset alone, and combined with GSE35896
#' @export
#' @examples
#' org = loadPilotOrganoidExpression(meanCenter=FALSE, Log=TRUE)
#' orgT = org[, grepl('t', colnames(org))]
#' out = combatToGSE35896(orgT)
#' print(names(out))
#' print(head(out$alone))

OLD_combatToGSE35896 = function(mat, plot = FALSE){

	# remove zero rows from mat
	zeros = rowSums(mat)==0
	mat = mat[!zeros,]

	y = loadGSE35896()$gex

	if (isTRUE(all.equal(mat, y))){
		print('combatToGSE35896(): no normalization performed, this is the reference dataset.')
		return( list('alone'=mat, 'together'=mat) )
	}

	# remove zeros from GSE35896
	zeros = rowSums(y)==0
	y = y[!zeros,]

	commonGenes = intersect(rownames(mat), rownames(y))
	if (isTRUE(all.equal(mat[commonGenes,],y[commonGenes,]))){
		return(list('alone' = mat))
	}
	matGSE = cbind(mat[commonGenes,], y[commonGenes,])

	zeros = rowSums(matGSE)==0
	matGSE = matGSE[!zeros,]

	if (plot){
		source('~//tools/generalPlottingTools.R')
		par(mfrow=c(1,2))
		MDS(matGSE, types = c(rep('data', ncol(mat)), rep('ref', ncol(y))), levelColors = c('red', 'blue'), title = 'before')
		lim = par('usr')
	}


	library(sva, quietly=TRUE)
	source('~//tools/M-ComBat/MComBatRScript.R') # get m-combat

	temp = M.COMBAT(as.matrix(matGSE), batch = c(rep('data', ncol(mat)), rep('ref', ncol(y))), center = 'ref', mod = as.matrix( rep(1,ncol(matGSE)) ) )

	if (plot){
		MDS(temp[, c((ncol(mat)+1):ncol(temp), 1:ncol(mat))], types = c(rep('ref', ncol(y)), rep('data', ncol(mat))),
			levelColors = c('red', 'blue'), xlim = c(lim[1], lim[2]), ylim = c(lim[3],lim[4]), title = 'after adjustment to GSE35896')
		par(mfrow=c(1,1))
	}


	return( list('alone'=temp[,1:ncol(mat)], 'together'=temp) )
}


#' M-combat with the GSE13294 and the GSE14333 datasets as reference.
#' @param mat - the matrix to be normalized to GSE13294/GSE14333
#' @return a list of two variables: the normalized dataset alone, and combined with GSE13294/GSE14333
#' @export
#' @examples
#' org = loadPilotOrganoidExpression(meanCenter=FALSE, Log=TRUE)
#' orgT = org[, grepl('t', colnames(org))]
#' out = combatToGSE13294_GSE14333(orgT)
#' print(names(out))
#' print(head(out$alone))

OLD_combatToGSE13294_GSE14333 = function(mat, plot = FALSE){

	# remove zero rows from mat
	zeros = rowSums(mat)==0
	mat = mat[!zeros,]

	y1 = loadGSE13294()$gex
	y2 = loadGSE14333()$gex

	if (isTRUE(all.equal(y1, mat)) | isTRUE(all.equal(y2, mat))){
		print('combatToGSE13294_GSE14333(): no normalization performed, this is the reference dataset.')
		return( list('alone'=mat, 'together'=mat) )
	}

	# combat the two datasets together
	library(sva, quietly=TRUE)

    genesToRemove1 = which(rowSums(y1)==0)
    genesToRemove2 = which(rowSums(y2)==0)

    if (length(genesToRemove1)>0){ y1 = y1[-genesToRemove1, ] }
    if (length(genesToRemove2)>0){ y2 = y2[-genesToRemove2, ] }

    print(paste('combatToGSE13294_GSE14333(): removed',length(genesToRemove1),'zero-count genes from GSE13294.',sep=' '))
	print(paste('combatToGSE13294_GSE14333(): removed',length(genesToRemove1),'zero-count genes from GSE14333.',sep=' '))
    commonGenes = intersect(rownames(y1), rownames(y2))
    y = as.matrix( cbind(y1[commonGenes,], y2[commonGenes,]) )

    # combat
    y = ComBat(y, batch = c(rep('1',ncol(y1)), rep('2',ncol(y2))), mod = c(rep(1,ncol(y))) )

	# remove zeros from GSE13294 & GSE14333, there shouldn't be any
	zeros = rowSums(y)==0
	y = y[!zeros,]

	commonGenes = intersect(rownames(mat), rownames(y))
	matGSE = cbind(mat[commonGenes,], y[commonGenes,])

	zeros = rowSums(matGSE)==0
	matGSE = matGSE[!zeros,]

	if (plot){
		source('~//tools/generalPlottingTools.R')
		par(mfrow=c(1,2))
		MDS(matGSE, types = c(rep('data', ncol(mat)), rep('ref', ncol(y))), levelColors = c('red', 'blue'), title = 'before')
		lim = par('usr')
	}


	library(sva, quietly=TRUE)
	source('~//tools/M-ComBat/MComBatRScript.R') # get m-combat

	temp = M.COMBAT(as.matrix(matGSE), batch = c(rep('data', ncol(mat)), rep('ref', ncol(y))), center = 'ref', mod = as.matrix( rep(1,ncol(matGSE)) ) )

	if (plot){
		MDS(temp[, c((ncol(mat)+1):ncol(temp), 1:ncol(mat))], types = c(rep('ref', ncol(y)), rep('data', ncol(mat))),
			levelColors = c('red', 'blue'), xlim = c(lim[1], lim[2]), ylim = c(lim[3],lim[4]), title = 'after adjustment to GSE13294_GSE14333')
		par(mfrow=c(1,1))
	}

	return( list('alone'=temp[,1:ncol(mat)], 'together'=temp) )
}

