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
#' dat1 = cbind(c(1,2,3,3), c(-5,7,3,8), c(10,-4,5,0))
#' dat2 = cbind(c(10,2.6,-30,39), c(-15,7.4,-3.2,10), c(5,4,-5,0.5))
#' rownames(dat1) = c('g1', 'g2', 'g3', 'g4')
#' rownames(dat2) = c('g1', 'g2', 'g3', 'g4')
#' print(combat_to_ref(dat1, dat2))

combat_to_ref = function(myDataset, refDataset){

  # make sure there are genes to intersect
  stopifnot(!is.null(rownames(myDataset)))
  stopifnot(!is.null(rownames(refDataset)))

 	# source the M-combat script from github
  require(RCurl)
  eval( parse(text = getURL('https://raw.githubusercontent.com/aeolianine/M-ComBat/master/MComBatRScript.R')) )
  # this loads the M.COMBAT function

  # combine two datasets
  myDataset = myDataset[rowSums(abs(myDataset))!=0, ]  # remove zero rows (ComBat will not converge)
  refDataset = refDataset[rowSums(abs(refDataset))!=0, ]  # remove zero rows (ComBat will not converge)

  commonGenes = intersect(rownames(myDataset), rownames(refDataset))
  if (length(commonGenes)==0){ print('The datasets do not share any genes. Cannot normalize.'); return(NULL) }
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
