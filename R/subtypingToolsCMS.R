#' Applying the random forest classifier of the CMS consortium (CMSclassifier::classifyCMS.RF)
#' @param mat - the matrix to be subtyped, genes by samples, genes are typically represented as gene symbols
#' @return - table of subtype information, columns are predicted CMS,
#'          nearest CMS and posterior probabilities of CMS assignment;
#'          rownames are samples
#' @export
#' @examples
#' mat = get_TCGA_synapse()
#' res = subtypeCMS.RF(mat)

subtypeCMS.RF = function(mat){

	# first make sure that there is more than one sample
	stopifnot( !is.null(dim(mat)) )
	stopifnot( ncol(mat)>1 )

	# convert gene names to entrez id nomenclature, any(is.na(as.numeric(gene names))) should be FALSE
	if ( !any(is.na(suppressWarnings(as.numeric( rownames(mat) ))))) {
	  # entrez identifiers
	  # do nothing
	  geneSymbol = FALSE
	} else if ( all(grepl('ENSG', rownames(mat))) ){
	  # ensembl ids
	  print('subtypeCMS.RF(): Ensembl-to-entrez id conversion is not implemented yet.')
	  geneSymbol = FALSE
	  return(NA)
	} else if ('KRAS' %in% rownames(mat) | 'FAP' %in% rownames(mat)){
    # these should be gene symbols, and there are probably better ways to test for this
	  geneSymbol = TRUE
	  matO = mat # keep the original
	  rownames(mat) = symbol2entrez(rownames(mat))[rownames(mat)]
	}

	if (!require(CMSclassifier)){
		# install CMSclassifier from github
		library(devtools)
		install_github("Sage-Bionetworks/CMSclassifier")
	}

	res = CMSclassifier::classifyCMS.RF(as.data.frame(mat), center = TRUE, minPosterior = 0.5)

	return(res)
}


#' Applying the single-sample classifier of the CMS consortium (CMSclassifier::classifyCMS.SSP)
#' @param mat - the matrix to be subtyped, genes by samples, genes are typically in gene symbol format
#' @return - table of samples and their predicted, nearest CMS assignment, as well as min/median and max
#'           correlations to predefined CMS1, CMS2, CMS3 and CMS4 centroids
#' @export
#' @examples
#' mat = get_TCGA_synapse()
#' res = subtypeCMS.SSP(mat)

subtypeCMS.SSP = function(mat){

  stopifnot( !is.null(dim(mat)) )
  stopifnot( ncol(mat)>1 )

  # convert gene names to entrez id nomenclature, any(is.na(as.numeric(gene names))) should be FALSE
  if ( !any(is.na(suppressWarnings(as.numeric( rownames(mat) ))))) {
    # entrez identifiers
    # do nothing
    geneSymbol = FALSE
  } else if ( all(grepl('ENSG', rownames(mat))) ){
    # ensembl ids
    print('subtypeCMS.RF(): Ensembl-to-entrez id conversion is not implemented yet.')
    geneSymbol = FALSE
    return(NA)
  } else if ('KRAS' %in% rownames(mat) | 'FAP' %in% rownames(mat)){
    # these should be gene symbols, and there are probably better ways to test for this
    geneSymbol = TRUE
    matO = mat # keep the original
    rownames(mat) = symbol2entrez(rownames(mat))[rownames(mat)]
  }

  if (!require(CMSclassifier)){
    # install CMSclassifier from github
    library(devtools)
    install_github("Sage-Bionetworks/CMSclassifier")
  }

  res = CMSclassifier::classifyCMS.SSP(mat, minCor=.15, minDelta=.06)

	return(res)

}

#' Print quick comparison stats on the consensus between the RF and SSP classifiers
#' This function uses the results of subtypeCMS.RF and subtypeCMS.SSP
#'
#' @param resRF - the results data frame output of subtypeCMS.RF
#' @param resSSP - the results data frame output of subtypeCMS.SSP
#' @return table (data.frame) of samples that disagree in assignment
#'         with their subtype and probabilities of assignment
#' @export
#' @examples
#' mat = get_TCGA_synapse()
#' resRF = subtypeCMS.RF(mat)
#' resSSP = subtypeCMS.SSP(mat)
#' compare_RF_to_SSP(resRF, resSSP)

compare_RF_to_SSP = function(resRF, resSSP){

  # compare resRF$RF.predictedCMS and resSSP$SSP.predictedCMS
  x1 = sum( resRF$RF.predictedCMS == resSSP$SSP.predictedCMS, na.rm=TRUE )
  print(paste0('compare_RF_to_SSP(): In ', x1, ' cases the RF and SSP classifiers *agree* and the assignment is significant (',
               format(100*x1/nrow(resRF), digits = 3), ' percent).'))

  x2 = sum( resRF$RF.nearestCMS == resSSP$SSP.nearestCMS & (is.na(resRF$RF.predictedCMS) | is.na(resSSP$SSP.predictedCMS)) )
  print(paste0('compare_RF_to_SSP(): In ', x2, ' cases the RF and SSP classifiers *agree* and the assignment is not significant (',
               format(100*x2/nrow(resRF), digits = 3), ' percent).'))

  print(paste0('compare_RF_to_SSP(): The total percentage of sample with the same assigned subtype is: ', format(100*(x1+x2)/nrow(resRF), digits = 3), '.' ))

  x3 = resRF$RF.predictedCMS != resSSP$SSP.predictedCMS
  print(paste0('compare_RF_to_SSP(): In ', sum(x3, na.rm = TRUE), ' cases the RF and SSP classifiers *disagree* and the assignment is significant (',
               format(100*sum(x3, na.rm = TRUE)/nrow(resRF), digits = 3), ' percent).'))

  x4 = resRF$RF.nearestCMS != resSSP$SSP.nearestCMS & (is.na(resRF$RF.predictedCMS) | is.na(resSSP$SSP.predictedCMS))
  print(paste0('compare_RF_to_SSP(): In ', sum(x4), ' cases the RF and SSP classifiers *disagree* and the assignment is not significant (',
               format(100*sum(x4)/nrow(resRF), digits = 3), ' percent).'))


  # samples with disagreement and significant assignment
  x3 = rownames(resRF)[x3]
  x3 = x3[!is.na(x3)]
  mismatchSig = paste(resRF[x3,'RF.predictedCMS'], resSSP[x3, 'SSP.predictedCMS'], sep = '-')
  print('compare_RF_to_SSP(): Typical mismatch between RF and SSP samples, significant assignments:')
  print(table(mismatchSig))

  x4 = rownames(resRF)[x4]
  mismatchSig = paste(resRF[x4,'RF.nearestCMS'], resSSP[x4, 'SSP.nearestCMS'], sep = '-')
  print('compare_RF_to_SSP(): Typical mismatch between RF and SSP samples, nonsignificant assignments:')
  print(table(mismatchSig))

  tab = cbind(resRF[c(x3,x4),], resSSP[c(x3,x4), ])
  # clean up
  rm(x1, x2, x3, x4, mismatchSig, resRF, resSSP)

	return(tab)
}

#' Compare the subtypeCMS.RF assigned nearest subtypes to the RF model (gold reference) training labels.
#' This is only possible for datasets that are used within the reference set.
#' Note: Currently only TCGA is implemented.
#' @param resRF - output of subtypeCMS.RF
#' @param whichDataset - currently only "TCGA" is a valid input
#' @return table (data.frame) of subtypes, training labels and posterior probabilities for samples that do not have matching labels;
#'         also a message with the percentage of matching sample labels
#' @export
#' @examples
#' mat = get_TCGA_synapse()
#' resRF = subtypeCMS.RF(mat)
#' compare_to_training_labels(resRF, whichDataset = 'TCGA')

compare_to_training_labels = function(resRF, whichDataset = NULL){

  stopifnot(!is.null(whichDataset))
  stopifnot(whichDataset %in% c('TCGA'))

  if (whichDataset == 'TCGA'){
    stopifnot( all(grepl('TCGA', rownames(resRF))) )
    extractLabels = finalModel$y[grepl('TCGA', names(finalModel$y))]
    names(extractLabels) = sapply(names(extractLabels), function(x) substr(x, nchar("tcga_rnaseqAll.")+1, nchar(x)), USE.NAMES=FALSE)

    commonSamples = intersect(names(extractLabels), rownames(resRF))
    # compare extractLabels[commonSamples] and resRF[commonSamples, 'RF.predictedCMS']
    x = extractLabels[commonSamples] == resRF[commonSamples, 'RF.nearestCMS']
    print(paste0('compare_to_training_labels(): ', format(100*sum(x)/length(x), digits = 3), ' labels match.' ))
    print('compare_to_training_labels(): Stats for samples that do not match ->')
    tab = cbind( extractLabels[commonSamples[!x]], resRF[commonSamples[!x], ] )
    tab = tab[, c(1,6:7,2:5)]
    colnames(tab)[1] = 'training label'
    print(tab)
  }

}
