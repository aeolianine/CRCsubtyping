#' This function is used within convertRownamesToEZID.
#' 
#' @param geneName - this is a character, the gene name to be translated.
#' @param dict - this is the map; the values are the desired translation, the names(dict) are the geneName(s)
#' @return - the value of dict[ind], such that names(dict)[ind] = geneName
#' @export
#' @examples
#' dict = c('a','b','c')
#' names(dict) = c('1','2','3')
#' geneSymbol2EZID('2',dict=dict)   # should return 'b'


geneSymbol2EZID = function(geneName, dict = map){

  	if (geneName %in% names(dict)){
  		ind = which(names(dict)==geneName)
    	# the [1] takes care of multiple id assignments (there shouldn't be any ...)
		ezid = unlist( as.list( dict[ind] ) )[1]   
		names(ezid) = NULL
		return(ezid)
  	} else {
  		# if this geneName is not mapped, simply return it
  		#print(paste(c('geneSymbol2EZID(): ',geneName, ' is not mapped.'), collapse=''))
		return( NA )
	}	

}

#' Convert Entrez ID to gene symbols
#' 
#' @param character vector of entrez ids
#' @return character vector of gene symbols
#' @export
#' @examples
#' ....

entrez2symbol = function(geneids){

  library(org.Hs.eg.db, quietly=TRUE)
  library(annotate, quietly=TRUE)
  return( getSYMBOL(geneids, data='org.Hs.eg') )

}


#' Convert gene symbol to Entrez ID
#' 
#' @param character vector of gene symbols
#' @return character vector of entrez ids
#' @export
#' @examples
#' ....

convertRownamesToEZID = function(geneNames){

  library(org.Hs.eg.db, quietly=TRUE)

  # geneNames is the rownames of the data matrix

  # load map
  map = org.Hs.egSYMBOL2EG
  mappedKeys = keys(map)
  mappedGenes = intersect(mappedKeys, geneNames)

  mapDict = as.list( links(map[mappedGenes])[,'gene_id'] )
  names(mapDict) = links(map[mappedGenes])[,'symbol']

  genesEZID = sapply(geneNames, geneSymbol2EZID, mapDict )
  genesEZID = unlist(genesEZID)
  names(genesEZID) = NULL
  
  return(genesEZID)
}


# param L, list of names to be cleaned from duplicates
# return list without duplicates
# removeDuplicates = function(L){
#	L[!(duplicated(L) | duplicated(L, fromLast = TRUE))]
#}


#' Applying the random forest classifier of the CMS consortium.
#' @param mat - the matrix to be subtyped, genes by samples, genes are typically in gene symbol format
#' @param plot - Boolean, TRUE or FALSE
#' @return - table of subtype information, columns are predicted CMS, nearest CMS and posterior probabilities of CMS assignment
#' @export

subtypeCMS.RF = function(mat, plot=FALSE){

	if (!require(CMSclassifier)){
		library(devtools, quietly=TRUE)
		install_github("Sage-Bionetworks/CMSclassifier")
		library(CMSclassifier)
	}

	# first make sure that there is more than one sample
	stopifnot( !is.null(dim(mat)) )
	stopifnot( ncol(mat)>1 )

	# make sure the matrix is not mean-centered yet ...
	stopifnot(mean(rowMeans(mat))>0.5)

	# change gene names from symbol to entrez id nomenclature
	entrezids = convertRownamesToEZID(rownames(mat))
	names(entrezids) = rownames(mat)
	entrezids = entrezids[!is.na(entrezids)]
	mat = mat[names(entrezids),]
	rownames(mat) = entrezids
	
	res = classifyCMS(as.data.frame(mat), method='RF')

	if (plot){
		temp = t(res$RF.details[,1:4])
		temp = temp[, order(temp['RF.CMS2.posteriorProb',],-temp['RF.CMS4.posteriorProb',])]
		barplot(temp, col = c('orange3','blue2','lightcoral','lightgreen'), 
										 legend=c('CMS1','CMS2','CMS3','CMS4'), space = 0, las = 2)
	
		if (FALSE){
		fit = finalModel
		imp = fit$importance[, c('CMS1','CMS2','CMS3','CMS4')]
		genes = rownames(imp)
		stopifnot( sort(genes) == sort(listModelGenes()) )

		geneanno=data.frame(CMS1 = imp[, 'CMS1'] > 0.01, 
							CMS2 = imp[, 'CMS2'] > 0.01,
							CMS3 = imp[, 'CMS3'] > 0.01,
							CMS4 = imp[, 'CMS4'] > 0.01, stringsAsFactors=FALSE)
		
		genes = intersect(genes, rownames(mat))
		temp = mat[genes,]
		rownames(temp) = entrez2symbol(genes)
		geneanno = geneanno[genes,]
		rownames(geneanno) = entrez2symbol(rownames(geneanno))
			
		source('~//tools/generalPlottingTools.R')
		multiFactorHeatmap(temp-apply(temp,1,mean), data.frame(nearestSubtype = res$nearestCMS), geneanno=geneanno)
		}

	}

	return(res)

}


#' Applying the single-sample classifier of the CMS consortium.
#' @param mat - the matrix to be subtyped, genes by samples, genes are typically in gene symbol format
#' @param plot - Boolean, TRUE or FALSE
#' @return - table of subtype information, columns are predicted CMS, nearest CMS and poster
#' @export

subtypeCMS.SSP = function(mat, plot=FALSE, plotWhichSamples = c()){

	library(CMSclassifier, quietly=TRUE)

	# make sure the matrix is not mean-centered yet ...
	stopifnot(mean(rowMeans(mat))>0.5)

	# change gene names from symbol to entrez id nomenclature
	entrezids = convertRownamesToEZID(rownames(mat))
	names(entrezids) = rownames(mat)
	entrezids = entrezids[!is.na(entrezids)]
	mat = mat[names(entrezids),]
	rownames(mat) = entrezids

	res = classifyCMS(mat, method='SSP')

	if (plot){

		source('~//tools/parsingTools.R')
		alldata=c()
		
		temp = res$SSP.details[, c('SSP.min.corToCMS1', 'SSP.median.corToCMS1', 'SSP.max.corToCMS1')]
		rownames(temp) = addSuffix(rownames(temp), '-CMS1')
		colnames(temp) = c('SSP.min.corToCMS', 'SSP.median.corToCMS', 'SSP.max.corToCMS')
		alldata = rbind(alldata, temp)

		temp = res$SSP.details[, c('SSP.min.corToCMS2', 'SSP.median.corToCMS2', 'SSP.max.corToCMS2')]
		rownames(temp) = addSuffix(rownames(temp), '-CMS2')
		colnames(temp) = c('SSP.min.corToCMS', 'SSP.median.corToCMS', 'SSP.max.corToCMS')
		alldata = rbind(alldata, temp)

		temp = res$SSP.details[, c('SSP.min.corToCMS3', 'SSP.median.corToCMS3', 'SSP.max.corToCMS3')]
		rownames(temp) = addSuffix(rownames(temp), '-CMS3')
		colnames(temp) = c('SSP.min.corToCMS', 'SSP.median.corToCMS', 'SSP.max.corToCMS')
		alldata = rbind(alldata, temp)

		temp = res$SSP.details[, c('SSP.min.corToCMS4', 'SSP.median.corToCMS4', 'SSP.max.corToCMS4')]
		rownames(temp) = addSuffix(rownames(temp), '-CMS4')
		colnames(temp) = c('SSP.min.corToCMS', 'SSP.median.corToCMS', 'SSP.max.corToCMS')
		alldata = rbind(alldata, temp)				

		if (length(plotWhichSamples)>0){
			keep = c()
			for (thisSample in plotWhichSamples){
				keep = append(keep, rownames(alldata)[grepl(thisSample, rownames(alldata))])
			}
			alldata = alldata[keep,]
		}
		alldata = alldata[sort(rownames(alldata)),]

		par(mfrow=c(1,1), oma = c(1,2,1,1))
		plot(alldata[,1], 1:nrow(alldata), col = 'blue', xlab = '{min, median, max}\n correlations to subtypes', xlim=c(min(alldata),max(alldata)), yaxt='n', ylab='')
		points(alldata[,3],1:nrow(alldata), col = 'red')
		for (row in 1:nrow(alldata)){
			lines(c(alldata[row,1],alldata[row,3]), c(row,row), lwd = 0.1, col = 'black')
		}
		points(alldata[,2],1:nrow(alldata), col = 'orange')
		axis(side=2, labels = rownames(alldata), at = 1:nrow(alldata), cex.axis=0.5, las=2)



		if (FALSE){
		# CMS2 ................
		temp = res$SSP.details[, c('SSP.min.corToCMS2', 'SSP.median.corToCMS2', 'SSP.max.corToCMS2')]
		plot(temp[,1], 1:nrow(temp), col = 'blue', xlab = '{min, median, max}\n correlations to CMS2', xlim=c(min(temp),max(temp)), yaxt='n', ylab='')
		points(temp[,3],1:nrow(temp), col = 'red')
		for (row in 1:nrow(temp)){
			lines(c(temp[row,'SSP.min.corToCMS2'],temp[row,'SSP.max.corToCMS2']), c(row,row), lwd = 0.1, col = 'black')
		}
		points(temp[,2],1:nrow(temp), col = 'orange')
		axis(side=2, labels = rownames(temp), at = 1:nrow(temp), cex.lab=1, las=2)

		# CMS3 ................
		temp = res$SSP.details[, c('SSP.min.corToCMS3', 'SSP.median.corToCMS3', 'SSP.max.corToCMS3')]
		plot(temp[,1], 1:nrow(temp), col = 'blue', xlab = '{min, median, max}\n correlations to CMS3', xlim=c(min(temp),max(temp)), yaxt='n', ylab='')
		points(temp[,3],1:nrow(temp), col = 'red')
		for (row in 1:nrow(temp)){
			lines(c(temp[row,'SSP.min.corToCMS3'],temp[row,'SSP.max.corToCMS3']), c(row,row), lwd = 0.1, col = 'black')
		}
		points(temp[,2],1:nrow(temp), col = 'orange')
		axis(side=2, labels = rownames(temp), at = 1:nrow(temp), cex.lab=1, las=2)

		# CMS4 ................
		temp = res$SSP.details[, c('SSP.min.corToCMS4', 'SSP.median.corToCMS4', 'SSP.max.corToCMS4')]
		plot(temp[,1], 1:nrow(temp), col = 'blue', xlab = '{min, median, max}\n correlations to CMS4', xlim=c(min(temp),max(temp)), yaxt='n', ylab='')
		points(temp[,3],1:nrow(temp), col = 'red')
		for (row in 1:nrow(temp)){
			lines(c(temp[row,'SSP.min.corToCMS4'],temp[row,'SSP.max.corToCMS4']), c(row,row), lwd = 0.1, col = 'black')
		}
		points(temp[,2],1:nrow(temp), col = 'orange')
		axis(side=2, labels = rownames(temp), at = 1:nrow(temp), cex.lab=1, las=2)


		par(mfrow=c(4,1))
		}
	}

	return(res)

}


compare_RF_to_SSP = function(mat){


	res.SSP = subtypeCMS.SSP(mat)
	res.RF = subtypeCMS.RF(mat)


	predicted = cbind(res.SSP$predictedCMS[sort(rownames(res.SSP$predictedCMS)), ], res.RF$predictedCMS[sort(rownames(res.RF$predictedCMS)), ])
	predicted = cbind(predicted, predicted[,1]==predicted[,2])
	colnames(predicted) = c('SSP', 'RF', 'match')
	rownames(predicted) = sort(rownames(res.SSP$predictedCMS))

	nearest = cbind(res.SSP$nearestCMS[sort(rownames(res.SSP$nearestCMS)), ], res.RF$nearestCMS[sort(rownames(res.RF$nearestCMS)), ])
	nearest = cbind(nearest, nearest[,1]==nearest[,2])
	colnames(nearest) = c('SSP', 'RF', 'match')
	rownames(nearest) = sort(rownames(res.SSP$nearestCMS))

	return(list(pred = predicted, near=nearest))
}


#' Subtype samples using the temporary Sage-CMS classifier. 
#'
#' @param mat - the input gene expression matrix, genes by samples (works better if gene-wise mean centered)
#' @param impute - Boolean, impute or not the missing genes
#' @return res - factor of CMSx values {CMS1, CMS2, CMS3, CMS4}. The factor/vector names are the sample names (or colnames(mat))
#' @export
#' @examples
#' ....


subtypeCMSOld = function(mat, impute=FALSE){

	library(CMSclassifier, quietly=TRUE)
	library(org.Hs.eg.db, quietly=TRUE)
	library(annotate, quietly=TRUE)

	CMSgenes = getSYMBOL(listModelGenes(), data='org.Hs.eg')
	names(CMSgenes) = NULL  # remove the ENTREZID gene names

	missingGenes = CMSgenes[ is.na( match(CMSgenes, rownames(mat)) ) ]

	if (length(missingGenes)>0){ 
		print('subtypeCMS: There are some classifier genes missing from the data:')
		print(missingGenes)

		if (impute){

			print('Imputing the missing genes ... this may take awhile.')
			for (gene in missingGenes){
				mat = rbind(mat, rep(NA, ncol(mat)))
  				rownames(mat)[nrow(mat)] = gene
			}
			library(organoidsProject, quietly=TRUE)
			y = loadTCGARNAseqData(meanCenter=TRUE, whichPlatform='IlluminaHiSeq')$tumor
			commonGenes = intersect(rownames(y), rownames(mat))

			#src = c( rep('data', ncol(mat)), rep('TCGA-RNAseq', ncol(y)) )
			temp = cbind(mat[commonGenes,],y[commonGenes,])
			#library(sva)
			#mat = ComBat( as.matrix(mat), batch = src, mod = c(rep(1,ncol(mat))) )
			library(SpatioTemporal, quietly=TRUE)
			out = SVDmiss(temp)
			temp = out$Xfill
			temp = temp[,1:ncol(mat)]

			temp = temp[CMSgenes,]
			print('The imputed genes:')
			print(rowMedians(as.matrix(temp[missingGenes,])))
			rownames( temp ) = convertRownamesToEZID(rownames( temp ))
			
			res = classifyCMS(temp)
			return(res)

		}
		return(list(NA, missingGenes))
	}

	# NKI data CMS-Sage subtypes:
	temp = mat
	rownames( temp ) = convertRownamesToEZID(rownames( temp ))
	res = classifyCMS(temp)
	
	return(res)
}


#' M-combat with the TCGA COAD RNA-seq dataset as reference.
#' @param mat - the matrix to be normalized to TCGA, with entries in log space
#' @return a list of two variables: the normalized dataset alone, and combined with TCGA
#' @export
#' @examples
#' org = loadPilotOrganoidExpression(meanCenter=FALSE, Log=TRUE)
#' orgT = org[, grepl('t', colnames(org))]
#' orgT = combatToTCGA(orgT)
#' print(head(orgT))

combatToTCGA = function(mat, plot=FALSE){

	# remove zero rows from mat	
	zeros = rowSums(mat)==0
	mat = mat[!zeros,]

	y = loadTCGARNAseqData(meanCenter=FALSE, Log=TRUE, removeZeros = TRUE)$tumor

	commonGenes = intersect(rownames(mat), rownames(y))
	if (isTRUE(all.equal(mat[commonGenes,],y[commonGenes,]))){
		return(list('alone' = mat))
	}

	matTCGA = cbind(mat[commonGenes,], y[commonGenes,])

	zeros = rowSums(matTCGA)==0
	matTCGA = matTCGA[!zeros,]
	
	if (plot){
		source('~//tools/generalPlottingTools.R')
		par(mfrow=c(1,2))
		MDS(matTCGA, types = c(rep('data', ncol(mat)), rep('ref', ncol(y))), levelColors = c('red', 'blue'), title = 'before')
		lim = par('usr')
	}
		
	library(sva, quietly=TRUE)
	source('~//tools/M-ComBat/MComBatRScript.R') # get m-combat

	temp = M.COMBAT(as.matrix(matTCGA), batch = c(rep('data', ncol(mat)), rep('TCGA', ncol(y))), center = 'TCGA', mod = as.matrix( rep(1,ncol(matTCGA)) ) )
	
	if (plot){
		MDS(temp[, c((ncol(mat)+1):ncol(temp), 1:ncol(mat))], types = c(rep('ref', ncol(y)), rep('data', ncol(mat))), 
			levelColors = c('red', 'blue'), xlim = c(lim[1], lim[2]), ylim = c(lim[3],lim[4]), title = 'after adjustment to TCGA')
		par(mfrow=c(1,1))
	}

	return( list('alone'=temp[,1:ncol(mat)], 'together'=temp) )
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

combatToGSE35896 = function(mat, plot = FALSE){

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

combatToGSE13294_GSE14333 = function(mat, plot = FALSE){

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



#' Combining the three main classifiers with some variations into a consensus table. [This has to be revised given the new normalization]
#'
#' @param matR - the input gene expression matrix, genes by samples; the non-mean-centered version, but logged (for RF and SSP)
#' @param plotTGFB - boolean, whether to plot TGFB expression or not, annotated by subtype
#' @return Table - table of samples by classifiers, populated by the corresponding subtypes, including a consensus column
#' @export
#' @examples
#' orgM = loadPilotOrganoidExpression(meanCenter=TRUE)
#' orgM = orgM[, grepl('t', colnames(orgM))]
#' orgR = loadPilotOrganoidExpression(meanCenter=FALSE)
#' orgR = orgR[, grepl('t', colnames(orgR))]
#' res = combineClassifiers(orgR)

combineClassifiers = function(matR, plot = FALSE, normalize = TRUE){

	samples = colnames(matR)
	details = c()
	if (normalize){ 
		out = combatToGSE35896(matR, plot=plot) 
	} else {
		out = list()
		out$'alone' = matR
	}

	# ...............................................
	matInmf = as.matrix( out$'alone' ) 
	matInmf = matInmf - apply(matInmf, 1, mean)

	inmfOne = subtypeOneAtATimeINMF(matInmf, distance = 'spearman')
	temp = unpackINMF(inmfOne)
	colnames(temp)[1] = 'inmfOne'
	details = temp[samples,]

	inmfOne = unpackINMFclustering(inmfOne)
	inmfOne = clust2list(inmfOne)

	print('inmfOne') # ..............................

	# ...............................................
	inmfAll = subtypeSamplesINMF(matInmf, distance = 'spearman')
	temp = unpackINMF(inmfAll)
	colnames(temp)[1] = 'inmfAll'
	details = cbind(details, temp[samples,])

	inmfAll = unpackINMFclustering(inmfAll)
	inmfAll = clust2list(inmfAll)

	print('inmfAll') # ..............................

	# ...............................................
	if ('together' %in% names(out)){
		matInmf = as.matrix( out$'together' )

		remove(out)
		matInmf = matInmf - apply(matInmf, 1, mean)	

		inmfWithGSE35896 = subtypeSamplesINMF(matInmf, distance = 'spearman')
		temp = unpackINMF(inmfWithGSE35896)
		colnames(temp)[1] = 'inmfWithGSE35896'
		details = cbind(details, temp[samples,])

		inmfWithGSE35896 = unpackINMFclustering(inmfWithGSE35896)
		inmfWithGSE35896 = clust2list(inmfWithGSE35896)

		print('inmf with GSE35896') # ..............................
	} 
	
	# for all the versions of the Sadanandam classifier
	# ...............................................
	if (normalize){ out = combatToGSE13294_GSE14333(matR, plot=plot) }
	# out$'alone' has been defined above
	matSad = as.matrix( out$'alone' )
	matSad = matSad - apply(matSad, 1, mean)

	sadOne = subtypeOneAtATimeSadanandam(matSad, distance = 'spearman')

	temp = as.matrix( sadOne$silhouette[, c('cluster', 'neighbor', 'sil_width')] )
	colnames(temp)[1] = 'sadOne'
	details = cbind(details, temp[samples,])

	sadOne = clust2list(sadOne$clustering)

	print('sadOne') # ..............................

	# ...............................................
	sadNear = nearestCentroidsSadanandam(as.matrix(matSad), distance = 'spearman')
	temp = as.matrix( sadNear$silhouette[, c('cluster', 'neighbor', 'sil_width')] )
	colnames(temp)[1] = 'sadNear'
	details = cbind(details, temp[samples,])

	sadNear = clust2list(sadNear$clustering)

	print('sadNear') # ..............................

	# ...............................................
	sadAll = subtypeSamplesSadanandam(as.matrix(matSad), distance = 'spearman')
	if (!is.na(sadAll$silhouette)){
		temp = as.matrix( sadAll$silhouette[, c('cluster', 'neighbor', 'sil_width')] )
	} else {
		l = clust2list(sadAll$clustering)
		temp = matrix(NA, nrow = length(l), ncol = 3)
		colnames(temp) = c('cluster', 'neighbor', 'sil_width')
		rownames(temp) = samples
		temp[,'cluster'] = l[samples]
	}
	colnames(temp)[1] = 'sadAll'
	details = cbind(details, temp[samples,])

	sadAll = clust2list(sadAll$clustering)

	print('sadAll') # ..............................

	# ...............................................
	if ('together' %in% names(out)){
		matSad = as.matrix( out$'together' )
		matSad = matSad - apply(matSad, 1, mean)

		sadAllwGSE = subtypeSamplesSadanandam(as.matrix(matSad), distance = 'spearman')
		if (!is.na(sadAllwGSE$silhouette)){
			temp = as.matrix( sadAllwGSE$silhouette[, c('cluster', 'neighbor', 'sil_width')] )
		} else {
			l = clust2list(sadAllwGSE$clustering)
			temp = matrix(NA, nrow = length(l), ncol = 3)
			colnames(temp) = c('cluster', 'neighbor', 'sil_width')
			rownames(temp) = samples
			temp[,'cluster'] = l[samples]
		}

		colnames(temp)[1] = 'sadAllwGSE'
		details = cbind(details, temp[samples,])

		sadAllwGSE = clust2list(sadAllwGSE$clustering)

		print('sad NMF style with GSE13294 and GSE14333') # ..............................
	}

	if ('inmfWithGSE35896' %in% ls() & 'sadAllwGSE' %in% ls()){
		Table = data.frame(inmfOne = inmfOne[samples], inmfAll = inmfAll[samples], 
					   inmfWithGSE35896 = inmfWithGSE35896[samples], sadOne = sadOne[samples], 
					   sadNear = sadNear[samples], sadAll = sadAll[samples], 
					   sadAllwGSE = sadAllwGSE[samples], stringsAsFactors=FALSE)
	} else if (!('inmfWithGSE35896' %in% ls()) & 'sadAllwGSE' %in% ls()) {
		Table = data.frame(inmfOne = inmfOne[samples], inmfAll = inmfAll[samples], 
					   sadOne = sadOne[samples], 
					   sadNear = sadNear[samples], sadAll = sadAll[samples], 
					   sadAllwGSE = sadAllwGSE[samples], stringsAsFactors=FALSE)
	} else if ('inmfWithGSE35896' %in% ls() & !('sadAllwGSE' %in% ls())){
		Table = data.frame(inmfOne = inmfOne[samples], inmfAll = inmfAll[samples], 
					   inmfWithGSE35896 = inmfWithGSE35896[samples], sadOne = sadOne[samples], 
					   sadNear = sadNear[samples], sadAll = sadAll[samples], 
					   stringsAsFactors=FALSE)
	} else {  # no normalization - no combination with other datasets
		Table = data.frame(inmfOne = inmfOne[samples], inmfAll = inmfAll[samples], 
					   sadOne = sadOne[samples], 
					   sadNear = sadNear[samples], sadAll = sadAll[samples], 
					   stringsAsFactors=FALSE)
	}

	# create a consensus column
	cmsLabels = list('1.1'='CMS4','1.2'='CMS1', '1.3'='CMS4', '2.1'='CMS3', '2.2'='CMS2',
				 	 'Stem.like'='CMS4', 'Goblet.like'='CMS3', 'TA'='CMS2', 'Enterocyte'='CMS2', 'Inflammatory'='CMS1',
				 	 'CMS1'='CMS1','CMS2'='CMS2','CMS3'='CMS3','CMS4'='CMS4',
				 	 '1'='CMS1','2'='CMS2','3'='CMS3','4'='CMS4')
	
	fraction = c()
	consensus = c()

	for (row in 1:nrow(Table)){
		x = as.character(Table[row,])
		x_cms = table( unlist(cmsLabels[ x ]) )
		ind_max = which(x_cms==max(x_cms))

		if (length(ind_max)>1){
			cons = paste(names(x_cms)[ind_max], collapse=',')
			f =  sum( names(x_cms)[ind_max[1]] == unlist( cmsLabels[x]) )/ncol(Table)
			fraction = append(fraction, f)
		} else if (length(ind_max)==1){
			cons = names(x_cms)[ind_max]
			f =  sum( names(x_cms)[ind_max] == unlist( cmsLabels[x]) )/ncol(Table)
			fraction = append(fraction, f)
		}
		consensus = append(consensus, cons)
	}

	Table = cbind(Table, consensus)
	colnames(Table)[ncol(Table)] = 'consensus'
	Table = cbind(Table, fraction)
	colnames(Table)[ncol(Table)] = 'consensus fraction'


	if (normalize){ out = combatToTCGA(matR, plot=plot)  } # for the RF and SSP classifiers
	mat = out$'alone'   # has been defined above

	cms.RF = subtypeCMS.RF(as.matrix(mat), plot = plot)
	cms.SSP = subtypeCMS.SSP(as.matrix(mat), plot = plot)
	nearestRF = cms.RF$nearestCMS
	nearestSSP = cms.SSP$nearestCMS

	RFdetails = cms.RF$RF.details
	SSPdetails = cms.SSP$SSP.details

	cms.RF = cms.RF$predictedCMS
	cms.SSP = cms.SSP$predictedCMS


	print('cms')

	Table = cbind(Table, cms.RF[samples,])
	colnames(Table)[ncol(Table)] = 'RF'
	Table = cbind(Table, cms.SSP[samples,])
	colnames(Table)[ncol(Table)] = 'SSP'

	temp = Table[,c('consensus','RF')]
	temp = temp[complete.cases(temp),]
	print(paste('"Consensus" vs RF', sum( as.character(temp[,1])!=as.character(temp[,2]) ), 'out of', nrow(temp), 'callable comparisons are wrong.' ) )

	temp = Table[,c('consensus','SSP')]
	temp = temp[complete.cases(temp),]
	print(paste('"Consensus" vs SSP', sum( as.character(temp[,1])!=as.character(temp[,2]) ), 'out of', nrow(temp), 'callable comparisons are wrong.' ) )

	temp = Table[,c('RF','SSP')]
	temp = temp[complete.cases(temp),]
	print(paste('RF vs SSP', sum( as.character(temp[,1])!=as.character(temp[,2]) ), 'out of', nrow(temp), 'callable comparisons are wrong.' ) )

	temp = Table[, c('consensus', 'RF', 'SSP')]
	cnt = 0
	for (row in 1:nrow(temp)){
		x1 = as.character(temp[row, 'consensus'])
		x2 = as.character(temp[row, 'SSP'])
		x3 = as.character(temp[row, 'RF'])

		if ( isTRUE(all.equal(x1,x2)) & isTRUE(all.equal(x1,x3)) ){ cnt = cnt + 1 }
	}
	print(paste(cnt, 'out of', nrow(temp), 'samples agree three-ways {consensus, SSP, RF}.'))

	temp = Table[, c('consensus', 'RF', 'SSP')]
	temp = temp[complete.cases(temp),]
	cnt = 0
	for (row in 1:nrow(temp)){
		x1 = as.character(temp[row, 'consensus'])
		x2 = as.character(temp[row, 'SSP'])
		x3 = as.character(temp[row, 'RF'])

		if ( isTRUE(all.equal(x1,x2)) & isTRUE(all.equal(x1,x3)) ){ cnt = cnt + 1 }
	}
	print(paste(cnt, 'out of', nrow(temp), 'samples agree three-ways {consensus, SSP, RF} including complete cases only.'))

	Table[,'SSP'] = as.character(Table[,'SSP'])
	Table[,'SSP'][is.na(Table[,'SSP'])] = 'unknown'


	Table = cbind(Table, nearestRF[samples,])
	colnames(Table)[ncol(Table)] = 'RF-nearest'
	Table = cbind(Table, nearestSSP[samples,])
	colnames(Table)[ncol(Table)] = 'SSP-nearest'

	rownames(Table) = samples

	# save the details from the RF and the SSP
	details = cbind(details, RFdetails[samples, ] )
	details = cbind(details, SSPdetails[samples, ] )


	if (plot){
		source('~//tools/generalPlottingTools.R')

		par(mfrow = c(3,2), mar=c(3,3,4,0.5), oma = c(3,1,1,1))

		types = unique( as.character(Table[,'consensus']) )
		types = sort(types)
		colors = rep('gray', length(types))
		colors[types == 'CMS1'] = 'orange3'
		colors[types == 'CMS2'] = 'blue2'
		colors[types == 'CMS3'] = 'lightcoral'
		colors[types == 'CMS4'] = 'lightgreen'

		plotGeneExpression(colSums(mat[c('TGFB1','TGFB2','TGFB3'), ]), types = as.character(Table[,'consensus']), colors = colors, title = 'consensus: TGFB1+TGFB2+TGFB3')
		plotGeneExpression(unlist(mat[c('GREM1'), ]), types = as.character(Table[,'consensus']), colors = colors, title = 'consensus: GREM1')

		types = unique( as.character(Table[,'RF-nearest']) )
		types = sort(types)
		colors = rep('gray', length(types))
		colors[types == 'CMS1'] = 'orange3'
		colors[types == 'CMS2'] = 'blue2'
		colors[types == 'CMS3'] = 'lightcoral'
		colors[types == 'CMS4'] = 'lightgreen'

		plotGeneExpression(colSums(mat[c('TGFB1','TGFB2','TGFB3'), ]), types = as.character(Table[,'RF-nearest']), colors = colors, title = 'RF-nearest: TGFB1+TGFB2+TGFB3')
		plotGeneExpression(unlist(mat[c('GREM1'), ]), types = as.character(Table[,'RF-nearest']), colors = colors, title = 'RF-nearest: GREM1')
		
		types = unique( as.character(Table[,'SSP-nearest']) )
		types = sort(types)
		colors = rep('gray', length(types))
		colors[types == 'CMS1'] = 'orange3'
		colors[types == 'CMS2'] = 'blue2'
		colors[types == 'CMS3'] = 'lightcoral'
		colors[types == 'CMS4'] = 'lightgreen'
		
		plotGeneExpression(colSums(mat[c('TGFB1','TGFB2','TGFB3'), ]), types = as.character(Table[,'SSP-nearest']), colors = colors, title = 'SSP-nearest: TGFB1+TGFB2+TGFB3')
		plotGeneExpression(unlist(mat[c('GREM1'), ]), types = as.character(Table[,'SSP-nearest']), colors = colors, title = 'SSP-nearest: GREM1')
		
		#if ('GREM2' %in% rownames(mat)){
		#	plotGeneExpression(unlist(mat[c('GREM2'), ]), types = as.character(Table[,'SSP']), colors = c('orange3','blue2','lightcoral','lightgreen','gray'), title = 'SSP: GREM2')
		#}
		par(mfrow = c(1,1))
	}


	rm(list = setdiff(ls(),c('Table','details')))
	return(list( 'Table' = Table, 'details' = details) )
}