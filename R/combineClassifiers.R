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

OLD_combineClassifiers = function(matR, plot = FALSE, normalize = TRUE){

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
