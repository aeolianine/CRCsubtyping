#' This functions shows how to download data from Synapse.
#' The example is tailored to table-format data: an expression matrix
#' Note: requires synapse login, can be automated by creating a .synapseConfig file
#'       with the following contents:
#'       [authentication]
#'       username: username
#'       password: password
#'
#' @param synapseID - the synapse identifier associated with that dataset
#' @param sep - the table separator, ex: '\t', ','; default is tab
#' @seealso \code{\link{get_GSE33113_synapse}}
#' @return data frame of the requested dataset
#' @export
#' @examples
#' tab = getSynapseTable('syn2363559')
#' print(tab[1:5,1:3])

getSynapseTable = function(synapseID, sep = '\t'){

	if (!require(synapseClient)){

		source('http://depot.sagebase.org/CRAN.R')
		# this is a wrapper for install.packages(, repos=), where repos is a Sage repository
		pkgInstall("synapseClient")
		library(synapseClient)
	}

	synapseLogin()

	sampleTable = read.table(synGet(synapseID)@filePath, sep=sep,
	                         header = TRUE, row.names = 1, check.names=FALSE)

}


#' Downloads the frozen-RMA-normalized of the GSE33113 expression from Synapse.
#' Note 1: synapse ID "syn2363559"
#' Note 2: requires synapse login, can be automated by creating a .synapseConfig file.
#' @param NULL
#' @seealso \code{\link{getSynapseTable}}
#' @return data frame of the expression matrix for the GSE33113 dataset
#' @export
#' @examples
#' exprs = get_GSE33113_synapse()
#' print(dim(exprs))
#' print(exprs[1:3,1:3])

get_GSE33113_synapse = function(){

	tab = getSynapseTable("syn2363559")
	# the first 6 samples are normal, not tumor
	tab = tab[, 7:ncol(tab)]

	# here summarize probes to genes
  print('get_GSE33113_synapse(): probe to gene summarization is not implemented yet.')
  return(tab)

}


get_GSE33113_GEO = function(){
  #Example on how to run fRMA on raw data.
  #https://github.com/Sage-Bionetworks/crcsc/blob/master/groups/F/normalization/frma_datasets.r

  print('get_GSE33113_GEO(): This function is not implemented yet.')
}


#' Download the merged (Illumina -HiSeq and -GA ) TCGA expression data from Synapse.
#' Note 1: synapse ID "syn2325328"
#' Note 2: requires synapse login, can be automated by creating a .synapseConfig file.
#'
#' @param NULL
#' @return data frame of the expression data, genes are rownames, samples as colnames
#' @seealso \code{\link{getSynapseTable}}
#' @export
#' @examples
#' tcga_expr = get_TCGA_synapse()
#' print(dim(tcga_exprs))
#' print(tcga_exprs[1:3,1:4])

get_TCGA_synapse = function(){
  tab = getSynapseTable("syn2325328")
}

