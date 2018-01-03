Example on how to run fRMA on raw data.
https://github.com/Sage-Bionetworks/crcsc/blob/master/groups/F/normalization/frma_datasets.r


getSynapseTable = function(synapseID){

	if (!require(synapseClient)){
		
		source('http://depot.sagebase.org/CRAN.R')
		# this is a wrapper for install.packages(, repos=), where repos is a Sage repository
		pkgInstall("synapseClient")  

		library(synapseClient)
	}

	synapseLogin()

	sampleTable = read.table(synGet(synapseID)@filePath, sep="\t",header = TRUE,row.names = 1,check.names=FALSE)

}

get_GSE33113_synapse = function(){

	tab = getSynapseTable("syn2363559")

}


get_GSE33113_GEO = function(){

}