rm(list = ls())

if (FALSE){
# iNMF gene signature stored as entrez identifiers

t1 = scan('12920_2012_343_MOESM2_ESM_entrez.txt', sep = '\t', skip=1)
t2 = scan('12920_2012_343_MOESM3_ESM_entrez.txt', sep = '\t', skip=1)
t1.1 = scan('12920_2012_343_MOESM4_ESM_entrez.txt', sep = '\t', skip=1)
t1.2 = scan('12920_2012_343_MOESM5_ESM_entrez.txt', sep = '\t', skip=1)
t1.3 = scan('12920_2012_343_MOESM6_ESM_entrez.txt', sep = '\t', skip=1)
t2.1 = scan('12920_2012_343_MOESM7_ESM_entrez.txt', sep = '\t', skip=1)
t2.2 = scan('12920_2012_343_MOESM8_ESM_entrez.txt', sep = '\t', skip=1)

inmf_signatures_entrez = list()
inmf_signatures_entrez[['1']] = t1
inmf_signatures_entrez[['2']] = t2
inmf_signatures_entrez[['1.1']] = t1.1
inmf_signatures_entrez[['1.2']] = t1.2
inmf_signatures_entrez[['1.3']] = t1.3
inmf_signatures_entrez[['2.1']] = t2.1
inmf_signatures_entrez[['2.2']] = t2.2

library(devtools)
devtools::use_data(inmf_signatures_entrez, overwrite=TRUE)
file.rename('inmf_signatures_entrez.rda', '../../data/inmf_signatures_entrez.rda')

remove(inmf_signatures_entrez, t1, t2, t1.1, t1.2, t1.3, t2.1, t2.2)


# iNMF gene signature stored as gene symbols

load("inmf_signatures_sym.rdata")

inmf_signatures_sym = list()
inmf_signatures_sym[['1']] = ts.az.clust1.sym
inmf_signatures_sym[['2']] = ts.az.clust2.sym
inmf_signatures_sym[['1.1']] = ts.az.clust1.1.sym
inmf_signatures_sym[['1.2']] = ts.az.clust1.2.sym
inmf_signatures_sym[['1.3']] = ts.az.clust1.3.sym
inmf_signatures_sym[['2.1']] = ts.az.clust2.1.sym
inmf_signatures_sym[['2.2']] = ts.az.clust2.2.sym


devtools::use_data(inmf_signatures_sym, overwrite=TRUE)
file.rename('inmf_signatures_sym.rda', '../../data/inmf_signatures_sym.rda')


remove(inmf_signatures_sym, ts.az.clust1.sym, ts.az.clust2.sym, ts.az.clust1.1.sym, ts.az.clust1.2.sym, 
		 ts.az.clust1.3.sym, ts.az.clust2.1.sym, ts.az.clust2.2.sym)
}


if (FALSE){
	# compare the Schlicker symbols to a current (2018) translation of entrez to sym
	# save the current symbols (corresponding to entrez id's)

	load('../../data/inmf_signatures_entrez.rda')  # list of size 5 "inmf_signatures_entrez"
	load('inmf_signatures_sym.rdata')

	library(biomaRt, quietly=TRUE)
	ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
	map1 = getBM(attributes=c('hgnc_symbol', 'entrezgene'), filters = 'entrezgene', values = inmf_signatures_entrez$`1`, mart = ensembl)
	print( length(intersect(map1[,1], ts.az.clust1.sym)) / length( union(map1[,1], ts.az.clust1.sym) ) )
	#all.equal( sort(map1[,1]), sort(ts.az.clust1.sym) )
	
	map2 = getBM(attributes=c('hgnc_symbol', 'entrezgene'), filters = 'entrezgene', values = inmf_signatures_entrez$`2`, mart = ensembl)
	print( length(intersect(map2[,1], ts.az.clust2.sym)) / length( union(map2[,1], ts.az.clust2.sym) ) )
	#all.equal( sort(map2[,1]), sort(ts.az.clust2.sym) )

	map1.1 = getBM(attributes=c('hgnc_symbol', 'entrezgene'), filters = 'entrezgene', values = inmf_signatures_entrez$`1.1`, mart = ensembl)
	print( length(intersect(map1.1[,1], ts.az.clust1.1.sym)) / length( union(map1.1[,1], ts.az.clust1.1.sym) ) )
	#all.equal( sort(map1.1[,1]), sort(ts.az.clust1.1.sym) )

	map1.2 = getBM(attributes=c('hgnc_symbol', 'entrezgene'), filters = 'entrezgene', values = inmf_signatures_entrez$`1.2`, mart = ensembl)
	print( length(intersect(map1.2[,1], ts.az.clust1.2.sym)) / length( union(map1.2[,1], ts.az.clust1.2.sym) ) )
	#all.equal( sort(map1.2[,1]), sort(ts.az.clust1.2.sym) )

	map1.3 = getBM(attributes=c('hgnc_symbol', 'entrezgene'), filters = 'entrezgene', values = inmf_signatures_entrez$`1.3`, mart = ensembl)
	print( length(intersect(map1.3[,1], ts.az.clust1.3.sym)) / length( union(map1.3[,1], ts.az.clust1.3.sym) ) )
	#all.equal( sort(map1.3[,1]), sort(ts.az.clust1.3.sym) )

	map2.1 = getBM(attributes=c('hgnc_symbol', 'entrezgene'), filters = 'entrezgene', values = inmf_signatures_entrez$`2.1`, mart = ensembl)
	print( length(intersect(map2.1[,1], ts.az.clust2.1.sym)) / length( union(map2.1[,1], ts.az.clust2.1.sym) ) )
	#all.equal( sort(map2.1[,1]), sort(ts.az.clust2.1.sym) )

	map2.2 = getBM(attributes=c('hgnc_symbol', 'entrezgene'), filters = 'entrezgene', values = inmf_signatures_entrez$`2.2`, mart = ensembl)
	print( length(intersect(map2.2[,1], ts.az.clust2.2.sym)) / length( union(map2.2[,1], ts.az.clust2.2.sym) ) )
	#all.equal( sort(map2.2[,1]), sort(ts.az.clust2.2.sym) )

	sym = list('1'=map1, '2'=map2, '1.1'=map1.1, '1.2'=map1.2, '1.3'=map1.3, '2.1'=map2.1, '2.2'=map2.2)
	inmf_signatures_sym_from_entrez = sym
	save(inmf_signatures_sym_from_entrez, file = 'inmf_signatures_sym_from_entrez.rdata')
	save(inmf_signatures_sym_from_entrez, file = '../../data/inmf_signatures_sym_from_entrez.rda')

	rm(sym, map1, map2, map1.1, map1.2, map1.3, map2.1, map2.2, ensembl)
	rm(ts.az.clust1.sym, ts.az.clust2.sym, ts.az.clust1.1.sym, ts.az.clust1.2.sym, ts.az.clust1.3.sym, ts.az.clust2.1.sym, ts.az.clust2.2.sym)
}


if (FALSE){
# Sadanandam gene signature, and centroids

t = read.delim('Sadanandam786genes.csv', sep = ';', skip=6, header = TRUE)
rownames(t) = t[, 'Genes']
t = t[, -1]

sadanandam_786_genes = t
devtools::use_data(sadanandam_786_genes, overwrite=TRUE)
file.rename('sadanandam_786_genes.rda', '../../data/sadanandam_786_genes.rda')

remove(t, sadanandam_786_genes)

stopifnot(length(ls()) == 0)

}