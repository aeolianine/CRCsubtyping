rm(list = ls())

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


# Sadanandam gene signature, and centroids

t = read.delim('Sadanandam786genes.csv', sep = ';', skip=6, header = TRUE)
rownames(t) = t[, 'Genes']
t = t[, -1]

sadanandam_786_genes = t
devtools::use_data(sadanandam_786_genes, overwrite=TRUE)
file.rename('sadanandam_786_genes.rda', '../../data/sadanandam_786_genes.rda')

remove(t, sadanandam_786_genes)

stopifnot(length(ls()) == 0)