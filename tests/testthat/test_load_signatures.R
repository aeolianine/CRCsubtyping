context('Subtype signatures')

test_that('The Schlicker signature was loaded properly.', {

  expect_equal(names(load_schlicker_signature()), c("1","2","1.1","1.2","1.3","2.1","2.2"), all = TRUE)
  expect_equal(length(load_schlicker_signature()$`1`), 1344)
  expect_equal(length(load_schlicker_signature()$`2`), 410)
  expect_equal(length(load_schlicker_signature()$`1.1`), 284)
  expect_equal(length(load_schlicker_signature()$`1.2`), 139)
  expect_equal(length(load_schlicker_signature()$`1.3`), 216)
  expect_equal(length(load_schlicker_signature()$`2.1`), 205)
  expect_equal(length(load_schlicker_signature()$`2.2`), 198)

  expect_equal(class(load_schlicker_signature()$`2.2`), 'character')

  expect_true("LOC100505806" %in% load_schlicker_signature()$`1.3`)

  expect_true("TWIST2" %in% load_schlicker_signature()$'1')
  expect_true("CORO2A" %in% load_schlicker_signature()$'2')
  expect_true("C20orf194" %in% load_schlicker_signature()$'1.1')
  expect_true("TNFRSF10A" %in% load_schlicker_signature()$'1.2')
  expect_true("KRTAP4-1" %in% load_schlicker_signature()$'1.3')
  expect_true("LOC100506100" %in% load_schlicker_signature()$'2.1')
  expect_true("KIAA0226L" %in% load_schlicker_signature()$'2.2')


  expect_error(load_schlicker_signature(noInput))
  expect_error(load_schlicker_signature(annotation = 'sgljshdgsj'))
  expect_equal(names(load_schlicker_signature(annotation = 'entrez')), c("1","2","1.1","1.2","1.3","2.1","2.2"), all = TRUE)

  expect_equal(class(load_schlicker_signature()), 'list')
})

test_that('The Sadanandam signature was loaded properly.', {

  expect_equal( class(load_sadanandam_signature()), 'data.frame' )

  expect_equal(colnames(load_sadanandam_signature()), c("Inflammatory", "Goblet.like", "Enterocyte", "TA", "Stem.like"), all = TRUE)

  expect_equal(nrow(load_sadanandam_signature()), 786)
  expect_equal(ncol(load_sadanandam_signature()), 5)
  expect_true('RORA' %in% rownames(load_sadanandam_signature()))
  expect_true('LRRC19' %in% rownames(load_sadanandam_signature()))
  expect_true('LOC729680' %in% rownames(load_sadanandam_signature()))
  expect_true('PPAPDC1A' %in% rownames(load_sadanandam_signature()))
  expect_true('HLA-DPA1' %in% rownames(load_sadanandam_signature()))
  expect_true('ST6GALNAC1' %in% rownames(load_sadanandam_signature()))

  expect_error(load_sadanandam_signature(dfkjghdfkj))
  expect_error(load_sadanandam_signature('dfkjghdhffkj'))

})


test_that('CMS genes are loaded properly.', {

  expect_error(load_CMS_genes(dfkhjdlfk))
  expect_error(load_CMS_genes('dfkhjdlfk'))

  expect_equal(class(load_CMS_genes())[1], 'matrix')
  expect_equal(ncol(load_CMS_genes()), 6)
  expect_equal(nrow(load_CMS_genes()), 273)

  expect_true("430" %in% rownames(load_CMS_genes()))
  expect_true("7078" %in% rownames(load_CMS_genes()))
  expect_true("4644" %in% rownames(load_CMS_genes()))

  expect_equal(colnames(load_CMS_genes()), c('CMS1', 'CMS2', 'CMS3', 'CMS4', 'MeanDecreaseAccuracy', 'MeanDecreaseGini'))

  genes = load_CMS_genes(geneSymbol = TRUE)
  expect_equal(class(genes)[1], 'matrix')
  expect_equal(ncol(genes), 6)
  expect_equal(nrow(genes), 273)
  expect('SULF1' %in% rownames(genes))
  expect('CDH11' %in% rownames(genes))
  expect('OGN' %in% rownames(genes))

  expect_equal(colnames(genes), c('CMS1', 'CMS2', 'CMS3', 'CMS4', 'MeanDecreaseAccuracy', 'MeanDecreaseGini'))

  remove(genes)
})

test_that('Signatures were intersected properly.', {

  expect_error(intersect_signature(1))
  expect_error(intersect_signature(1:2,2:4))
  expect_equal(length(intersect_signature(1:2,'a', verbose=FALSE)), 0)
  expect_equal(intersect_signature(as.character(1:2), as.character(2:4)), '2')

  expect_equal(intersect_signature(c('a','b','c'), list('1'=c('a','b','d','e'),'2'=c('0','b','d','f'))),
                                                  list('1'=c('a','b'),'2'=c('b')))

  expect_equal( class(intersect_signature('TMEM37', rownames(load_sadanandam_signature()), verbose = TRUE)), class(rownames(load_sadanandam_signature())) )
  expect_equal(length(intersect_signature('TMEM37', rownames(load_sadanandam_signature()), verbose = TRUE)), 1)


  sig = load_schlicker_signature()
  genes = c("SLC16A6",  "ST8SIA4",  "APBB1IP",  "KGFLP2",      # 1
            "PRR11", "CASC5", "SHROOM3", "KIAA1804", "RORC",   # 2
            "PDE5A", "MIR143HG", "ARHGEF26", "TSPAN2", "FAM126A",  # 1.1
            "CDCP1", "DDX58", "USP18", "HPSE", "SIGLEC1", "OGFRL1",   # 1.2
            "LOC100507312", "PLEKHG4", "NFYA", "FOXA3", "RPS6KA6")   # 1.3
  newsig = intersect_signature(genes, sig, verbose = FALSE)

  expect_equal( class(newsig), class(sig) )
  expect_equal( length(newsig), length(sig) )
  expect_equal( names(newsig), names(sig) )

  expect_true(length(newsig$'1')==11)
  expect_true(length(newsig$'2')==6)
  expect_true(length(newsig$'1.1')==5)
  expect_true(length(newsig$'1.2')==6)
  expect_true(length(newsig$'1.3')==5)

  expect_true(length(newsig$'2.1')==1)
  expect_true(length(newsig$'2.2')==1)


  remove(sig, newsig, genes)
})


test_that('Entrez to symbol conversion works.', {

  genesInEntrez = c("5350", "7070", "196051", "1000", "1959", "1289", "25937", "5118", "1290", "59272")
  genesInSymbol = c("PLN", "THY1", "PLPP4", "CDH2", "EGR2", "COL5A1", "WWTR1", "PCOLCE", "COL5A2", "ACE2")
  names(genesInSymbol) = genesInEntrez

  expect_equal(entrez2symbol(genesInEntrez), genesInSymbol)
  expect_equal(class(entrez2symbol(c("5350", "7070", "196051"))), 'character')

  expect_error(entrez2symbol(dfhbdfhbdh))
  expect_true(is.na(entrez2symbol('dfhbdfhbdh')))

  remove(genesInEntrez, genesInSymbol)

})


test_that('Symbol to entrez conversion works.', {

  expect_equal(length(symbol2entrez(c('APC', 'KRAS', 'TP53', 'SMAD4', 'TGFB1', 'TGFB2', 'TGFB3', 'SMAD7'))), 8)
  expect_equal(class(symbol2entrez(c('APC', 'KRAS', 'TP53', 'SMAD4', 'TGFB1', 'TGFB2', 'TGFB3', 'SMAD7'))), 'character')

  gene = 'APCKRASTP53SMAD4'
  names(gene) = 'APCKRASTP53SMAD4'
  expect_true(is.na(symbol2entrez(c('APCKRASTP53SMAD4'))))

  result = c("324", "3845", "7157", "4089", "7040", "7042", "7043", "4092")
  names(result) = c('APC', 'KRAS', 'TP53', 'SMAD4', 'TGFB1', 'TGFB2', 'TGFB3', 'SMAD7')
  expect_equal(symbol2entrez(c('APC', 'KRAS', 'TP53', 'SMAD4', 'TGFB1', 'TGFB2', 'TGFB3', 'SMAD7')), result)

  remove(gene, result)
})
