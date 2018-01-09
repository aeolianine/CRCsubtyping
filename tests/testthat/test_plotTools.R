context('Plotting tools')

test_that('CMS colors have been defined properly.', {
  colors = CMScolors()

  expect_true( all( colors == c('orange3','blue2','lightcoral','lightgreen') ) )
  expect_true( all( names(colors) == c('CMS1', 'CMS2', 'CMS3', 'CMS4') ) )

})

test_that('The Sadanandam centroids are plotted properly.', {
  expect_error( plotSadanandamSubtypeCentroids('lkhedklhb') )
  pdf('test.pdf')
  out = plotSadanandamSubtypeCentroids()
  dev.off()
  expect_type(out, 'list')
  expect_true('colDendrogram' %in% names(out))
  expect_true('call' %in% names(out))

  remove(out)
  file.remove('test.pdf')
})

test_that('plotSubtypesCMS.RF works.', {

  expect_error(plotSubtypesCMS.RF())
  expect_error(plotSubtypesCMS.RF(dfhdfhn))
  expect_error(plotSubtypesCMS.RF('hdh'))
  expect_error(plotSubtypesCMS.RF(c(1,2,33,6)))
  expect_error(plotSubtypesCMS.RF( cbind(c(1,2,33,6), c(3,4,5,2)) ))

  mat = get_TCGA_synapse()
  res = subtypeCMS.RF(mat)
  png('test.png')
  plotSubtypesCMS.RF(res)
  dev.off()

  expect_true(file.exists('test.png'))

  # clean up
  rm(mat, res)
  file.remove('test.png')

})


test_that('plotSubtypesCMS.SSP works.', {

  expect_error(plotSubtypesCMS.SSP())
  expect_error(plotSubtypesCMS.SSP(dfhdfhn))
  expect_error(plotSubtypesCMS.SSP('hdh'))
  expect_error(plotSubtypesCMS.SSP(c(1,2,33,6)))
  expect_error(plotSubtypesCMS.SSP( cbind(c(1,2,33,6), c(3,4,5,2)) ))

  mat = get_TCGA_synapse()
  res = subtypeCMS.SSP(mat)
  png('test.png')
  plotSubtypesCMS.SSP(res)
  dev.off()

  expect_true(file.exists('test.png'))

  # clean up
  rm(mat, res)
  file.remove('test.png')

})


test_that('plot_mCombat_effect() works.', {

  n = 30

  mat1 = rbind(runif(n), runif(n), runif(n))
  colnames(mat1) = paste0('mat1-',1:n)
  mat2 = rbind(runif(n), runif(n), runif(n)) + 11
  colnames(mat2) = paste0('mat2-',1:n)

  rownames(mat1) = c('geneA', 'geneB', 'geneC')
  rownames(mat2) = rownames(mat1)

  matNorm = combatToRef(mat1, mat2)
  png('test.png')
  plot_mCombat_effect(mat1, mat2, matNorm)
  dev.off()

  expect_true(file.exists('test.png'))
  file.remove('test.png')

  # clean up
  rm(n, mat1, mat2, matNorm)
})
