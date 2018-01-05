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
