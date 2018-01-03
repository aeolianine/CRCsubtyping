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

test_that('', {

  plotSubtypesCMS.RF


})
