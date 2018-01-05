context('CMS subtyping functions.')

test_that('subtypeCMS.RF() works correctly', {

  expect_error(subtypeCMS.RF('dhddrh'))
  expect_error(subtypeCMS.RF(dhddrh))
  expect_error(subtypeCMS.RF(c(1,1,3)))
  expect_error(subtypeCMS.RF())

  mat = get_TCGA_synapse()
  res = subtypeCMS.RF(mat)
  expect_equal(class(res), 'data.frame')
  expect_equal(ncol(res), 6)
  expect_equal(colnames(res), c("RF.CMS1.posteriorProb", "RF.CMS2.posteriorProb", "RF.CMS3.posteriorProb",
                                "RF.CMS4.posteriorProb", "RF.nearestCMS", "RF.predictedCMS" ))
  expect_equal(ncol(mat), nrow(res))

  # clean up
  rm(mat, res)
})



test_that('subtypeCMS.SSP() works correctly', {

  expect_error(subtypeCMS.SSP('dhddrh'))
  expect_error(subtypeCMS.SSP(dhddrh))
  expect_error(subtypeCMS.SSP(c(1,1,3)))
  expect_error(subtypeCMS.SSP())

  mat = get_TCGA_synapse()
  res = subtypeCMS.SSP(mat)
  expect_equal(class(res), 'data.frame')
  expect_equal(ncol(res), 14)
  expect_equal(colnames(res), c("SSP.min.corToCMS1", "SSP.min.corToCMS2", "SSP.min.corToCMS3", "SSP.min.corToCMS4",
                                "SSP.median.corToCMS1","SSP.median.corToCMS2", "SSP.median.corToCMS3", "SSP.median.corToCMS4",
                                "SSP.max.corToCMS1", "SSP.max.corToCMS2", "SSP.max.corToCMS3", "SSP.max.corToCMS4",
                                "SSP.nearestCMS", "SSP.predictedCMS"))
  expect_equal(ncol(mat), nrow(res))
  # clean up
  rm(mat, res)

})
