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
