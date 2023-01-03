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


test_that('compare_RF_to_SSP() works correctly', {

  expect_error( compare_RF_to_SSP() )
  expect_error( compare_RF_to_SSP('rgsrhsr') )
  expect_error( compare_RF_to_SSP(seghdjyk) )
  expect_error( compare_RF_to_SSP(c(1,2,4,6)) )

  mat = get_TCGA_synapse()
  resRF = subtypeCMS.RF(mat)
  resSSP = subtypeCMS.SSP(mat)

  expect_error( compare_RF_to_SSP(mat) )
  expect_error( compare_RF_to_SSP(resRF) )

  #expect_message(compare_RF_to_SSP(resRF, resSSP))
  out = compare_RF_to_SSP(resRF, resSSP)
  expect_equal(dim(out), c(170, 20))
  expect_true('RF.CMS1.posteriorProb' %in% colnames(out))
  expect_true('SSP.median.corToCMS4' %in% colnames(out))
  expect_true('RF.nearestCMS' %in% colnames(out))
  expect_true('SSP.predictedCMS' %in% colnames(out))

  # clean up
  rm(mat, resRF, resSSP, out)
})


test_that('compare_to_training_labels() works correctly', {

  expect_error(compare_to_training_labels())
  expect_error(compare_to_training_labels(dhdhd))
  expect_error(compare_to_training_labels('weghjyt'))
  expect_error(compare_to_training_labels(c(1,2,4,56)))
  expect_error(compare_to_training_labels(cbind(c(1,2,4,56), c(4,6,8,3))))
  expect_error(compare_to_training_labels(cbind(c(1,2,4,56), c(4,6,8,3)), whichDataset = 'TCGA'))

  mat = get_TCGA_synapse()
  resRF = subtypeCMS.RF(mat)

  expect_error(compare_to_training_labels(resRF))
  expect_error(compare_to_training_labels(resRF, whichDataset = 'some other dataset'))
  out = compare_to_training_labels(resRF, whichDataset = 'TCGA')

  expect_equal(class(out), 'data.frame')
  expect_equal(ncol(out), 7)
  expect_true(all(colnames(out) == c('training label', 'RF.nearestCMS', 'RF.predictedCMS',
                                   'RF.CMS1.posteriorProb', 'RF.CMS2.posteriorProb',
                                   'RF.CMS3.posteriorProb', 'RF.CMS4.posteriorProb')))
  # clean up
  rm(mat, resRF)

})
