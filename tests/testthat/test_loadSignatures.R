context('Subtype signatures')

test_that('The Schlicker signature was loaded properly.', {

  expect_equal(names(loadSchlickerSignature()), c("1","2","1.1","1.2","1.3","2.1","2.2"), all = TRUE)
  expect_equal(length(loadSchlickerSignature()$`1`), 1344)
  expect_equal(length(loadSchlickerSignature()$`2`), 410)
  expect_equal(length(loadSchlickerSignature()$`1.1`), 284)
  expect_equal(length(loadSchlickerSignature()$`1.2`), 139)
  expect_equal(length(loadSchlickerSignature()$`1.3`), 216)

  expect_true("LOC100505806" %in% loadSchlickerSignature()$`1.3`)


  expect_error(loadSchlickerSignature(noInput))
  expect_error(loadSchlickerSignature(annotation = 'sgljshdgsj'))
  expect_equal(names(loadSchlickerSignature(annotation = 'entrez')), c("1","2","1.1","1.2","1.3","2.1","2.2"), all = TRUE)

  expect_equal(class(loadSchlickerSignature()), 'list')
})

test_that('The Sadanandam signature was loaded properly.', {

  expect_equal( class(loadSadanandamSignature()), 'data.frame' )

  expect_equal(colnames(loadSadanandamSignature()), c("Inflammatory", "Goblet.like", "Enterocyte", "TA", "Stem.like"), all = TRUE)
  expect_equal(nrow(loadSadanandamSignature()), 786)
  expect_true('RORA' %in% rownames(loadSadanandamSignature()))

  expect_error(loadSadanandamSignature(dfkjghdfkj))
  expect_error(loadSadanandamSignature('dfkjghdhffkj'))

})



