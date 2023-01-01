context('Normalization tools')

test_that('combat_to_ref() works as intended.', {

  expect_error(combat_to_ref())
  expect_error(combat_to_ref(dhdhnd))
  expect_error(combat_to_ref('fgdhdhnd'))
  expect_error(combat_to_ref(c(112,2,5,7,3)))

  mat = get_TCGA_synapse()
  out = cmdscale(as.dist(1-cor(mat, method = 'spearman')))
  group1 = rownames(out)[out[,1] > 0.015]
  group2 = rownames(out)[out[,1] <= -0.015]

  matNorm = combat_to_ref(mat[,group1], mat[,group2])
  expect_equal(dim(matNorm), dim(mat[,c(group1,group2)]))
  expect_equal(class(matNorm)[1], 'matrix')

  # clean up
  rm(mat, out, group1, group2, matNorm)
})
