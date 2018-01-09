context('Normalization tools')

test_that('combatToRef works as intended.', {

  expect_error(combatToRef())
  expect_error(combatToRef(dhdhnd))
  expect_error(combatToRef('fgdhdhnd'))
  expect_error(combatToRef(c(112,2,5,7,3)))

  mat = get_TCGA_synapse()
  out = cmdscale(as.dist(1-cor(mat, method = 'spearman')))
  group1 = rownames(out)[out[,1] > 0.015]
  group2 = rownames(out)[out[,1] <= -0.015]

  matNorm = combatToRef(mat[,group1], mat[,group2])
  expect_equal(dim(matNorm), dim(mat[,c(group1,group2)]))
  expect_equal(class(matNorm), 'matrix')

  # clean up
  rm(mat, out, group1, group2, matNorm)
})
