context('Downloading and transforming data.')

test_that('getSynapseTable works properly.', {

	if (file.exists('~/.synapseConfig')){

		expect_error(getSynapseTable(ljndfjn))   # non-existant input
		expect_error(getSynapseTable())          # no input
		# expect_message(getSynapseTable('syn2165691')): it does throw a message, but expect_message does not see it
	}


})


test_that('get_GSE33113_synapse works properly.', {

	out = get_GSE33113_synapse()
	expect_equal(class(out), 'data.frame')
	expect_equal(dim(out), c(54675, 90))
	expect_true( all(grepl('at', rownames(out))) )

	remove(out)
})


test_that('get_GSE33113_GEO works properly.', {
})


test_that('get_TCGA_synapse works properly.', {

	out = get_TCGA_synapse()
	expect_equal(class(out), 'data.frame')
	expect_equal(dim(out), c(20293, 577))
	expect_true('KRAS' %in% rownames(out))
	expect_true('GREM1' %in% rownames(out))
	expect_true( all(grepl('TCGA', colnames(out))) )

	remove(out)

})

