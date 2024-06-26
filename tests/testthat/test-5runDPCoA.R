context("addDPCoA")
test_that("addDPCoA", {
    skip_if_not(require("ade4", quietly = TRUE))
    data(esophagus, package="mia")
    #
    esophagus <- addDPCoA(esophagus)
    expect_named(reducedDims(esophagus),"DPCoA")
    expect_true(is.matrix(reducedDim(esophagus,"DPCoA")))
    expect_equal(dim(reducedDim(esophagus,"DPCoA")),c(3,2))
    red <- reducedDim(esophagus,"DPCoA")
    expect_equal(names(attributes(red)),
                 c("dim","dimnames","eig","sample_red","feature_weights",
                   "sample_weights" ))
    expect_equal(dim(attr(red,"sample_red")),c(58,2))
    
    # ERRORs
    expect_error(
        addDPCoA(esophagus, assay.type = "test", tree.name = "phylo", ncomponents = 2, ntop = NULL,
                 subset.row = NULL, scale = FALSE, transposed = FALSE)
    )
    expect_error(
        addDPCoA(esophagus, assay.type = 1, tree.name = "phylo", ncomponents = 2, ntop = NULL,
                 subset.row = NULL, scale = FALSE, transposed = FALSE)
    )
    expect_error(
        addDPCoA(esophagus, assay.type = TRUE, tree.name = "phylo", ncomponents = 2, ntop = NULL,
                 subset.row = NULL, scale = FALSE, transposed = FALSE)
    )
    expect_error(
        addDPCoA(esophagus, assay.type = "counts", tree.name = "test", ncomponents = 2, ntop = NULL,
                 subset.row = NULL, scale = FALSE, transposed = FALSE)
    )
    expect_error(
        addDPCoA(esophagus, assay.type = "counts", tree.name = 1, ncomponents = 2, ntop = NULL,
                 subset.row = NULL, scale = FALSE, transposed = FALSE)
    )
    expect_error(
        addDPCoA(esophagus, assay.type = "counts", tree.name = "phylo", ncomponents = TRUE, ntop = NULL,
                 subset.row = NULL, scale = FALSE, transposed = FALSE)
    )
    expect_error(
        addDPCoA(esophagus, assay.type = "counts", tree.name = "phylo", ncomponents = "test", ntop = NULL,
                 subset.row = NULL, scale = FALSE, transposed = FALSE)
    )
    expect_error(
        addDPCoA(esophagus, assay.type = "counts", tree.name = "phylo", ncomponents = "test", ntop = "test",
                 subset.row = NULL, scale = FALSE, transposed = FALSE)
    )
    expect_error(
        addDPCoA(esophagus, assay.type = "counts", tree.name = "phylo", ncomponents = 1.3, ntop = "test",
                 subset.row = NULL, scale = FALSE, transposed = FALSE)
    )
    expect_error(
        addDPCoA(esophagus, name = c("test", "test2"), assay.type = "counts",
                 tree.name = "phylo", ncomponents = 1.3, ntop = "test",
                 subset.row = NULL, scale = FALSE, transposed = FALSE)
    )
    expect_error(
        addDPCoA(esophagus, name = 1, assay.type = "counts", tree.name = "phylo", ncomponents = 1.3, ntop = "test",
                 subset.row = NULL, scale = FALSE, transposed = FALSE)
    )
    
    data(GlobalPatterns, package="mia")
    tse <- mergeSEs(esophagus, GlobalPatterns)
    # expect_warning(addDPCoA(tse))
    # expect_warning(addDPCoA(tse, tree.name = "phylo.1"))
})
