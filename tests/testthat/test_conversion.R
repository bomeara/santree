test_that("conversion from phylo works", {
  phy <- ape::rcoal(5)
  output <- convert_phylo_to_sankey(phy)
  expect_true(inherits(output, "data.frame"))
})
