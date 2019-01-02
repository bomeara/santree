test_that("conversion from phylo works", {
  phy <- ape::rcoal(5)
  output <- convert_phylo_to_sankey(phy)
  expect_true(inherits(output, "data.frame"))
})

test_that("get_diversity_for_lineage works", {
  diversities <- get_diversity_for_lineage("Dinosauria")
  expect_true(nrow(diversities)>5)
  expect_true(max(diversities$implied_in_bin)>5)
  diversities <- get_diversity_for_lineage("Unicorns")
  expect_true(is.null(diversities))
})
