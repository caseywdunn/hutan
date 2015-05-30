# Per http://r-pkgs.had.co.nz/tests.html

library(hutan)
context("general")

data( siphonophore_ml )
data( siphonophore_constraint )

zc <- zero_constrained( siphonophore_ml, siphonophore_constraint )
edge_diff = siphonophore_ml$edge.length - zc$edge.length
mod_edges = which( edge_diff > 0 )

test_that("correct edges are modified by zero_constrained", {
  expect_equal( mod_edges, c(12, 22, 57) )
})
