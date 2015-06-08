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

test_that("monophyly is correctly evaluated",{
	expect_equal( is_monophyletic( siphonophore_ml, c("Hippopodius_hippopus_Pacific","Hippopodius_hippopus_Atlantic") ), TRUE )
	expect_equal( is_monophyletic( siphonophore_ml, c("Hippopodius_hippopus_Pacific","Clausophyes_ovata") ), FALSE )
	expect_equal( is_monophyletic( siphonophore_ml, c("Diphyes_dispar","Abylopsis_tetragona","Chelophyes_appendiculata") ), TRUE ) # Check across root to make sure evaluation is unrooted
})

cut_nodes = c( 74, 98 ) # 74 is Forskalia, 98 is Cystonectae
subtrees = decompose( siphonophore_ml, cut_nodes )
test_that("tree is successfully decomposed",{
	expect_equal( length( subtrees ), 3 )
})