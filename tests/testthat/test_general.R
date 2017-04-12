# Per http://r-pkgs.had.co.nz/tests.html

library(hutan)
library(ape)
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
subtrees = decompose_tree( siphonophore_ml, cut_nodes )
test_that("tree is successfully decomposed",{
	expect_equal( length( subtrees ), 3 )
})

# Create tree
treetext = "((Prayidae_D27SS7@1152004:0.310638840467,(Hippopodius_hippopus@867239:0.242845544461,(Craseoa_lathetica@180617:1.21823973672E-6,Prayidae_D27D2@922608:0.013252166117)n38047893-N-100:0.143340517912)n38047889-N-83:0.0697164264107)n38047887-Y-100:0.155909099947,(Kephyes_ovata@1402412:0.0624504448762,Chuniphyes_multidentata@278022:0.247764507969)n38047885-N-100:0.155909099947)n38048123-N;"
gene_tree = read.tree( text=treetext )

# Identify duplicated nodes 
duplications = grep( "-Y", gene_tree$node.label ) + length( gene_tree$tip.label )

# Get their immediate descendants, which define the clades we want to excise
to_prune = gene_tree$edge[,2][ gene_tree$edge[,1] %in% duplications ]

subgenetrees = decompose_tree( gene_tree, to_prune )

test_that("tree with node names is successfully decomposed",{
	expect_equal( length( subgenetrees ), 3 )
})

test_that("can calculate difference_from_calibrated",{
	expect_true( abs(difference_from_calibrated( gene_tree, quiet=TRUE ) - 0.3048597 ) < 1e-06 )
})

test_that("result of difference_from_calibrated is 0 on already calibrated tree",{
	expect_true( abs(difference_from_calibrated( chronos(gene_tree), quiet=TRUE ) ) < 1e-02 )
})

test_that("the correct number of depths are produced",{
	depths = distance_from_tip(gene_tree)
	expect_equal( length(depths), 11 )
})

test_that("can get the edges that connect a tip to the root",{
	edges = ancestral_edges( siphonophore_ml, 2 )
	expect_true( setequal( edges, c( 3, 2 ) ) )
})

test_that("can get the edges that connect a tip to the root",{
	edges = ancestral_edges( siphonophore_ml, 10 )
	expect_true( setequal( edges, c( 30, 29, 28, 26, 24, 23, 22, 12, 11, 10, 9, 8, 7, 6, 5, 4, 2 ) ) )
})

test_that("can get the edges that connect an internal node to the root",{
	edges = ancestral_edges( siphonophore_ml, 65 )
	expect_true( setequal( edges, c( 10, 9, 8, 7, 6, 5, 4, 2 ) ) )
})

test_that("can get the connecting edges for sister species",{
	edges = connecting_edges( siphonophore_ml, 31, 32)
	expect_true( setequal( edges, c( 71, 72 ) ) )
})

test_that("can get the connecting edges for tips that span root",{
	edges = connecting_edges( siphonophore_ml, 1, 32)
	expect_true( setequal( edges, c( 1, 72, 70, 69, 67, 66, 58, 57, 11, 10, 9, 8, 7, 6, 5, 4, 2 ) ) )
})

test_that("can get the connecting edges for ancestor and descendent nodes",{
	edges = connecting_edges( siphonophore_ml, 65, 67)
	expect_true( setequal( edges, c( 11, 12 ) ) )
})