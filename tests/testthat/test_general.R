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
subtrees = decompose( siphonophore_ml, cut_nodes )
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

#subgenetrees = decompose( gene_tree, to_prune )

#test_that("tree with node names is successfully decomposed",{
#	expect_equal( length( subgenetrees ), 3 )
#})