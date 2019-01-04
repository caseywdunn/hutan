#' Decomposes a single tree into a series of subtrees designated by internal
#' node numbers
#'
#' @param phy The tree to be decomposed, as an ape phylo object
#' @param x A vector of internal node numbers. The tree phy will be cut on each branch
#' that subtends each of these nodes.
#' @return A list of phylo objects
#' @export
decompose_tree <- function( phy, x ){

	# Create a vactor that indicates which subtree each tip will be in
	partitions = rep( "r", length(phy$tip.label) )

	for (node in x){
		labelstring = paste( "_", node, sep="")
		tips = tip_descendants( phy, node )
		partitions[ tips ] = paste( partitions[ tips ], labelstring, sep="" )
	}

	subtrees = lapply( unique( partitions ), function(partition) safe.drop.tip( phy, which( partitions != partition ) ) )

	return( subtrees )

}

#' Drops specified tips from a phylogeny. Like ape's drop.tip(), but it works when only a single tip is
#' to be retained.
#'
#' @param phy The tree, as an ape phylo object
#' @param tip A vector of tip numbers to be removed.
#' @return The reduced tree, as a phylo object
#' @export
safe.drop.tip <- function( phy, tip ){
	keep = 1:length( phy$tip.label )
	keep = keep[ ! keep %in% tip ]

	if( length( keep ) == 0 ){
		return( NULL )
	}
	else if ( length( keep ) == 1 ){
		# Create a one-taxon tree and return it
		text = paste( "(", phy$tip.label[ keep ], ");" )
		return( ape::read.tree( text=text ) )
	}
	else{
		return( ape::drop.tip( phy, tip ) )
	}
}

#' Cuts a single tree on the branch subtending a specified node
#'
#' @param phy The tree to be cut, as an ape phylo object
#' @param x An internal node number. The tree phy will be cut on the branch
#' that subtends this nodes.
#' @return A list of phylo objects that are the subtrees
#' @export
cut_tree <- function( phy, x ){
	# Get the number of the branch that subtends the node
	cut_branches = which( phy$edge[,2] == x )

	if ( length( cut_branches ) < 1 ){
		ape::plot.phylo(phy)
		ape::nodelabels()
		stop( paste ( "No edge subtending requested node! Requested node: ", x, sep='' ) )
	}

	if ( length( cut_branches ) > 1 ){
		stop( "More than one subtending branch!" )
	}

	tips1 = bipartition_for_edge_by_label( cut_branches[1], phy )
	tips2 = flip_bipartition( phy, tips1 )

	return( list( ape::drop.tip( phy, tips1 ), ape::drop.tip( phy, tips2 ) ) )

}

#' Generates a tree with a single resolved bipartition between two sets of
#' tip names. Useful for generating constraint trees.
#'
#' @param tips1 A vector of tip names for clade 1
#' @param tips2 A vector of tip names for clade 2
#' @return A phylo object
#' @examples
#' library( ape )
#' tips1 = c("a", "b", "c")
#' tips2 = c("d", "e", "f")
#' ctree = generate_constaint_tree( tips1, tips2 )
#' @export
generate_constaint_tree <- function( tips1, tips2 ){
	clade1 = paste( tips1, collapse=",")
	clade2 = paste( tips2, collapse=",")
	treetext = paste( c("((", clade1, "),(", clade2, "));"), collapse="")
	phy = read.tree( text=treetext )
	return( phy )
}

#' Generates the "zero-constrained" tree described by Susko 2014
#' (http://dx.doi.org/10.1093/molbev/msu039)
#'
#' @param phy_resolved A fully resolved phylogeny stored as a phylo object, e.g. an ML
#' tree.
#' @param phy_constraint A partially resolved constraint tree.
#' @param epsilon The value to replace the branch length with
#' @return A phylo object containing a tree that is the same as phy_resolved, except that
#' the length of edges that are incompatible with phy_constraint are replaced with epsilon.
#' @examples
#' data( siphonophore_ml )
#' data( siphonophore_constraint )
#' zc <- zero_constrained( siphonophore_ml, siphonophore_constraint )
#' @export
zero_constrained <- function ( phy_resolved, phy_constraint, epsilon=0.000001 ){
	phy_resolved$edge.length[ ! compatible_edges(  phy_resolved, phy_constraint ) ] <- epsilon
	return( phy_resolved )
}

#' Identify the edges in one phylo object that are compatible with the edges in another
#' phylo object. Requires the same tip labels for each tree.
#'
#' @param phy1 The tree under consideration
#' @param phy2 The tree to be compared to
#' @return A boolean vector corresponding to the edges in phy1. Each element is FALSE if
#' the edge is iscompatible with phy2, or TRUE if compatible.
#' @export
compatible_edges <- function( phy1, phy2 ){

	# First check to see that the two trees have the same tips
	if( setequal(phy1$tip.label, phy2$tip.label) == FALSE ) {
		stop("Trees do not have the same tip names.")
	}

	bi1 = get_bipartitions( phy1 )
	bi2 = get_bipartitions( phy2 )

	compat = lapply( bi1, is_compatible_with_set, bi_list=bi2, phy=phy1 )

	return ( unlist( compat ) )
}

#' Test if a set of tips, specified as a vector of tip labels, forms a monophyletic
#' group in a given tree. The test is unrooted, i.e. the group can span the root.
#'
#' @param phy The tree under consideration
#' @param x A vector of the labels of the tips in question
#' @return A boolean, TRUE if the tips form a monophyletic group.
#' @export
is_monophyletic <- function( phy, x ){
	bi_list = get_bipartitions( phy )
	return ( is_compatible_with_set( x, bi_list, phy ) )
}

#' Check if bipartition bi is compatible with the bipartitions in bi_list.
#' Each bipartition is defined as a vector of the names of the tips on one side of the
#' bipartition.
#'
#' @param bi The query bipartition.
#' @param bi_list A list of the bipartitions to be compared against.
#' @param phy A phylo object describing a tree that includes all tips under investigation.
#' This is used to infer the other half of each bipartition.
#' @return TRUE if bi is compatible with all bipartition in bi_list, otherwise FALSE.
#' @export
is_compatible_with_set <- function( bi, bi_list, phy ) {

	compatible <- lapply( bi_list, are_bipartitions_compatible, bi2=bi, phy=phy )
	return ( all( unlist( compatible ) ) )
}

#' Check if two bipartitions drawn from trees with the same tips are compatible with each other.
#' Each bipartition is defined as a vector of the names of the tips on one side of the
#' bipartition.
#'
#' @param bi1 The first bipartition.
#' @param bi2 The second bipartition.
#' @param phy A phylo object describing a tree that includes all tips under investigation.
#' This is used to infer the other half of each bipartition.
#' @return TRUE if bi1 is compatible with bi2, otherwise FALSE.
#' @export
are_bipartitions_compatible <- function( bi1, bi2, phy ) {

	labels = phy$tip.label

	# Take the left side of the bipartition, as given
	left1  = bi1

	# Take the right side of the bipartition as all taxa not in the left side
	right1 = labels[ ! labels %in% left1 ]

	# Do the same for the second bipartition
	left2  = bi2
	right2 = labels[ ! labels %in% left2 ]

	# Bipartition 1 is compatible with Bipartition 2 if either side of Bipartition 1
	# is a subset of either side of Bipartition 2

	left1_compat  = all( left1 %in% left2 ) | all( left1 %in% right2 )
	right1_compat = all( right1 %in% left2 ) | all( right1 %in% right2 )

	compatible = left1_compat | right1_compat

	return( compatible )
}

#' Get tips and labels of a phylo object.
#'
#' @param phy A phylo object.
#' @return A vector of all the tips, annotated with their names
#' @export
tips <- function(phy) {

	t <- phy$edge[ ! phy$edge[,2] %in% phy$edge[,1] ,2]
	t <- t[order(t)]

	names(t) <- phy$tip.label

	return(t)
}

#' Get all the descendants of a given node in a tree.
#'
#' @param phy A phylo object that specifies the tree.
#' @param a The number of a node in phy.
#' @param keep_node If FALSE, do not include a in the result.
#' @return A vector of nodes (specified by number) that are descendants of a. Includes
#' internal and tip nodes.
#' @export
descendants <- function(phy, a, keep_node=FALSE) {
	# returns a vector of all the descendants of node a, including tips and internal nodes
	# Based on a breadth-first search
	q <- c(a)
	d <- c()
	while (length(q)>0){
		# dequeue first item
		n <- q[1]
		q <- q[q != n]

		# add it to the descendants vector
		d <- append(d, n)

		# add it's descendants to the queue
		q <- append(q, phy$edge[phy$edge[,1] == n,2])
	}

	# Remove the original source node from the list
	if ( ! keep_node ){
		d <- d[d != a ]
	}
	return(d)
}

#' Get all the tips that are descendants of a given node in a tree.
#'
#' @param phy A phylo object that specifies the tree.
#' @param a The number of a node in phy.
#' @return A vector of tip nodes (specified by number) that are descendants of a. If a is
#' a tip, it is the sole element of this vector.
#' @export
tip_descendants <- function(phy, a) {
	# returns a vector of all the tips that are descendants of a
	t <- tips(phy)
	return( t[ t %in% descendants( phy, a, keep_node=TRUE ) ] )
}

#' Get a bipartition, described as a vector of tip numbers, from a specified tree and
#' edge number.
#'
#' @param phy A phylo object that specifies the tree.
#' @param edge The number of the edge that defines the bipartition.
#' @return A vector of tip nodes (specified by numbers) that define one half of the
#' bipartition (the other half is the set of tip nodes that are not in this vector).
#' @export
bipartition_for_edge <- function( phy, edge ){

	# Not certain which of the two nodes that make up the edge is ancestor and which is
	# descendant. Descendant will have fewer descendant tips, so pick the node with the
	# fewest descendants.

	left_node  <- phy$edge[edge,1]
	right_node <- phy$edge[edge,2]

	left_tips <- tip_descendants( phy, left_node )
	right_tips <- tip_descendants( phy, right_node )

	if ( length( left_tips ) < length( right_tips ) ){
		return( left_tips )
	}
	else {
		return( right_tips )
	}

}

#' Get a bipartition, described as a vector of tip labels, from a specified tree and
#' edge number.
#'
#' @param phy A phylo object that specifies the tree.
#' @param edge The number of the edge that defines the bipartition.
#' @return A vector of tip nodes (specified by labels) that define one half of the
#' bipartition (the other half is the set of tip nodes that are not in this vector).
#' @export
bipartition_for_edge_by_label <- function( edge, phy ){

	return( phy$tip.label[ bipartition_for_edge( phy, edge ) ]	 )

}

#' Given a tree and a bipartition, described as a vector of tip labels on one side of
#' of the bipartition, return the same bipartition but defined by the tip labels on
#' the other side of the bipartition.
#'
#' @param phy A phylo object that specifies the tree.
#' @param bi The bipartition.
#' @return A vector of tip nodes (specified by labels) that define one half of the
#' bipartition (the other half is the set of tip nodes that are provided as bi).
#' @export
flip_bipartition <- function( phy, bi ){

	return( phy$tip.label[ ! phy$tip.label %in% bi ] )

}


#' Get a list of all the bipartitions in a tree.
#'
#' @param phy A phylo object that specifies the tree.
#' @return A list of bipartitions for the tree. The order of the list corresponds to the
#' edges in phy$edge. Bipartitions are specified as a vector of the tip labels that make
#' up one half of the bipartition.
#' @export
get_bipartitions <- function( phy ){
	# Takes a tree, returns a list of
	edge_nums = as.list( 1:nrow( phy$edge ) )

	return( lapply( edge_nums, bipartition_for_edge_by_label, phy=phy ) )

}

#' Given two trees phy1 and phy2 with the same topology and tip labels,
#' get a vector that indicates which node numbers in phy2 correspond to
#' the nodes in phy1
#'
#' @param phy1 A phylo object
#' @param phy2 A phylo object
#' @return A numeric vector in the order of nodes in phy1, providing
#' corresponding node numbers from phy2
#' @export
get_corresponding_nodes <- function( phy1, phy2 ){
	stopifnot( setequal(phy1$tip.label, phy2$tip.label) )

	tip_order = match( phy1$tip.label, phy2$tip.label )

	return( c(tip_order) )

}

#' Assesses how much phy deviates from an ultrametric tree
#'
#' @param phy A phylo object
#' @param model The model used for fitting. "discrete" is used by default for speed
#' @param ... Additional chronos arguments
#' @return The sum of absolute changes in branch length required
#' to make an ape::chronos time calibrated tree, normalized by the total branch
#' length of the calibrated tree. The higher the value, the more the
#' tree deviates from the calibrated tree.
#' @export
difference_from_calibrated <- function( phy, model="discrete", ... ){

	calibrated = ape::chronos( phy, model=model, ...)
	scaling = sum(calibrated$edge.length) / sum(phy$edge.length)
	phy$edge.length = phy$edge.length * scaling
	diff = sum(abs(phy$edge.length - calibrated$edge.length))
	return( diff/sum(calibrated$edge.length) )

}

#' Repartitions lengths along edges that descend from root node so
#' that they are equal. Useful after rooting operations that result
#' in branches with 0 length
#'
#' @param phy A phylo object
#' @return A phylo object with modified edge lengths
#' @export
slide_root_edges <- function( phy ){

	original_length = sum(phy$edge.length)
	root_node = phy$edge[,1][ ! phy$edge[,1] %in% phy$edge[,2] ]
	root_node = unique( root_node )

	# Make sure only one root node is identified
	stopifnot( length(root_node) == 1 )

	# Identify edges that descend from root
	from_root = phy$edge[,1] == root_node

	# Make sure there are only two edges that descend from root
	stopifnot( sum(from_root) == 2 )

	total_length = sum(phy$edge.length[from_root])

	new_length = total_length / 2
	phy$edge.length[from_root] = new_length

	# Make sure total branch length is unchanged
	stopifnot( abs(sum(phy$edge.length) - original_length) < 1e-06 )
	return(phy)

}

#' For each node in the tree (including internal nodes and tips)
#' get the shortest distance to a descendant node. Values for tip
#' nodes should be 0.
#'
#' @param phy A phylo object
#' @return A vector with elements corresponding to each node
#' @export
distance_from_tip = function(phy){

	max_node = max(phy$edge)

	# Make sure nodes are consecutively numbered
	stopifnot( setequal(unique(phy$edge), 1:max_node) )

	# Get the pairwise distances between all nodes
	distance_matrix = ape::dist.nodes(phy)

	ages = sapply(1:max_node, function(x) {
		tips = hutan::tip_descendants(phy, x)
		distances = distance_matrix[x, tips]
		return( min(distances) )
	})

	# Make sure the tips have value 0
	stopifnot( all( ages[1:length(phy$tip.label)]== 0) )

	return(ages)
}

#' For a node in a tree, get a vector of the edges between the
#' node and the root
#'
#' @param phy A phylo object
#' @param node A node number
#' @return A vector of edge numbers
#' @export
ancestral_edges = function(phy, node){
	edges = c()
	focal_node = node

	# While there is an edge with child focal node,
	# add the edge and then sqitch the focal node to its
	#parent.
	while( focal_node %in% phy$edge[ ,2 ] ){
		new_edge = which( phy$edge[ ,2 ] == focal_node )
		edges = c( edges, new_edge )
		focal_node = phy$edge[ new_edge, 1 ]
	}

	return(edges)
}

#' For two nodes in a tree, get a vector of the edges that
#' connect them
#'
#' @param phy A phylo object
#' @param node_a Number of the first node
#' @param node_b Number of the second node
#' @return A vector of edge numbers
#' @export
connecting_edges = function(phy, node_a, node_b){

	# For each node, get the set of edges that connect
	# it to the root
	edges_a = ancestral_edges( phy, node_a )
	edges_b = ancestral_edges( phy, node_b )

	# The path between the nodes is the set of edges that
	# are present in one but not the other, ie xor. This
	# can be obtained by taking the union and then removing
	# the intersection
	edges = setdiff(
			union( edges_a, edges_b ),
			intersect( edges_a, edges_b )
		)

	return(edges)
}

#' Extends each terminal branch by specified length
#'
#' @param phy A phylogeny in ape::phylo format
#' @param x Amount to extend each branch by
#' @return A phylogeny in ape::phylo format
#' @export
extend_terminal_branches = function( phy, x ){

	is_terminal = phy$edge[ , 2 ] %in% 1:length( phy$tip.label )
	phy$edge.length[ is_terminal ] = phy$edge.length[ is_terminal ] + x
	return( phy )

}

### Extended pic functions

#' For each internal node, calculate the difference in state values
#' at the two child nodes
#'
#' @param phy An ape::phylo object
#' @param states A vector with length equal to nodes with the state at each node
#' @return A vector of values corresponding to each internal node with the value difference
#' @importFrom magrittr %>%
#' @export
calc_diffs = function( phy, states ){
	# Get an ordered vector of internal node numbers
	internal_nodes = phy$edge[ , 1 ] %>% unique() %>% sort()

	# Create a matrix with one column per internal node, and
	# two rows - each with the number of a child node
	child_nodes = sapply( internal_nodes, function(n){
		children = phy$edge[ n == phy$edge[,1], 2 ]

		# make sure each internal node has only two children, ie tree is bifurcating
		if( length(children) != 2 ){
			stop("ERROR: Tree is not bifurcating.")
		}
		return( children )
	} )

	# Create a matrix with one column per internal node, and
	# two rows - each with the character state of a child node
	child_states = states[child_nodes]
	dim(child_states) = c(2, length(child_states)/2)
	colnames(child_states) = internal_nodes

	# For each internal node, get the difference between
	# states a descendent nodes
	child_diff = child_states[1, ] - child_states[2, ]
	return(child_diff)
}

#' For each internal node, simulate the difference between state
#' values at the node's children
#'
#' @param phy An ape::phylo object
#' @param model_parameters Parameter estimates for evolutionary model
#' @param model_method The model of trait evolution. Can be one of c("BM","OU")
#' @return The simulated child differences for each internal node
sim_diffs = function( phy, model_parameters, model_method="BM" ){

	states = NA

	if( model_method == "BM"){
		# expects that model_parameters generated by geiger::fitContinuous()

		states = phytools::fastBM(
			phy,
			a = model_parameters$opt$z0,
			sig2 = model_parameters$opt$sigsq,
			bounds = c(-Inf,Inf),
			internal = TRUE
		)

	} else if ( model_method == "OU"){
		# expects that model_parameters generated by geiger::fitContinuous()
		# http://blog.phytools.org/2013/11/new-ou-simulator-in-fastbm.html

		states = phytools::fastBM(
			phy,
			a = model_parameters$opt$z0,
			sig2 = model_parameters$opt$sigsq,
			alpha = model_parameters$opt$alpha,
			theta = model_parameters$opt$z0, # This is done by fastBM anyway, but making it explicit supresses a warning
			bounds = c(-Inf,Inf),
			internal = TRUE
		)

	}


	child_diff = calc_diffs( phy, states )

	return(child_diff)
}

#' Estimate the extended phylogenetic independent contrast. Rather
#' than normalize differences across nodes by branch lengths, differences
#' are normalized be the expected difference obtained from replicate
#' simulations. This allows for greater model flexibility.
#'
#' @param x A numeric vector with one trait value per tip
#' @param phy An ape::phylo object
#' @param var.contrasts logical, indicates whether the expected variances of the contrasts should be returned
#' @param model_method The model of trait evolution. Can be one of c("BM", "OU")
#' @param model_parameters Model parameters from fitContinuous. Will estimate if not provided.
#' @param n_replicates The number of simulations used to estimate the expected differences
#' @return A vector of phylogenetically independent contrasts
#' @importFrom magrittr %>%
#' @export
picx = function( x, phy, var.contrasts=FALSE, model_method="BM", model_parameters=NA, n_replicates = 200 ){

	names( x ) = phy$tip.label

	if( is.na(model_parameters) ){
		model_parameters = geiger::fitContinuous( phy, x, model=model_method, ncores=1 )
	}

	expected_diff =
		replicate( n_replicates, sim_diffs( phy, model_parameters, model_method ) ) %>%
		abs() %>%
		rowMeans( na.rm=TRUE ) # Estimates sometimes have NA

	# Replace 0s with the smallest possible float
	expected_diff[ expected_diff==0 ] = .Machine$double.eps

	if ( model_method=="BM" ){
		ancestral_states = ape::reconstruct( x, phy, method="ML" )$ace
	} else if ( model_method=="OU" ) {
		ancestral_states = ape::reconstruct( x, phy, method="GLS_OUS" )$ace
	}

	states = c( x, ancestral_states )

	actual_diffs = calc_diffs( phy, states )
	contrasts = actual_diffs / expected_diff
	names( contrasts ) = phy$node.label

	if ( var.contrasts ){
		contrasts = cbind( contrasts, expected_diff )
		rownames( contrasts ) = phy$node.label
		colnames( contrasts ) = c( "contrasts", "variance" )
	}

	return( contrasts )
}
