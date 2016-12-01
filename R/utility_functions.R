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

#' Check if two bipartitions drawn from trees with the same tips are compatible with eachother. 
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
#' to make an ape::chronos time calbrated tree, normalized by the total branch 
#' length of the calibrated tree. The higher the value, the more the 
#' tree deviates from the calibrated tree.
difference_from_calibrated <- function( phy, model="discrete", ... ){
	
	calibrated = ape::chronos( phy, model=model, ...)
	scaling = sum(calibrated$edge.length) / sum(phy$edge.length)
	phy$edge.length = phy$edge.length * scaling
	diff = sum(abs(phy$edge.length - calibrated$edge.length))
	return( diff/sum(calibrated$edge.length) )

}
