library(ape)

#' Generates the "zero-constrained" tree described by Susko 2014 
#' (http://dx.doi.org/10.1093/molbev/msu039)
#' 
#' @param phy_resolved A fully resolved phylogeny stored as a phylo object, e.g. an ML 
#' tree.
#' @param phy_constraint A partially resolved constraint tree.
#' @param epsilon 
#' @return A phylo object containing a tree that is the same as phy_resolved, except that 
#' the length of edges that are incompatible with phy_constraint are replaced with epsilon.
#' @examples
#' zc <- zero_constrained( ml, constraint )
zero_constrained <- function ( phy_resolved, phy_constraint, epsilon=0.000001 ){
	phy_resolved$edge.length[ incompatible_edges(  phy_resolved, phy_constraint ) ] <- epsilon
	return( phy_resolved )
}

#' Identify the edges in one phylo object that are incompatible with the edges in another 
#' phylo object. Requires the same tip labels for each tree.
#' 
#' @param phy1 The tree under consideration
#' @param phy2 The tree to be compared to
#' @return A boolean vector corresponding to the edges in phy1. Each element is FALSE if 
#' the edge is compatible with phy2, or TRUE if incompatible.
incompatible_edges <- function( phy1, phy2 ){
	
	# First check to see that the two trees have the same tips
	if( setequal(phy1$tip.label, phy2$tip.label) == FALSE ) {
		stop("Trees do not have the same tip names.")
	}

	bi1 = get_bipartitions( phy1 )
	bi2 = get_bipartitions( phy2 )
	
	incompat = lapply( bi1, is_incompatible_with_set, bi_list=bi2, phy=phy1 )
	
	return ( unlist( incompat ) )
}

#' Check if bipartition bi is incompatible with the bipartitions in bi_list. 
#' Each bipartition is defined as a vector of the names of the tips on one side of the
#' bipartition.
#' 
#' @param bi The query bipartition.
#' @param bi_list A list of the bipartitions to be compared against.
#' @param phy A phylo object describing a tree that includes all tips under investigation. 
#' This is used to infer the other half of each bipartition.
#' @return TRUE if bi is incompatible with any bipartition in bi_list, otherwise FALSE.
is_incompatible_with_set <- function( bi, bi_list, phy ) {

	incompatible <- lapply( bi_list, are_bipartitions_incompatible, bi2=bi, phy=phy )
	return ( any( unlist( incompatible ) ) )
}

#' Check if two bipartitions drawn from trees with the same tips are incompatible with eachother. 
#' Each bipartition is defined as a vector of the names of the tips on one side of the
#' bipartition.
#' 
#' @param bi1 The first bipartition.
#' @param bi2 The second bipartition.
#' @param phy A phylo object describing a tree that includes all tips under investigation. 
#' This is used to infer the other half of each bipartition.
#' @return TRUE if bi1 is incompatible with bi2, otherwise FALSE.
are_bipartitions_incompatible <- function( bi1, bi2, phy ) {
	
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
	
	return( ! compatible )
}

#' Get tips and labels of a phylo object.
#' 
#' @param phy A phylo object.
#' @return A vector of all the tips, annotated with their names
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
bipartition_for_edge_by_label <- function( edge, phy ){
	
	return( phy$tip.label[ bipartition_for_edge( phy, edge ) ]	 )

}

#' Get a list of all the bipartitions in a tree.
#' 
#' @param phy A phylo object that specifies the tree.
#' @return A list of bipartitions for the tree. The order of the list corresponds to the 
#' edges in phy$edge. Bipartitions are specified as a vector of the tip labels that make 
#' up one half of the bipartition.
get_bipartitions <- function( phy ){
	# Takes a tree, returns a list of 
	edge_nums = as.list( 1:nrow( phy$edge ) )
	
	return( lapply( edge_nums, bipartition_for_edge_by_label, phy=phy ) )
	
}

#######################################
# Examples

#' Generates example output of zero_constrained based on the Agalmatidae sensu stricto + 
#' Bargmannia constraint of Dunn et al. 2005 (http://dx.doi.org/10.1080/10635150500354837)
#'  
#' @return A vector of three phylo objects: the ml tree, the constraint tree, and the zero-constrained tree
#' @examples
#' trees <- run_example()
#' plot( trees[ 3 ] )
run_example <- function(){
	ml = read.tree(text="(Diphyes_dispar:0.22692492488504159565,(Chelophyes_appendiculata:0.56688200733654858787,((((((((((Cordagalma_cordiforme:0.19153120394744349575,(Vogtia_glabra:0.27178535158124850213,(Vogtia_pentacantha:0.18558670643435176695,(Hippopodius_hippopus_Pacific:0.08157590402101695670,Hippopodius_hippopus_Atlantic:0.01856755234157863796):0.20024032399709384977):0.02183508876652192829):0.30663054909489212418):0.19891026771712405630,(((Forskalia_asymmetrica:0.11490532525079964754,(Forskalia_formosa:0.05160455192736665420,((Forskalia_edwardsi_Atlantic:0.00000224734284353031,(Forskalia_edwardsi_Pacific_1:0.00000224734284353031,Forskalia_edwardsi_Pacific_2:0.00000224734284353031):0.00412525089215428752):0.02248486623902519813,Forskalia_tholoides:0.06167773790809999435):0.03847579686021841000):0.03082645504597005676):0.15497716566979480302,Physophora_hydrostatica:0.12238075433245980406):0.01872790657600384434,((Nanomia_bijuga_Pacific:0.04453911678446929867,Nanomia_bijuga_Atlantic:0.07460365969277675535):0.05107949913978521855,((Halistemma_rubrum_Pacific:0.00366093532817418112,(Halistemma_rubrum_Atlantic:0.00449577602737188942,Halistemma_rubrum_Med:0.00000224734284353031):0.00549201394213728183):0.03450707156613055243,(Agalma_clausi:0.01476660529730974836,((Agalma_elegans_Atlantic:0.00000224734284353031,Agalma_elegans_Pacific:0.00000224734284353031):0.03187418352818011807,(Agalma_okeni:0.02178679449485416975,(Athorybia_rosacea_Atlantic:0.00144774322266859627,Athorybia_rosacea_Pacific:0.00317491068777582975):0.03308186825951565935):0.01104841643220617897):0.00572757494939845225):0.03507214271150678464):0.01606723659496930162):0.06274123574699491668):0.01123163359043143072):0.02400526861166018405,((((Apolemia_sp_3:0.00120035375267486810,(Apolemia_sp_4:0.00288674009202094213,Apolemia_sp_1:0.00833049218938709292):0.01432486320729266065):0.00864728897649716562,Apolemia_sp_2:0.00415089566189395538):0.06059773022817797200,((Staurocladia_wellingtoni:0.15847319991869418532,((Porpita_porpita:0.08496722554500601987,Velella_velella:0.08114608521919654260):0.24701392679291822585,Hydra:0.42168906611701839626):0.01733183950441345711):0.09762306408444636208,(Physalia_physalis:0.20650480923349964768,(Rhizophysa_filiformis:0.08283380891531452739,Rhizophysa_eysenhardti:0.06193625069341131967):0.02538129554939630439):0.08782382463981050491):0.04448643879911705346):0.02255685932353277551,(Bargmannia_amoena:0.03312246762809323286,Bargmannia_elongata:0.03326302226225569258):0.06804410003878090529):0.01576627217692200614):0.00607225751780783932,(Stephalia_dilata:0.05105558432743991654,(Erenna_sp:0.03011526982724788137,Stephanomia_amphytridis:0.02850631655327572414):0.01106247400834093828):0.00640929348350334233):0.02956976282976982173,(Praya_dubia:0.11260157379132591793,(Nectopyramis_natans:0.07865777450988088726,Nectadamas_diomedeae:0.06231942694425601087):0.02856339672510353755):0.00849491257061945476):0.05454390570365009322,(Craseoa_lathetica:0.14875777146793811578,(Gymnopraia_lapislazula:0.30385840514508777321,Rosacea_flaccida:0.25478495719644539408):0.06101911279191868376):0.06798001289189997687):0.05667870530279224550,Clausophyid_sp_1:0.10694144017219688048):0.10133553961269042842,(Chuniphyes_multidentata:0.18503295070735162331,Clausophyes_ovata:0.17371321160888641977):0.08470486431729384869):0.12306942347463049880,Sphaeronectes_gracilis:0.32163279354293178303):0.16697742036412610567,(Muggiaea_atlantica:0.18174272643378777681,(Lensia_conoidea:0.19411997246707121678,(Sulculeolaria_quadrivalvis_Atlantic:0.09853151892235978426,Sulculeolaria_quadrivalvis_Pacific:0.03959874271153360908):0.17712679190138830299):0.09689228082319385760):0.05065274456105777617):0.03687128412829689117):0.03610131213107283660,Abylopsis_tetragona:0.20667527268400712193):0.0;")

	constraint = read.tree(text="((Nanomia_bijuga_Atlantic,Nanomia_bijuga_Pacific,Agalma_clausi,Agalma_elegans_Atlantic,Agalma_elegans_Pacific,Agalma_okeni,Athorybia_rosacea_Atlantic,Athorybia_rosacea_Pacific,Halistemma_rubrum_Atlantic,Halistemma_rubrum_Med,Halistemma_rubrum_Pacific,Bargmannia_amoena,Bargmannia_elongata),Apolemia_sp_1,Apolemia_sp_2,Apolemia_sp_3,Apolemia_sp_4,Abylopsis_tetragona,Chelophyes_appendiculata,Chuniphyes_multidentata,Clausophyes_ovata,Clausophyid_sp_1,Cordagalma_cordiforme,Craseoa_lathetica,Diphyes_dispar,Erenna_sp,Forskalia_asymmetrica,Forskalia_edwardsi_Atlantic,Forskalia_edwardsi_Pacific_1,Forskalia_edwardsi_Pacific_2,Forskalia_formosa,Forskalia_tholoides,Gymnopraia_lapislazula,Hippopodius_hippopus_Atlantic,Hippopodius_hippopus_Pacific,Hydra,Lensia_conoidea,Muggiaea_atlantica,Nectadamas_diomedeae,Nectopyramis_natans,Physalia_physalis,Physophora_hydrostatica,Porpita_porpita,Praya_dubia,Rhizophysa_eysenhardti,Rhizophysa_filiformis,Rosacea_flaccida,Sphaeronectes_gracilis,Staurocladia_wellingtoni,Stephalia_dilata,Stephanomia_amphytridis,Sulculeolaria_quadrivalvis_Atlantic,Sulculeolaria_quadrivalvis_Pacific,Velella_velella,Vogtia_glabra,Vogtia_pentacantha);")

	zc <- zero_constrained( ml, constraint )
	return( c( ml, constraint, zc ) )

}