# hutan: An R package for manipulation of phylogenetic trees

There are many excellent R packages for phylogenetic analyses 
(http://cran.r-project.org/web/views/Phylogenetics.html). Tools for 
topological analysis and manipulation of trees are, however, 
relatively limitted. Such methods are required to preprocess data, 
implement analysis methods, and summarize phylogenetic results.
This package provides a variety of functions for topological 
manipulations and assessments of phylogenetic trees, building on 
those already present in the 
[ape](http://cran.r-project.org/web/packages/ape/index.html) p
ackage.


## Documentation

In R, run `help(hutan)`.

## Example usage

### Create a zero-constrained tree

The following creates a zero-constrained tree as described by Susko 2014 
(http://dx.doi.org/10.1093/molbev/msu039):

	library(hutan)
	data( siphonophore_ml )
	data( siphonophore_constraint )

	zc <- zero_constrained( siphonophore_ml, siphonophore_constraint )
	plot( zc )

## Citing

R has been developed in support of multiple projects, but doesn't yet 
have a dedicated manuscript that describes it. For now, please cite 
[this preprint](http://dx.doi.org/10.1101/005264).

## Installing

In R, run the following:

    library(devtools)
    install_bitbucket('caseywdunn/hutan')

## Development

This package is built with the excellent devtools 
(https://github.com/hadley/devtools).

Development typically involves `cd`ing to the package directory, launching R, 
and running some combination of the following: 
	
	options(error=traceback) # Get line numbers for errors
    library(devtools)
    load_all()
    document()
    test()
    
