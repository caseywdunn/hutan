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

See the manual, 
[hutan-manual.pdf](https://bitbucket.org/caseywdunn/hutan/src/master/hutan-manual.pdf)

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

To find out how to cite hutan, run the following in R:

    citation("hutan")

## Installing

### From the git repository

First, install the [devtools](https://github.com/hadley/devtools) package. Then, run the following in R:

    library(devtools)
    install_bitbucket('caseywdunn/hutan')

## Development

This package is built with the excellent devtools 
(https://github.com/hadley/devtools). Extensive explanations on using devtools 
to create R packages is available in the book 
[R Packages](http://r-pkgs.had.co.nz/).

Development typically involves `cd`ing to the package directory, launching R, 
and running some combination of the following: 
	
	options(error=traceback) # Get line numbers for errors
    library(devtools)
    load_all()
    test()
    document()
    check() # A wrapper for R CMD check, see http://r-pkgs.had.co.nz/check.html#check
    buid() # Create package bundle, including executed vignettes

To regenerate the pdf manual, run the following shell command in the package directory:

    R CMD Rd2pdf . --force --output=hutan-manual.pdf
    
