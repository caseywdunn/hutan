# hutan is an R package with a collection of tools for phylogenetic tree manipulation

Run `dev_help('hutan')` or `help(hutan)` for more information.


## Installing

    library(devtools)
    install_bitbucket('caseywdunn/hutan')

## Documentation


## Example usage

### Create a zero-constrained tree

The following creates a zero-constrained tree as described by Susko 2014 (http://dx.doi.org/10.1093/molbev/msu039)

	library(hutan)
	data( siphonophore_ml )
	data( siphonophore_constraint )

	zc <- zero_constrained( siphonophore_ml, siphonophore_constraint )
	plot( zc )

## Citing



## Development

This package is built with devtools: 
https://github.com/hadley/devtools

More on package development with devtools can be found at the 
[devtools site](https://github.com/hadley/devtools/wiki/development).

Development typically involves `cd`ing to the package directory, launching R, 
and running some combination of the following: 
	
	options(error=traceback) # Get line numbers for errors
    library(devtools)
    load_all()
    document()
    test()
    dev_help('hutan')
    
