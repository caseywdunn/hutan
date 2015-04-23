# hutan is an R package with a collection of tools for phylogenetic tree manipulation

Run `dev_help('hutan')` or `help(hutan)` for more information.


## Installing

    library(devtools)
    install_bitbucket('caseywdunn')
    install_bitbucket('hutan')

## Development

This package has been built with devtools: 
https://github.com/hadley/devtools

More on package development with devtools can be found at the 
[devtools site](https://github.com/hadley/devtools/wiki/development).

Development typically involves `cd`ing to the package directory, launching R, 
and running some combination of the following: 
	
	options(error=traceback) # Get line numbers for errors
    library(devtools)
    load_all()
    document()
    dev_example('hutan')
    dev_help('hutan')
    
