# Load libraries
library(ape)

# An example of how to produce a zero-constraint tree based on the Agalmatidae 
# sensu stricto + Bargmannia constraint of Dunn et al. 2005 
# (http://dx.doi.org/10.1080/10635150500354837)

data( siphonophore_ml )
data( siphonophore_constraint )

zc <- zero_constrained( siphonophore_ml, siphonophore_constraint )
plot( zc )