library(femR)

## load domain data and generate mesh object
data("unit_square", package="femR")

mesh = Mesh(unit_square)
class(mesh)
plot(mesh)

# create Functional Space
fe_order = 2
Vh <- FunctionSpace(mesh, fe_order)

exact_solution <- function(points){
  return( sin(2. * pi * points[,1]) * sin(2. * pi * points[,2]) )
}

## define differential operator in its strong formulation
f <- Function(Vh)

alpha <- 0.9
L <- -laplace(f) + alpha*f*(1-f) # problema alpha * (1-f) * f
## forcing term
u <- function(points){
  return(8.*pi^2* sin( 2.* pi * points[,1]) * sin(2.*pi* points[,2]) ) 
}
## Dirichlet BC
g <- function(points){
  return(rep(0, times=nrow(points)))
}