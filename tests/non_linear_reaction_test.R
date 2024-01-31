devtools::document()
#library(femR)

## load domain data and generate mesh object
data("unit_square", package="femR")

mesh = Mesh(unit_square)
class(mesh)
plot(mesh)

# create Functional Space
fe_order = 1
Vh <- FunctionSpace(mesh, fe_order)

exact_solution <- function(points){
  return( 3*points[,1]^2 + 2*points[,2]^2 )
}

## define differential operator in its strong formulation
u <- Function(Vh)

mu <- 1.
alpha <- 1.
DifferentialOp <- -mu*laplace(u) + alpha*u*(1-u) # problema alpha * (1-u) * u

## forcing term
forcing <- function(points){
  return(-9*points[,1]^4 - 12 * points[,1]^2*points[,2]^2 + 
          3*points[,1]^2 - 4*points[,2]^4 + 2*points[,2]^2 - 10*mu) 
}

## Dirichlet BC
g <- exact_solution

## Pde constructor
pde <- Pde(DifferentialOp = DifferentialOp, 
           forcing = forcing,
           boundary_condition = g)

## solve problem
pde$solve()

## compute L2 norm of the error
u_ex <- as.matrix(exact_solution(pde$dofs_coordinates()))
error.L2 <- sqrt(sum(pde$mass() %*% (u_ex - u$coefficients())^2))
cat("L2 error = ", error.L2, "\n")

u_ex <- Function(Vh)
u_ex$set_coefficients(as.matrix(exact_solution(Vh$basis()$dofs_coordinates())))                      

contour(u_ex)
contour(u)
