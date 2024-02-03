if(!require("pacman", quietly = T)) install.packages("pacman") 

pacman::p_load("fdaPDE", "ReacTran", "devtools")

if(!require("femR", quietly = T)){
   devtools::install_github(repo="ema-tambu/pacs-project-femR", ref="develop")
}
library(femR)
rmse <- function(x,y){return(sqrt(mean( (x-y)^2)))}

## forcing term
f <- function(points){
  return(-9*points[,1]^4 - 12 * points[,1]^2*points[,2]^2 + 
           3*points[,1]^2 - 4*points[,2]^4 + 2*points[,2]^2 - 10*mu) 
}

#
exact <- function(points){
  return( 3*points[,1]^2 + 2*points[,2]^2 )
}

mu <- 1.
alpha <- 1.
#Dy    <- Dx <- -1.   # diffusion coeff, X- and Y-direction

N = c(16, 32, 64)
errors.l2 <- list("deSolve" = rep(0, times = length(N)),
                  "femR_1" = rep(0, times = length(N)),
                  "femR_2" = rep(0, times = length(N)))

times <- list("deSolve" = rep(0, times = length(N)),
              "femR_1" = rep(0, times = length(N)),
              "femR_2" = rep(0, times = length(N)))
h <- rep(0, times = length(N))
nnodes <- rep(0, times = length(N))

# for(i in 1:length(N))
i = 3
{
  cat("------------ ", N[i], "x",N[i] , " ------------\n")
  x.grid    <- setup.grid.1D(x.up = 0, x.down = 1, N = N[i])
  y.grid    <- setup.grid.1D(x.up = 0, x.down = 1, N = N[i])
  grid2D    <- setup.grid.2D(x.grid, y.grid)
  
  h[i] = max(grid2D$dx, grid2D$dy)
  cat("nodes (deSolve)", N[i]^2,"\n")
  D.grid    <- setup.prop.2D(value = mu, y.value = mu, grid = grid2D)
  v.grid    <- setup.prop.2D(value = 0, grid = grid2D) 
  A.grid    <- setup.prop.2D(value = 1, grid = grid2D)
  VF.grid   <- setup.prop.2D(value = 1, grid = grid2D)
  
  C.x.up    <- exact(cbind(rep(x.grid$x.up, times=N[i]), y.grid$x.mid))# rep(1., times=length(y.grid$x.mid))#exact(rep(x.grid$x.up, times=N[i]), y.grid$x.mid)
  C.y.up    <- exact(cbind(x.grid$x.mid, rep(y.grid$x.up, times=N[i])))#rep(0., times=length(x.grid$x.mid))#exact(x.grid$x.mid, rep(y.grid$x.up, times=N[i]))
  
  C.x.down  <- exact(cbind(rep(x.grid$x.down, times=N[i]), y.grid$x.mid)) #rep(0., times=length(y.grid$x.mid))#exact(rep(x.grid$x.down, times=N[i]), y.grid$x.mid)
  C.y.down  <- exact(cbind(x.grid$x.mid, rep(y.grid$x.down, times=N[i]))) #rep(0., times=length(x.grid$x.mid))#exact(x.grid$x.mid, rep(y.grid$x.down, times=N[i]))
  
  forcing = matrix(nrow=N[i], ncol=N[i], data=0)
  u.ex = matrix(nrow=N[i], ncol=N[i], data=0)
  for(k in 1:N[i]){
    for(l in 1:N[i]){
      forcing[k,l] = f(cbind(grid2D$x.mid[k], grid2D$y.mid[l]))
      u.ex[k,l] = exact(cbind(grid2D$x.mid[k], grid2D$y.mid[l]))
    }
  }
  
  Diff2Db <- function (t, y, parms)  {
    
    Y  <- matrix(nrow = N[i], ncol = N[i], data = y)
    
    dY <- tran.2D(C = Y, 
                  C.x.up = C.x.up, C.x.down = C.x.down,
                  C.y.up = C.y.up, C.y.down = C.y.down,
                  grid = grid2D, 
                  D.grid = D.grid, 
                  A.grid = A.grid, 
                  VF.grid = VF.grid, 
                  v.grid = v.grid)$dC
    
    dY <- dY - alpha*Y*(1-Y) + forcing 
    
    return (list(dY))
  }
  
  y = matrix(data = rep(1,times=N[i]*N[i]))
  y = matrix(data = rnorm(N[i]*N[i], mean = 1, sd =0.5))
  
  start_ = Sys.time()
  Y <- steady.2D(y=y,
                 dimens = c(N[i],N[i]), 
                 time = 0, 
                 func = Diff2Db, parms=NULL, lrw=1e8)
  times$deSolve[i] = as.numeric(difftime(Sys.time(), start_, units="secs"))                 
  cat("deSolve ", times$deSolve[i], " s\n")
  
  u.deSolve <- matrix(Y$y, nrow=N[i], ncol=N[i])
  errors.l2$deSolve[i] = rmse(Y$y, as.vector(u.ex))
  error.deSolve = abs(u.deSolve-u.ex)

  # femR - ORDER 1 ---------------------
  # mesh <- build_mesh(unit_square, 
  #                            maximum_area = h[i]^2, 
  #                            minimum_angle=20)
  grid2D <- expand.grid(x.grid$x.int, y.grid$x.int)
  mesh <- fdaPDE::create.mesh.2D(grid2D)
  nnodes[i] <- nrow(mesh$nodes)
  square <- list(nodes= mesh$nodes, elements= mesh$triangles, boundary= mesh$nodesmarkers)
  
  mesh <- Mesh(square)
  cat("nodes", nrow(mesh$nodes()),"\n")
  Vh <- FunctionSpace(mesh, 1L)
  u <- Function(Vh)
  DifferentialOp <- -mu*laplace(u) + alpha*u*(1-u)
  
  pde <- Pde(DifferentialOp = DifferentialOp, 
             forcing = f,
             boundary_condition = exact)
  
  ## solve problem
  tic_ <- Sys.time() 
  pde$solve()
  times$femR_1[i] <- as.numeric(difftime(Sys.time(), tic_, units="secs"))
  cat("femR (order 1) ", times$femR_1[i], " s\n")
  
  u.ex <- as.matrix(exact(pde$dofs_coordinates()))
  errors.l2$femR_1[i] = rmse(u$coefficients(), u.ex)

}

# contour(u)

u.femR <- matrix(nrow = N[i], ncol= N[i])
U <- u.femR

x.grid    <- setup.grid.1D(x.up = 0, x.down = 1, N = N[i])
y.grid    <- setup.grid.1D(x.up = 0, x.down = 1, N = N[i])
grid2D    <- setup.grid.2D(x.grid, y.grid)
for(k in 1:N[i]){
  for(l in 1:N[i]){
    u.femR[k,l] = u$eval(cbind(grid2D$x.mid[k], grid2D$y.mid[l]))
    U[k,l] = exact(cbind(grid2D$x.mid[k], grid2D$y.mid[l]))
  }
}

filled.contour(u.femR, main="femR")
filled.contour(U, main="exact")
filled.contour(u.deSolve, main="deSolve")
