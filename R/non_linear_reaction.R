.NonLinearReaction <- R6Class("NonLinearReaction",
                              inherit = .DifferentialOpCtr)
setOldClass(c("NonLinearReaction" , "DifferentialOp"))

# param * (1 -  Function) param costante , F Function
# ideale: alpha * (A + B * F) 
# IDEA -> Binomial eredita da DifferentialOp quindi tutte le operazioni sono disponibili ...
# a questo punto nel prodotto per DifferentialOp.. metti un if/else per caso e1==Binomial/e2==Binomial -> NonLinearReaction !
# in questo modo in futuro potrai gestire K(x,u)*u beta(x,u) :)
.Binomial <- R6Class("Binomial", inherit = .DifferentialOpCtr,
                     public = list(
                       initialize = function(params = c(1.,1.,1.), Function){
                         super$initialize(tokens = "binomial", 
                                          params = list(params),
                                          Function = Function)
                       }
                     )
)

setOldClass(c("Binomial", "NonLinearReaction"))

# (A-F) // (F-A) 
#' @export
`-.Function` <- function(e1, e2){
  if(is(e1, "numeric") & is(e2, "Function"))
    .Binomial$new(params = c(1.,1.,-1), Function = e2) 
  else if(is(e2, "numeric") & is(e1, "Function"))
    .Binomial$new(params = c(1.,-1.,1.), Function = e1)
}

# (A+F) // (F+A) 
#' @export
`+.Function` <- function(e1, e2){
  if(is(e1, "numeric") & is(e2, "Function"))
    .Binomial$new(Function = e2) #
  else if(is(e2, "numeric") & is(e1, "Function"))
    .Binomial$new(Function = e1)
}

# casi (c*f)*[alpha*(A-B*F)] 
#      (c*f)*[alpha*(A+B*F)] in operators.R :)