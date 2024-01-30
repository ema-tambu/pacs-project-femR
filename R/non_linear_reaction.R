# param * (1 +  F) param costante , F Function
.Binomial <- R6Class("Binomial",
                     private = list(
                      param_ = 1.0, 
                      Function_ = "Function"
                     ),
                     public = list(
                       initialize = function(param, Function){
                         private$param_ <- param
                         private$Function_ <- Function
                       },
                       Function = function(){
                         private$Function_
                       },
                       param = function(){
                         private$param_
                       }
                     )
)

#' @export
`-.Function` <- function(e1, e2){
  if(is(e1, "numeric") & is(e2, "Function"))
    .Binomial$new(e1, e2)
  else if(is(e2, "numeric") & is(e1, "Function"))
    .Binomial$new(e2, e1)
}

.NonLinearReaction <- R6Class("NonLinearReaction",
                              inherit = .DifferentialOpCtr)
#' @export
`*.Binomial` <- function(e1, e2){
  if(is(e1, "ReactionOperator") & is(e2, "Binomial") & 
     identical(e1$Function(), e2$Function())){
    .NonLinearReaction$new(params=e1$params[[1]]*e2$param(), # (alpha*f) * [alpha2*(1-f)]
                           tokens = "non_linear_reaction",
                           Function=e1$Function()) # alpha * (1-f) * f
  }else if(is(e1, "Binomial") & is(e2, "Function") & 
           identical(e1$Function(), e2$Function())){
    .NonLinearReaction$new(params=e1$params[[1]]*e2$param(), # alpha*(1-f)*f)
                           tokens = "non_linear_reaction",
                           Function=e1$Function()) # alpha * (1-f) * f
  }
}
