# (c + F) c costante , F Function
.Binomial <- R6Class("Binomial",
                     private = list(
                      c_ = 0.0,
                      Function_ = "Function"
                     ),
                     public = list(
                       initialize = function(c, Function){
                         private$c_ <- c
                         private$Function_ <- Function
                       },
                       Function(){
                         private$Function_
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
