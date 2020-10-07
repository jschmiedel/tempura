#' calls the optimzation routine to estimate dG models
#'
#' @param start_par vector of starting parameters for modeling
#' @param parlist list with all model parameters
#' @param varlist list with all variables
#' @param maxit maximum number of iterations for optim
#'
#' @return returns a model in data.table format
#' @import data.table
#' @import Matrix
#' @export
#'
dg_model_optim <- function(
	start_par,
	parlist,
	varlist,
  maxit
){

	## call optimizer
	dg_model <- stats::optim(
      par = start_par,
      fn = dg_model_fn,
      gr = dg_model_gr,
      method = "L-BFGS-B",
      lower = parlist[["lower_bounds"]],
      upper = parlist[["upper_bounds"]],
      control = list(maxit = maxit),
      parlist = parlist,
      varlist = varlist
    )

    ## return model
    return(dg_model)
}
