#' calls the optimzation routine to estimate dG models
#'
#' @param start_par vector of starting parameters for modeling
#' @param parlist list with all model parameters
#' @param varlist list with all variables
#' @param maxit maximum number of iterations for optim
#' @param trace_optim logical, if TRUE, shows trace = 3, default: FALSE
#'
#' @return returns a model in data.table format
#' @import data.table
#' @import Matrix
#' @export
#'
dg__model_optim <- function(
	start_par,
	parlist,
	varlist,
  maxit,
  trace_optim = FALSE
){

  lower_bound <- upper_bound <- NULL

  if (trace_optim == TRUE) {
    trace = 3
  } else {
    trace = 0
  }
	## call optimizer
	dg_model <- stats::optim(
      par = start_par,
      fn = dg__model_fn,
      gr = dg__model_gr,
      method = "L-BFGS-B",
      lower = parlist[["dt_par"]][, lower_bound],
      upper = parlist[["dt_par"]][, upper_bound],
      control = list(maxit = maxit, trace = trace),
      parlist = parlist,
      varlist = varlist
    )

    ## return model
    return(dg_model)
}
