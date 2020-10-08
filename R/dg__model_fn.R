#' the dg model function to be minimized
#'
#' @param par vector of starting parameters for modeling
#' @param parlist list with all model parameters
#' @param varlist list with all variables
#'
#' @import data.table
#' @import Matrix
#' @export
#'
dg__model_fn <- function(
	par,
	parlist,
	varlist
) {

    ## ensure fixed parameters are fixed
  	for (i in seq_along(parlist[["fixed_par"]])) {
    	par[names(par) == names(parlist[["fixed_par"]])[i]] <- parlist[["fixed_par"]][i]
  	}


  	## precalc ddg vectors
  	if (parlist[["no_folded_states"]] == 1) {
      f_ddg_var <- list(varlist[["varxmut"]] %*% par[grep("f_ddg", names(par))])
  	} else {
  		f_ddg_var <- list(
        varlist[["varxmut"]] %*% par[grep("fA_ddg", names(par))],
        varlist[["varxmut"]] %*% par[grep("fB_ddg", names(par))]
      )
  	}
  	b_ddg_var <- varlist[["varxmut"]] %*% par[grep("b_ddg", names(par))]


    ## extract global pars to speed up lookups
    global_par <- par[!grepl("ddg", names(par))]


  	## initalize msd
  	msd <- 0


  	## folding phenotype
  	if (varlist[["no_abd_datasets"]] > 0) {
  		for (i in 1:varlist[["no_abd_datasets"]]) {
        f_fitpred <- convert_dg2foldingfitness(
          f_ddg = f_ddg_var,
          f_dgwt = global_par[grep(paste0("^f[AB]?", i, "_dgwt"), names(global_par))],
          f_fitwt = global_par[grep(paste0("f", i, "_fitwt"), names(global_par))],
          f_fit0 = global_par[grep(paste0("f", i, "_fit0"), names(global_par))],
          fitness_scale = varlist[["fitness_scale"]],
          no_folded_states = parlist[["no_folded_states"]]
        )

  			msd <- msd + sum(((varlist[["variant_data"]][, unlist(.SD), .SDcols = paste0("f", i, "_fitness")] -
  					f_fitpred) /
  					varlist[["variant_data"]][, unlist(.SD), .SDcols = paste0("f", i, "_sigma")])^2,
      			na.rm = T)
  		}
  	}


  	## binding phenotype
    if (varlist[["no_bind_datasets"]] > 0) {
      for (i in 1:varlist[["no_bind_datasets"]]) {
        b_fitpred <- convert_dg2bindingfitness(
          b_ddg = b_ddg_var,
          f_ddg = f_ddg_var,
          b_dgwt = global_par[grep(paste0("b", i, "_dgwt"), names(global_par))],
          f_dgwt = global_par[grep(paste0("^bf[AB]?", i, "_dgwt"), names(global_par))],
          b_fitwt = global_par[grep(paste0("b", i, "_fitwt"), names(global_par))],
          b_fit0 = global_par[grep(paste0("b", i, "_fit0"), names(global_par))],
          fitness_scale = varlist[["fitness_scale"]],
          no_folded_states = parlist[["no_folded_states"]]
        )

        msd <- msd + sum(((varlist[["variant_data"]][, unlist(.SD), .SDcols = paste0("b", i, "_fitness")] -
            b_fitpred) /
            varlist[["variant_data"]][, unlist(.SD), .SDcols = paste0("b", i, "_sigma")])^2,
            na.rm = T)
      }
    }


  	## return mean square deviation
  	return(msd)
}