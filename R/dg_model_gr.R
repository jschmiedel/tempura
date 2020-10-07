#' gradient for the dg model function to be minimized
#'
#' @param par vector of starting parameters for modeling
#' @param parlist list with all model parameters
#' @param varlist list with all variables
#'
#' @import data.table
#' @import Matrix
#' @export
#'
dg_model_gr <- function(
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


  	## folding phenotype
    gradient_f = list()
    for (i in 1:varlist[["no_abd_datasets"]]) {
        gradient_f[[i]] <- convert_dg2foldinggradient(
            f_ddg = f_ddg_var,
            f_dgwt = global_par[grep(paste0("^f[AB]?", i, "_dgwt"), names(global_par))],
            f_fitwt = global_par[grep(paste0("f", i, "_fitwt"), names(global_par))],
            f_fit0 = global_par[grep(paste0("f", i, "_fit0"), names(global_par))],
            fitness = varlist[["variant_data"]][, unlist(.SD), .SDcols = paste0("f", i, "_fitness")],
            w = varlist[["variant_data"]][, unlist(.SD), .SDcols = paste0("f", i, "_sigma")],
            mutxvar = varlist[["mutxvar"]],
            fitness_scale = varlist[["fitness_scale"]],
            no_folded_states = parlist[["no_folded_states"]]
        )
    }


    ## binding phenotype
    gradient_b = list()
    for (i in 1:varlist[["no_bind_datasets"]]) {
        gradient_b[[i]] <- convert_dg2bindinggradient(
            b_ddg = b_ddg_var,
            f_ddg = f_ddg_var,
            b_dgwt = global_par[grep(paste0("b", i, "_dgwt"), names(global_par))],
            f_dgwt = global_par[grep(paste0("^bf[AB]?", i, "_dgwt"), names(global_par))],
            b_fitwt = global_par[grep(paste0("b", i, "_fitwt"), names(global_par))],
            b_fit0 = global_par[grep(paste0("b", i, "_fit0"), names(global_par))],
            fitness = varlist[["variant_data"]][, unlist(.SD), .SDcols = paste0("b", i, "_fitness")],
            w = varlist[["variant_data"]][, unlist(.SD), .SDcols = paste0("b", i, "_sigma")],
            mutxvar = varlist[["mutxvar"]],
            fitness_scale = varlist[["fitness_scale"]],
            no_folded_states = parlist[["no_folded_states"]]
        )
    }


    ## list and sum gradients from individual contributions
    f_idx = 1:varlist[["no_abd_datasets"]]
    b_idx = 1:varlist[["no_bind_datasets"]]
    if (parlist[["no_folded_states"]] == 1) {
        gr_f_ddg <- rowSums(sapply(f_idx, FUN = function(f_idx){as.numeric(gradient_f[[f_idx]][["f_ddg"]])})) +
                rowSums(sapply(b_idx, FUN = function(b_idx){as.numeric(gradient_b[[b_idx]][["f_ddg"]])}))

        gr_f_dgwt <- c()
        for (i in f_idx) {
          gr_f_dgwt[i] <- gradient_f[[i]][["f_dgwt"]]
        }

        gr_bf_dgwt <- c()
        for (i in b_idx) {
          gr_bf_dgwt[i] <- gradient_b[[i]][["f_dgwt"]]
        }
    } else {
        gr_fA_ddg <- rowSums(sapply(f_idx, FUN = function(f_idx){as.numeric(gradient_f[[f_idx]][["fA_ddg"]])})) +
                rowSums(sapply(b_idx, FUN = function(b_idx){as.numeric(gradient_b[[b_idx]][["fA_ddg"]])}))
        gr_fB_ddg <- rowSums(sapply(f_idx, FUN = function(f_idx){as.numeric(gradient_f[[f_idx]][["fB_ddg"]])})) +
                rowSums(sapply(b_idx, FUN = function(b_idx){as.numeric(gradient_b[[b_idx]][["fB_ddg"]])}))
        gr_f_ddg <- c(gr_fA_ddg, gr_fB_ddg)

        gr_fA_dgwt <- c()
        gr_fB_dgwt <- c()
        for (i in f_idx) {
          gr_fA_dgwt[i] <- gradient_f[[i]][["fA_dgwt"]]
          gr_fB_dgwt[i] <- gradient_f[[i]][["fB_dgwt"]]
        }
        gr_f_dgwt <- c(gr_fA_dgwt, gr_fB_dgwt)

        gr_bfA_dgwt <- c()
        gr_bfB_dgwt <- c()
        for (i in b_idx) {
          gr_bfA_dgwt[i] <- gradient_b[[i]][["fA_dgwt"]]
          gr_bfB_dgwt[i] <- gradient_b[[i]][["fB_dgwt"]]
        }
        gr_bf_dgwt <- c(gr_bfA_dgwt, gr_bfB_dgwt)
    }

    gr_b_ddg <- rowSums(sapply(b_idx, FUN = function(b_idx){as.numeric(gradient_b[[b_idx]][["b_ddg"]])}))

    gr_f_fitwt <- c()
    gr_f_fit0 <- c()
    for (i in f_idx) {
      gr_f_fitwt[i] <- gradient_f[[i]][["f_fitwt"]]
      gr_f_fit0[i] <- gradient_f[[i]][["f_fit0"]]
    }

    gr_b_dgwt <- c()
    gr_b_fitwt <- c()
    gr_b_fit0 <- c()
    for (i in b_idx) {
      gr_b_dgwt[i] <- gradient_b[[i]][["b_dgwt"]]
      gr_b_fitwt[i] <- gradient_b[[i]][["b_fitwt"]]
      gr_b_fit0[i] <- gradient_b[[i]][["b_fit0"]]
    }

    # collect gradients
    gradient <- c(
      gr_f_ddg,
      gr_b_ddg,
      gr_f_dgwt, gr_f_fitwt, gr_f_fit0,
      gr_b_dgwt, gr_bf_dgwt, gr_b_fitwt, gr_b_fit0
    )
    names(gradient) <- names(par)


    ## set gradients from fixed par to 0
    for (i in seq_along(parlist[["fixed_par"]])) {
        gradient[names(gradient) == names(parlist[["fixed_par"]])[i]] <- 0
    }


    ## return gradient
    return(gradient)

}
