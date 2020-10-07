#' Prepare parameters for a specific dG model
#'
#'
#' @param dataset_folder absolute path to the dataset folder, is created if non-existent
#' @param model_name name of the model that should be computed on the dataset
#' @param no_folded_states integer, '1': one folded, binding-competent state, '2': binding-incompetent and binding-competent folded states, default = 1
#' @param dg_extremes float, lower and upper bounds for dG parameters, default = 15
#' @param fixed_par list of parameters that are being fixed at their starting values for model fitting
#'
#' @return writes a .RData file to $base_folder/#dataset_name/#model_name/parameter_list.RData containing the parlist list with all necessary parameters to compute a specific dG model
#' @import data.table
#' @export
prepare_dg_model <- function(
	dataset_folder,
	model_name,
	no_folded_states = 1,
	dg_extremes = 15,
	fixed_par = c()
) {

	## load varlist
    load(file = file.path(dataset_folder, "data/fitness_dataset.RData"))

	## create dataset folders, if not already present
	dir.create(file.path(dataset_folder))
	dir.create(file.path(dataset_folder, model_name))


	## define parameter names for dG modeling
	if (no_folded_states == 1) {
		f_mut_ddg <- paste0(varlist[["mut_list"]][, mutation], "_f_ddg")
		f_global <- c(paste0("f", 1:varlist[["no_abd_datasets"]], "_dgwt"),
				paste0("f", 1:varlist[["no_abd_datasets"]], "_fitwt"),
				paste0("f", 1:varlist[["no_abd_datasets"]], "_fit0"))
		b_global <- c(paste0("b", 1:varlist[["no_bind_datasets"]], "_dgwt"),
				paste0("bf", 1:varlist[["no_bind_datasets"]], "_dgwt"),
				paste0("b", 1:varlist[["no_bind_datasets"]], "_fitwt"),
				paste0("b", 1:varlist[["no_bind_datasets"]], "_fit0"))
	} else { #label two folded states as A and B, where B is the binding-competent state
		f_mut_ddg <- c(paste0(varlist[["mut_list"]][, mutation], "_fA_ddg"),
				paste0(varlist[["mut_list"]][, mutation], "_fB_ddg"))
		f_global <- c(paste0("fA", 1:varlist[["no_abd_datasets"]], "_dgwt"),
				paste0("fB", 1:varlist[["no_abd_datasets"]], "_dgwt"),
				paste0("f", 1:varlist[["no_abd_datasets"]], "_fitwt"),
				paste0("f", 1:varlist[["no_abd_datasets"]], "_fit0"))
		b_global <- c(paste0("b", 1:varlist[["no_bind_datasets"]], "_dgwt"),
				paste0("bfA", 1:varlist[["no_bind_datasets"]], "_dgwt"),
				paste0("bfB", 1:varlist[["no_bind_datasets"]], "_dgwt"),
				paste0("b", 1:varlist[["no_bind_datasets"]], "_fitwt"),
				paste0("b", 1:varlist[["no_bind_datasets"]], "_fit0"))
	}
	b_mut_ddg <- paste0(varlist[["mut_list"]][, mutation], "_b_ddg")

	par_names <- c(f_mut_ddg, b_mut_ddg, f_global, b_global)


    ## bounds for solver
    lower_bounds <- rep(-dg_extremes, length(par_names))
    upper_bounds <- rep(dg_extremes, length(par_names))
    if (varlist[["fitness_scale"]] == "lin") {
    	lower_bounds[grep("fitwt", par_names)] <- 0.9 #assume wild-type fitness within 10% of expected value
    	lower_bounds[grep("fit0", par_names)] <- 0

    	upper_bounds[grep("fitwt", par_names)] <- 1.1 #assume wild-type fitness within 10% of expected value
    	upper_bounds[grep("fit0", par_names)] <- 0.9 #lethal fitness always smaller than fitwt
    } else {
    	lower_bounds[grep("fitwt", par_names)] <- -0.1
    	lower_bounds[grep("fit0", par_names)] <- -Inf

    	upper_bounds[grep("fitwt", par_names)] <- 0.1
    	upper_bounds[grep("fit0", par_names)] <- -0.1
    }

    ## mean and sd to draw starting parameters values for solver
    start_par_mean_sd <- matrix(0, nrow = length(par_names), ncol = 2) #all ddg values are not sampled
    # dgwt values are drawn from N(0, 1)
    start_par_mean_sd[grep("dgwt", par_names), ] <- rep(c(0, 1), each = length(grep("dgwt", par_names)))

    # fitwt values are drawn from N(fitwt, 0.05)
    if (varlist[["fitness_scale"]] == "lin") {
    	start_par_mean_sd[grep("fitwt", par_names), ] <- rep(c(1, 0.05), each = length(grep("fitwt", par_names)))
    } else {
    	start_par_mean_sd[grep("fitwt", par_names), ] <- rep(c(0, 0.05), each = length(grep("fitwt", par_names)))
    }

    # fit0 values are drawn from N(mean(lower_half_fitness_distribution), sd(lower_half_fitness_distribution))
    for (i in 1:varlist[["no_abd_datasets"]]) {
    	q50 <- varlist[["variant_data"]][,stats::quantile(.SD, 0.5, na.rm = T), .SDcols = paste0("f", i, "_fitness")]
    	mean_q50 <- varlist[["variant_data"]][get(paste0("f", i, "_fitness")) < q50, mean(unlist(.SD)), .SDcols = paste0("f", i, "_fitness")]
    	sd_q50 <- varlist[["variant_data"]][get(paste0("f", i, "_fitness")) < q50, stats::sd(unlist(.SD)), .SDcols = paste0("f", i, "_fitness")]
    	start_par_mean_sd[grep(paste0("f", i, "_fit0"), par_names), ] <- c(mean_q50, sd_q50)
    }
    for (i in 1:varlist[["no_bind_datasets"]]) {
    	q50 <- varlist[["variant_data"]][,stats::quantile(.SD, 0.5, na.rm = T), .SDcols = paste0("b", i, "_fitness")]
    	mean_q50 <- varlist[["variant_data"]][get(paste0("b", i, "_fitness")) < q50, mean(unlist(.SD)), .SDcols = paste0("b", i, "_fitness")]
    	sd_q50 <- varlist[["variant_data"]][get(paste0("b", i, "_fitness")) < q50, stats::sd(unlist(.SD)), .SDcols = paste0("b", i, "_fitness")]
    	start_par_mean_sd[grep(paste0("b", i, "_fit0"), par_names), ] <- c(mean_q50, sd_q50)
    }


    ## collect all parameters
    parlist <- list(par_names = par_names,
				lower_bounds = lower_bounds,
				upper_bounds = upper_bounds,
				start_par_mean_sd = start_par_mean_sd,
				no_folded_states = no_folded_states,
				fixed_par = fixed_par)

	## save varlist to .RData file
	save(parlist, file = file.path(dataset_folder, model_name, "parameter_list.RData"))
}
