#' Prepare parameters for a specific dG model
#'
#'
#' @param dataset_folder absolute path to the dataset folder, is created if non-existent
#' @param model_name name of the model that should be computed on the dataset
#' @param no_folded_states integer, '1': one folded, binding-competent state, '2': binding-incompetent and binding-competent folded states, default = 1
#' @param dg_extremes float, lower and upper bounds for dG parameters, default = 15
#' @param fixed_par list of parameters that are being fixed at their starting values for model fitting
#'
#' @return writes a .RData file to $dataset_folder/$model_name/parameter_list.RData containing the parlist list with all necessary parameters to compute a specific dG model
#' @import data.table
#' @export
dg_prepare_model <- function(
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

	dt_par <- data.table(parameter = c(f_mut_ddg, b_mut_ddg, f_global, b_global))


    ## bounds for solver
    dt_par[, lower_bound := -dg_extremes]
    dt_par[, upper_bound := dg_extremes]
    if (varlist[["fitness_scale"]] == "lin") {
    	dt_par[grep("fitwt", parameter), lower_bound := 0.9] #assume wild-type fitness within 10% of expected value
        dt_par[grep("fit0", parameter), lower_bound := 0]

        dt_par[grep("fitwt", parameter), upper_bound := 1.1] #assume wild-type fitness within 10% of expected value
        dt_par[grep("fit0", parameter), upper_bound := 0.9] #lethal fitness always smaller than fitwt
    } else {
    	dt_par[grep("fitwt", parameter), lower_bound := -0.1]
        dt_par[grep("fit0", parameter), lower_bound := -Inf]

        dt_par[grep("fitwt", parameter), upper_bound := 0.1]
        dt_par[grep("fit0", parameter), upper_bound := -0.1]
    }

    ## mean and sd to draw starting parameters values for solver
    dt_par[, start_par_mean := 0]
    dt_par[, start_par_sd := 0] #all ddg values are not sampled
    # dgwt values are drawn from N(0, 1)
    dt_par[grep("dgwt", parameter), ':=' (start_par_mean = 0, start_par_sd = 1)]

    # fitwt values are drawn from N(fitwt, 0.05)
    if (varlist[["fitness_scale"]] == "lin") {
        dt_par[grep("fitwt", parameter), ':=' (start_par_mean = 1, start_par_sd = 0.05)]
    } else {
        dt_par[grep("fitwt", parameter), ':=' (start_par_mean = 0, start_par_sd = 0.05)]
    }

    # fit0 values are drawn from N(mean(lower_half_fitness_distribution), sd(lower_half_fitness_distribution))
    for (i in 1:varlist[["no_abd_datasets"]]) {
    	q50 <- varlist[["variant_data"]][,stats::quantile(.SD, 0.5, na.rm = T), .SDcols = paste0("f", i, "_fitness")]
    	mean_q50 <- varlist[["variant_data"]][get(paste0("f", i, "_fitness")) < q50, mean(unlist(.SD)), .SDcols = paste0("f", i, "_fitness")]
    	sd_q50 <- varlist[["variant_data"]][get(paste0("f", i, "_fitness")) < q50, stats::sd(unlist(.SD)), .SDcols = paste0("f", i, "_fitness")]
        dt_par[grep(paste0("f", i, "_fit0"), parameter), ':=' (start_par_mean = mean_q50, start_par_sd = sd_q50)]
    }
    for (i in 1:varlist[["no_bind_datasets"]]) {
    	q50 <- varlist[["variant_data"]][,stats::quantile(.SD, 0.5, na.rm = T), .SDcols = paste0("b", i, "_fitness")]
    	mean_q50 <- varlist[["variant_data"]][get(paste0("b", i, "_fitness")) < q50, mean(unlist(.SD)), .SDcols = paste0("b", i, "_fitness")]
    	sd_q50 <- varlist[["variant_data"]][get(paste0("b", i, "_fitness")) < q50, stats::sd(unlist(.SD)), .SDcols = paste0("b", i, "_fitness")]
    	dt_par[grep(paste0("b", i, "_fit0"), parameter), ':=' (start_par_mean = mean_q50, start_par_sd = sd_q50)]
    }


    ## key dt_par in alphabetical order of parameter names
    data.table::setkey(dt_par, parameter)


    ## collect all parameters
    parlist <- list(dt_par = dt_par,
				no_folded_states = no_folded_states,
				fixed_par = fixed_par)

	## save varlist to .RData file
	save(parlist, file = file.path(dataset_folder, model_name, "parameter_list.RData"))
}
