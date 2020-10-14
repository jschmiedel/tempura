#' Estimate free energies (dG) for a model from a dataset
#'
#'
#' @param dataset_folder absolute path to the dataset folder, is created if non-existent
#' @param model_name name of the model that should be computed on the dataset
#' @param which_test_set integer, test_set to exclude from model training, default = 0, i.e. will be chosen depending on iteration parameter
#' @param iteration integer, to calucate a set of models and get reproducible results
#' @param return_model logical, if TRUE, returns a data.table with fitted values for all model parameters, default = FALSE
#' @param maxit integer, maximum number of iterations by optim algorithm, default = 1e4
#' @param trace_optim logical, if TRUE, shows trace = 3, default: FALSE
#'
#' @return writes the model parameters as .txt file to $dataset_folder/$model_name/tmp/dg_model_$testset_$iteration
#' @import data.table
#' @import Matrix
#'
#' @export
#'
dg_estimate <- function(
	dataset_folder,
	model_name,
	which_test_set = 0,
	iteration = 1,
    return_model = FALSE,
    maxit = 1e4,
    trace_optim = FALSE
) {

    test_set <- parlist <- start_par_mean <- start_par_sd <- parameter <- lower_bound <- upper_bound <- NULL

    ## load preprocessed files
    # load varlist
    load(file = file.path(dataset_folder, "data/fitness_dataset.RData"))
    # load parlist
    load(file = file.path(dataset_folder, model_name, "parameter_list.RData"))


	## set seed for reproducible model results
    set.seed(iteration)


    ## exclude test_set variants
    # if which_test_set == 0, the script adaptively chooses the test_set and its model iteration
    # this facilitates running tempura in batch-mode from the command line
    if (which_test_set == 0) {
        # get all test_sets
        test_sets <- sort(unique(varlist[["variant_data"]][, test_set]))
        print(paste0("train/test sets in dataset: ", paste0(test_sets, collapse = " ")))
        # assign test_set given the modulus of iteration / # test_sets (+1 to fit test_set enumeration)
        which_test_set <- test_sets[ifelse((iteration %% length(test_sets)) != 0, iteration %% length(test_sets), test_sets[length(test_sets)])]
        # adjust iteration parameter (given iteration parameters runs through repeated list of all test sets)
        iteration <- ceiling(iteration / length(test_sets))
    }

    print(paste0("test set excluded: ", which_test_set))
    print(paste0("iteration: ", iteration))


    ## exclude test set variants
    train_set <- varlist[["variant_data"]][, which(test_set != which_test_set)]
    varlist[["variant_data"]] <- varlist[["variant_data"]][train_set]
    varlist[["varxmut"]] <- varlist[["varxmut"]][train_set, ]
    # transpose varxmut matrix
    varlist[["mutxvar"]] <- Matrix::t(varlist[["varxmut"]])


    ## sample start parameters according to start_par_mean_sd
    start_par <- stats::rnorm(
    	nrow(parlist[["dt_par"]]),
    	mean = parlist[["dt_par"]][, start_par_mean],
    	sd = parlist[["dt_par"]][, start_par_sd]
    )
    names(start_par) <- parlist[["dt_par"]][, parameter]


    ## fix parameters
    for (i in seq_along(parlist[["fixed_par"]])) {
      start_par[names(start_par) == names(parlist[["fixed_par"]])[i]] <- parlist[["fixed_par"]][i]
    }


    ## enforce lower/upper bounds
    start_par[start_par < parlist[["dt_par"]][, lower_bound]] = parlist[["dt_par"]][start_par < lower_bound, lower_bound]
    start_par[start_par > parlist[["dt_par"]][, upper_bound]] = parlist[["dt_par"]][start_par > upper_bound, upper_bound]


    ## call optimzer
    dg_model <- dg__model_optim(
        start_par = start_par,
        parlist = parlist,
        varlist = varlist,
        maxit = maxit,
        trace_optim = trace_optim
    )


    ## convert output to data.table
    dg_model <- data.table::data.table(
        t(dg_model[["par"]]),
        dg_model[["value"]],
        dg_model[["convergence"]],
        which_test_set,
        iteration
    )
    names(dg_model) <- c(
        parlist[["dt_par"]][, parameter],
        "objective",
        "convergence",
        "test_set",
        "iteration"
    )


    ## write model to file
	# create tmp dir
	base::dir.create(file.path(dataset_folder, model_name, "tmp"))
	# write model
	utils::write.table(dg_model,
		file = file.path(dataset_folder, model_name, paste0("tmp/dg_model_testset", which_test_set, "_it", iteration,".txt")),
            quote = F,
            row.names = F
    )


    ## return fitted values
    if (return_model == TRUE) {
        return(dg_model)
    }
}
