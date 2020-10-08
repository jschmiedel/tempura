#' Bootstrap procedure to estimate the uncertainty of free energies parameters of a model
#'
#'
#' @param dataset_folder absolute path to the dataset folder, is created if non-existent
#' @param model_name name of the model that should be computed on the dataset
#' @param iteration integer, to calucate a set of models and get reproducible results
#' @param return_model logical, if TRUE, returns a data.table with fitted values for all model parameters, default = FALSE
#' @param maxit integer, maximum number of iterations by optim algorithm, default = 1e4
#'
#' @return writes the boostrapped model parameters as .txt file to $dataset_folder/$model_name/tmp/dg_boostrap_$iteration
#' @import data.table
#' @import Matrix
#'
#' @export
#'
dg_bootstrap <- function(
    dataset_folder,
    model_name,
    iteration = 1,
    return_model = FALSE,
    maxit = 1e4
) {

    ## load preprocessed files
    # load varlist
    load(file = file.path(dataset_folder, "data/fitness_dataset.RData"))
    # load parlist
    load(file = file.path(dataset_folder, model_name, "parameter_list.RData"))
    # load collected models
    load(file = file.path(dataset_folder, model_name, "model_results.RData"))


    ## set seed for reproducible model results
    set.seed(iteration)


    ## get parameters from average dG model
    start_par <- model_results[["avg_model"]][, value]
    names(start_par) <- model_results[["avg_model"]][, parameter]

    # fix parameters (pro forma)
    for (i in seq_along(parlist[["fixed_par"]])) {
      start_par[names(start_par) == names(parlist[["fixed_par"]])[i]] <- parlist[["fixed_par"]][i]
    }
    # & enforce bounds
    start_par[start_par < parlist[["lower_bounds"]]] = parlist[["lower_bounds"]][start_par < parlist[["lower_bounds"]]]
    start_par[start_par > parlist[["upper_bounds"]]] = parlist[["upper_bounds"]][start_par > parlist[["upper_bounds"]]]


    ## bootstrap variants in dataset
    bootset <- sample(nrow(varlist[["variant_data"]]), replace = T)

    variant_data <- varlist[["variant_data"]][bootset]
    varxmut <- varlist[["varxmut"]][bootset,]
    # transpose varxmut matrix
    mutxvar <- Matrix::t(varxmut)


    ## sample fitness values
    fitval <- grep("fitness", names(variant_data), value = T)
    for (i in fitval) {
        variant_data[!is.na(get(i)), eval(i) := stats::rnorm(n = .N, mean = unlist(.SD[,1]), sd = unlist(.SD[,2])),
            .SDcols = c(i, paste0(strsplit(i, "_")[[1]][1], "_sigma"))]
    }


    ## write modified data (back) to varlist
    varlist[["variant_data"]] <- variant_data
    varlist[["varxmut"]] <- varxmut
    varlist[["mutxvar"]] <- mutxvar


    ## call optimzer
    dg_model <- dg__model_optim(
        start_par = start_par,
        parlist = parlist,
        varlist = varlist,
        maxit = maxit
    )


    ## convert output to data.table
    dg_boot <- data.table::data.table(
        t(dg_model[["par"]]),
        dg_model[["value"]],
        dg_model[["convergence"]],
        iteration
    )
    names(dg_boot) <- c(
        parlist[["par_names"]],
        "objective",
        "convergence",
        "iteration"
    )


    ## write to tmp/
    utils::write.table(dg_boot,
        file = file.path(dataset_folder, model_name, paste0("tmp/dg_bootstrap_it", iteration,".txt")),
            quote = F,
            row.names = F
    )


    ## return fitted values
    if (return_model == TRUE) {
        return(dg_boot)
    }
}
