#' Produces some basic analyses plots
#'
#' @param dataset_folder absolute path to the dataset folder, is created if non-existent
#' @param model_name name of the model that should be computed on the dataset
#'
#'  @import data.table
#'
#' @export
#'

dg_basic_analyses <- function(
    dataset_folder,
    model_name
){

    ggplot2::theme_set(ggplot2::theme_bw(base_size = 8))

    # load varlist
    load(file = file.path(dataset_folder, "data/fitness_dataset.RData"))
    # load parlist
    load(file = file.path(dataset_folder, model_name, "parameter_list.RData"))
    # load estimated models
    load(file = file.path(dataset_folder, model_name, "model_results.RData"))


    ## predict fitness from ddGs
    model_results[["variant_data"]] <- dg__fitness_from_model(
        dg_model0 = model_results[["avg_model"]],
        variant_data0 = model_results[["variant_data"]],
        varlist = varlist,
        parlist = parlist,
        calc_performance = FALSE,
        per_testset = FALSE,
        overwrite = FALSE,
        bfit_ddg0 = TRUE
    )

    ## plot inter-relationship between estimated ddg values and fitness values



    ## plot relationship between different fitness measurements if effects where only driven by folding
    idx <- combn(varlist[["no_abd_datasets"]],2)
    pidx <- 1
    for (i in seq_along(ncol(idx))) {
        plot_list[[pidx]] <- ggplot2::ggplot(model_results[["variant_data"]]) +
            ggplot2::geom_point(ggplot2::aes_string(
                x = paste0("f", idx[1, i], "_fitness"),
                y = paste0("f", idx[2, i], "_fitness"))) +
            ggplot2::geom_point(ggplot2::aes_string(
                x = paste0("f", idx[1, i], "_pred"),
                y = paste0("f", idx[2, i], "_pred")))
    }



}
