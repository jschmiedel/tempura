#' Produces some basic analyses plots
#'
#' @param dataset_folder absolute path to the dataset folder, is created if non-existent
#' @param model_name name of the model that should be computed on the dataset
#' @param color_type c("", "residue_type", "ligand_dist", "rsa") how to color variants; any other than the default "" will exclude variants with multiple mutations and the wild-type variant
#' @param max_variants maximual number of variants to plot, default = 3e3, will be downsampled if # variants > max_variants
#'
#' @import data.table
#'
#' @export
#'

dg_basic_analyses <- function(
    dataset_folder,
    model_name,
    color_type = "",
    max_variants = 3e3
){

    ggplot2::theme_set(ggplot2::theme_bw(base_size = 8))

    # load varlist
    load(file = file.path(dataset_folder, "data/fitness_dataset.RData"))
    # load parlist
    load(file = file.path(dataset_folder, model_name, "parameter_list.RData"))
    # load estimated models
    load(file = file.path(dataset_folder, model_name, "model_results.RData"))


    ## predict fitness from ddGs, esp. binding fitness assuming either ddG_folding = 0 or ddG_binding = 0
    X <- dg__fitness_from_model(
        dg_model0 = model_results[["avg_model"]],
        variant_data0 = model_results[["variant_data"]],
        varlist = varlist,
        parlist = parlist,
        calc_performance = FALSE,
        per_testset = FALSE,
        overwrite = TRUE,
        bfit_ddg0 = TRUE
    )
    model_results[["variant_data"]] <- X[["variant_data"]]


    ## plot inter-relationship between estimated ddg values and fitness values
    # compute ddGs for all variants
    if (parlist[["no_folded_states"]] == 1) {
        model_results[["variant_data"]][, f_ddg := as.numeric(varlist[["varxmut"]] %*% model_results[["avg_model"]][grep("f_ddg", parameter), value])]
    } else {
        model_results[["variant_data"]][, fA_ddg := as.numeric(varlist[["varxmut"]] %*% model_results[["avg_model"]][grep("fA_ddg", parameter), value])]
        model_results[["variant_data"]][, fB_ddg := as.numeric(varlist[["varxmut"]] %*% model_results[["avg_model"]][grep("fB_ddg", parameter), value])]
    }
    model_results[["variant_data"]][, b_ddg := as.numeric(varlist[["varxmut"]] %*% model_results[["avg_model"]][grep("b_ddg", parameter), value])]


    vd <- model_results[["variant_data"]][!is.na(f1_fitness) & !is.na(b1_fitness)]
    if (color_type == "") {
        vd[, color_type := "all vars"]
    }
    am <- model_results[["avg_model"]]

    if (parlist[["no_folded_states"]] == 1) {
        df_dist <- ggplot2::ggplot(data = vd) +
            ggplot2::geom_density(ggplot2::aes(x = f_ddg + am[grep("^f1_dgwt", parameter), value], color = color_type)) +
            ggplot2::geom_vline(xintercept = am[grep("^f1_dgwt", parameter), value],linetype = 2) +
            ggplot2::scale_color_brewer(palette = "Set1") +
            ggplot2::theme(legend.position = c(0.75, 0.75)) +
            ggplot2::labs(x = "dG folding", color = "")

        df_ffitness <- ggplot2::ggplot(data = vd) +
            ggplot2::geom_density2d(ggplot2::aes(x = f_ddg + am[grep("^f1_dgwt", parameter), value],
                  y = f1_fitness), color = 'black') +
            ggplot2::geom_point(ggplot2::aes(x = f_ddg + am[grep("^f1_dgwt", parameter), value],
                  y = f1_fitness, color = color_type)) +
            ggplot2::geom_line(ggplot2::aes(x = f_ddg + am[grep("^f1_dgwt", parameter), value],
                  y = f1_pred),color = "black") +
            # scale_y_continuous(breaks = seq(-6, 1, 1)) +
            # scale_x_continuous(breaks = seq(-10,6,2)) +
            ggplot2::geom_hline(yintercept = am[grep("^f1_fit", parameter), c(value)],linetype = 2) +
            ggplot2::geom_vline(xintercept = am[grep("^f1_dgwt", parameter), value],linetype = 2) +
            ggplot2::labs(x = "dG folding", y = "predicted folding fitness", color = "") +
            ggplot2::theme(legend.position = "none") +
            ggplot2::scale_color_brewer(palette = "Set1")

        df_bfitness <- ggplot2::ggplot(data = vd) +
            ggplot2::geom_density2d(ggplot2::aes(x = f_ddg + am[grep("^bf1_dgwt", parameter), value],
                  y = b1_fitness), color = 'black') +
            ggplot2::geom_point(ggplot2::aes(x = f_ddg + am[grep("^bf1_dgwt", parameter), value],
                  y = b1_fitness, color = color_type)) +
            ggplot2::geom_line(ggplot2::aes(x = f_ddg + am[grep("^bf1_dgwt", parameter), value],
                   y = b1_pred_bddg0), color = "black") +
            # ggplot2::scale_y_continuous(breaks = seq(-6, 1, 1)) +
            # ggplot2::scale_x_continuous(breaks = seq(-10,6,2)) +
            ggplot2::geom_hline(yintercept = am[grep("^b1_fit", parameter), c(value)],linetype = 2) +
            ggplot2::geom_vline(xintercept = am[grep("^bf1_dgwt", parameter), value],linetype = 2) +
            ggplot2::labs(x = "dG folding", y = "binding fitness", color = "") +
            ggplot2::theme(legend.position = "none") +
            ggplot2::scale_color_brewer(palette = "Set1")

        df_db <- ggplot2::ggplot(vd,
              ggplot2::aes(x = f_ddg + am[grep("^f1_dgwt", parameter), value],
                  y = b_ddg + am[grep("^b1_dgwt", parameter), value], color = color_type)) +
            ggplot2::geom_density2d() +
            ggplot2::geom_point(alpha=0.2) +
            # ggplot2::scale_x_continuous(breaks = seq(-6,6,1)) +
            # ggplot2::scale_y_continuous(breaks = seq(-6,6,1)) +
            ggplot2::scale_color_brewer(palette = "Set1") +
            ggplot2::geom_vline(xintercept = am[grep("^f1_dgwt", parameter), value],linetype = 2) +
            ggplot2::geom_hline(yintercept = am[grep("^b1_dgwt", parameter), value],linetype = 2) +
            ggplot2::labs(x = "dG folding", y = "dG binding", color = "") +
            ggplot2::theme(legend.position = "none")
    }

    db_dist <- ggplot2::ggplot() +
        ggplot2::geom_density(data = vd,
          ggplot2::aes(x = b_ddg + am[grep("^b1_dgwt", parameter), value], color = color_type)) +
        ggplot2::geom_vline(xintercept = am[grep("^b1_dgwt", parameter), value],linetype = 2) +
        # ggplot2::scale_x_continuous(breaks = seq(-6,6,1)) +
        ggplot2::scale_color_brewer(palette = "Set1") +
        ggplot2::coord_flip() +
        ggplot2::labs(x = "dG binding") +
        ggplot2::theme(legend.position = "none")

    db_bfitness <- ggplot2::ggplot(data = vd) +
        ggplot2::geom_density2d(ggplot2::aes(x = b_ddg + am[grep("^b1_dgwt", parameter), value],
              y = b1_fitness), color = 'black') +
        ggplot2::geom_point(ggplot2::aes(x = b_ddg + am[grep("^b1_dgwt", parameter), value],
              y = b1_fitness, color = color_type)) +
        ggplot2::geom_line(ggplot2::aes(x = b_ddg + am[grep("^b1_dgwt", parameter), value],
              y = b1_pred_fddg0),color = "black") +
        # ggplot2::scale_y_continuous(breaks = seq(-6, 1, 1)) +
        # ggplot2::scale_x_continuous(breaks = seq(-6,6,1)) +
        ggplot2::scale_color_brewer(palette = "Set1") +
        ggplot2::geom_hline(yintercept = am[grep("^b1_fit", parameter), c(value)],linetype = 2) +
        ggplot2::geom_vline(xintercept = am[grep("^b1_dgwt", parameter), value],linetype = 2) +
        ggplot2::labs(x = "dG binding", y = "binding fitness", color = "") +
        ggplot2::theme(legend.position = "none")

    if (parlist[["no_folded_states"]] == 1) {
        p <- gridExtra::grid.arrange(df_dist, df_ffitness, df_bfitness, df_db, db_dist, db_bfitness,
                       nrow = 2)
        ggplot2::ggsave(p,
            file = file.path(dataset_folder, model_name, "results/dG_relationship_fitness.pdf"),
            width = 10,
            height = 7)
    }





    ## plot relationship between different fitness measurements if effects where only driven by folding
    plot_list = list()
    pidx <- 1

    # abundance fitness datasets against one another
    idx <- utils::combn(varlist[["no_abd_datasets"]],2)
    for (i in 1:ncol(idx)) {
        plot_list[[pidx]] <- ggplot2::ggplot(data = model_results[["variant_data"]][
                    !is.na(get(paste0("f", idx[1, i], "_fitness"))) &
                    !is.na(get(paste0("f", idx[2, i], "_fitness")))]) +
            ggplot2::geom_point(ggplot2::aes_string(
                x = paste0("f", idx[1, i], "_fitness"),
                y = paste0("f", idx[2, i], "_fitness"))) +
            ggplot2::geom_line(ggplot2::aes_string(
                x = paste0("f", idx[1, i], "_pred"),
                y = paste0("f", idx[2, i], "_pred")), color = "red")
        pidx <- pidx + 1
    }

    # abundance versus binding fitness datasets
    for (fi in 1:varlist[["no_abd_datasets"]]) {
        for (bi in 1:varlist[["no_bind_datasets"]]) {
            plot_list[[pidx]] <- ggplot2::ggplot(model_results[["variant_data"]][
                    !is.na(get(paste0("b", bi, "_fitness"))) &
                    !is.na(get(paste0("f", fi, "_fitness")))]) +
                ggplot2::geom_point(ggplot2::aes_string(
                    x = paste0("b", bi, "_fitness"),
                    y = paste0("f", fi, "_fitness"))) +
                ggplot2::geom_line(ggplot2::aes_string(
                    x = paste0("b", bi, "_pred_bddg0"),
                    y = paste0("f", fi, "_pred")), color = "red")
            pidx <- pidx + 1
        }
    }

    # binding fitness datasets against one another
    if (varlist[["no_bind_datasets"]] > 1) {
        idx <- utils::combn(varlist[["no_bind_datasets"]],2)

        for (i in seq_along(ncol(idx))) {
            plot_list[[pidx]] <- ggplot2::ggplot(model_results[["variant_data"]][
                    !is.na(get(paste0("b", idx[1, i], "_fitness"))) &
                    !is.na(get(paste0("b", idx[2, i], "_fitness")))]) +
                ggplot2::geom_point(ggplot2::aes_string(
                    x = paste0("b", idx[1, i], "_fitness"),
                    y = paste0("b", idx[2, i], "_fitness"))) +
                ggplot2::geom_line(ggplot2::aes_string(
                    x = paste0("b", idx[1, i], "_pred"),
                    y = paste0("b", idx[2, i], "_pred")), color = "red")
            pidx <- pidx + 1
        }
    }

    p <- gridExtra::grid.arrange(grobs = plot_list,
        nrow = ceiling(sqrt(length(plot_list))),
        ncol = ceiling(sqrt(length(plot_list))))
    ggplot2::ggsave(p,
        file = file.path(dataset_folder, model_name, "results/fitness_scatter_dg_relationship.pdf"),
        width = 4 * ceiling(sqrt(length(plot_list))),
        height = 3.5 * ceiling(sqrt(length(plot_list))))



}
