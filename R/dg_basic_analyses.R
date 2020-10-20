#' Produces some basic analyses plots
#'
#' @param dataset_folder absolute path to the dataset folder, is created if non-existent
#' @param model_name name of the model that should be computed on the dataset
#' @param color_type variable with which to color variants; can be Nmut or any column_name from structural_properties file; any other than the default "" will exclude variants with multiple mutations and the wild-type variant
#' @param max_variants maximual number of variants to plot, default = 3e3, will be downsampled if # variants > max_variants
#' @param datasets_ab index of datasets to use (in case multiple abundance or binding datasets were supplied)
#'
#' @import data.table
#'
#' @export
#'

dg_basic_analyses <- function(
    dataset_folder,
    model_name,
    color_type = "",
    max_variants = 3e3,
    datasets_ab = c(1, 1)
){

    varlist <- parlist <- f_ddg <- parameter <- value <- fA_ddg <- fB_ddg <- b_ddg <- aa_subs <- Pos <- NULL

    ggplot2::theme_set(ggplot2::theme_bw(base_size = 8))
    col_purple = "#9161A8"
    col_orange = "#F7941E"

    # load varlist
    load(file = file.path(dataset_folder, "data/fitness_dataset.RData"))
    # load parlist
    load(file = file.path(dataset_folder, model_name, "parameter_list.RData"))
    # load estimated models
    load(file = file.path(dataset_folder, model_name, "model_results.RData"))

    if (file.exists(file.path(dataset_folder, "data/structural_properties.txt"))) {
        structural_properties <- fread(file.path(dataset_folder, "data/structural_properties.txt"))
    }

    ## predict fitness from ddGs, esp. binding fitness assuming either ddG_folding = 0 or ddG_binding = 0
    X <- dg__fitness_from_model(
        dg_model0 = model_results[["avg_model"]],
        variant_data0 = model_results[["variant_data"]],
        varlist = varlist,
        parlist = parlist,
        calc_performance = FALSE,
        per_testset = FALSE,
        overwrite = TRUE,
        predfit_ddgvar = TRUE
    )
    model_results[["variant_data"]] <- X[["variant_data"]]

    ## compute ddGs for all variants
    if (parlist[["no_folded_states"]] == 1) {
        model_results[["variant_data"]][, f_ddg := as.numeric(varlist[["varxmut"]] %*% model_results[["avg_model"]][grep("f_ddg", parameter), value])]
    } else {
        model_results[["variant_data"]][, fA_ddg := as.numeric(varlist[["varxmut"]] %*% model_results[["avg_model"]][grep("fA_ddg", parameter), value])]
        model_results[["variant_data"]][, fB_ddg := as.numeric(varlist[["varxmut"]] %*% model_results[["avg_model"]][grep("fB_ddg", parameter), value])]
    }
    model_results[["variant_data"]][, b_ddg := as.numeric(varlist[["varxmut"]] %*% model_results[["avg_model"]][grep("b_ddg", parameter), value])]

    ## define helper pars and include structural properties
    vd <- model_results[["variant_data"]]
    if (color_type == "") {
        vd[, color_type := factor("all vars")]
    } else if (color_type == "Nmut") {
        vd[grepl("[0-9]", aa_subs), color_type := factor(length(grep("_", aa_subs)) + 1), aa_subs]
    } else if (exists("structural_properties") == TRUE) {
        vd[!grepl("_", aa_subs) & grepl("[0-9]", aa_subs), Pos := as.integer(paste0(strsplit(aa_subs,"")[[1]][1:(nchar(aa_subs)-1)], collapse = "")), aa_subs]
        vd <- merge(vd,
            structural_properties[, list(Pos, color_type = factor(unlist(.SD[,1]))), .SDcols = color_type],
            by = "Pos")
    } else {
        print("error: no structural properties file provided to color variants")
    }
    am <- model_results[["avg_model"]]


    ## plot inter-relationship between estimated ddg values and fitness values
    # plot f_ddg distribution and relationship with f_fitness, b_ddg and b_fitness
    vd_plot <- vd[!is.na(get(paste0("f", datasets_ab[1], "_fitness"))) &
                  !is.na(get(paste0("b", datasets_ab[2], "_fitness"))) &
                  !is.na(color_type)]

    if (parlist[["no_folded_states"]] == 1) {
        df_dist <- ggplot2::ggplot(data = vd_plot)
        if (is.factor(vd_plot$color_type)) {
            df_dist <- df_dist +
                ggplot2::geom_density(ggplot2::aes(x = f_ddg + am[grep(paste0("^f", datasets_ab[1], "_dgwt"), parameter), value],
                                            color = color_type)) +
                ggplot2::scale_color_brewer(palette = "Set1")
        } else {
            df_dist <- df_dist +
                ggplot2::geom_density(ggplot2::aes(x = f_ddg + am[grep(paste0("^f", datasets_ab[1], "_dgwt"), parameter), value]), color = 'black')
        }
       df_dist <- df_dist +
            ggplot2::geom_vline(xintercept = am[grep(paste0("^f", datasets_ab[1], "_dgwt"), parameter), value],linetype = 2) +
            ggplot2::theme(legend.position = c(0.75, 0.75)) +
            ggplot2::labs(x = "dG folding", color = "")

        df_ffitness <- ggplot2::ggplot(data = vd_plot,
                ggplot2::aes(x = f_ddg + am[grep(paste0("^f", datasets_ab[1], "_dgwt"), parameter), value])) +
            ggplot2::geom_density2d(ggplot2::aes_string(y = paste0("f", datasets_ab[1], "_fitness")), color = 'black') +
            ggplot2::geom_point(ggplot2::aes_string(y = paste0("f", datasets_ab[1], "_fitness"), color = "color_type")) +
            ggplot2::geom_line(ggplot2::aes_string(y = paste0("f", datasets_ab[1], "_pred")), color = "black") +
            # scale_y_continuous(breaks = seq(-6, 1, 1)) +
            # scale_x_continuous(breaks = seq(-10,6,2)) +
            ggplot2::geom_hline(yintercept = am[grep(paste0("^f", datasets_ab[1], "_fit"), parameter), c(value)],linetype = 2) +
            ggplot2::geom_vline(xintercept = am[grep(paste0("^f", datasets_ab[1], "_dgwt"), parameter), value],linetype = 2) +
            ggplot2::labs(x = "dG folding",
                y = paste0("folding fitness ", datasets_ab[1]),
                color = "") +
            ggplot2::theme(legend.position = "none")
        if (is.factor(vd_plot$color_type)) {
            df_ffitness <- df_ffitness + ggplot2::scale_color_brewer(palette = "Set1")
        } else {
            df_ffitness <- df_ffitness + ggplot2::scale_color_gradient(low = col_orange, high = col_purple)
        }


        df_bfitness <- ggplot2::ggplot(data = vd_plot,
                ggplot2::aes(x = f_ddg + am[grep(paste0("^bf", datasets_ab[2], "_dgwt"), parameter), value])) +
            ggplot2::geom_density2d(ggplot2::aes_string(y = paste0("b", datasets_ab[2], "_fitness")), color = 'black') +
            ggplot2::geom_point(ggplot2::aes_string(y = paste0("b", datasets_ab[2], "_fitness"), color = "color_type")) +
            ggplot2::geom_line(ggplot2::aes_string(y = paste0("b", datasets_ab[2], "_pred_bddg0")), color = "black") +
            # ggplot2::scale_y_continuous(breaks = seq(-6, 1, 1)) +
            # ggplot2::scale_x_continuous(breaks = seq(-10,6,2)) +
            ggplot2::geom_hline(yintercept = am[grep(paste0("^b", datasets_ab[2], "_fit"), parameter), c(value)],linetype = 2) +
            ggplot2::geom_vline(xintercept = am[grep(paste0("^bf", datasets_ab[2], "_dgwt"), parameter), value],linetype = 2) +
            ggplot2::labs(x = "dG folding",
                y = paste0("binding fitness ", datasets_ab[2]),
                color = "") +
            ggplot2::theme(legend.position = "none")
        if (is.factor(vd_plot$color_type)) {
            df_bfitness <- df_bfitness + ggplot2::scale_color_brewer(palette = "Set1")
        } else {
            df_bfitness <- df_bfitness + ggplot2::scale_color_gradient(low = col_orange, high = col_purple)
        }

        df_db <- ggplot2::ggplot(vd_plot,
              ggplot2::aes(x = f_ddg + am[grep(paste0("^f", datasets_ab[1], "_dgwt"), parameter), value],
                  y = b_ddg + am[grep(paste0("^b", datasets_ab[2], "_dgwt"), parameter), value], color = color_type)) +
            ggplot2::geom_density2d() +
            ggplot2::geom_point(alpha = 0.4) +
            # ggplot2::scale_x_continuous(breaks = seq(-6,6,1)) +
            # ggplot2::scale_y_continuous(breaks = seq(-6,6,1)) +
            ggplot2::geom_vline(xintercept = am[grep(paste0("^f", datasets_ab[1], "_dgwt"), parameter), value],linetype = 2) +
            ggplot2::geom_hline(yintercept = am[grep(paste0("^b", datasets_ab[2], "_dgwt"), parameter), value],linetype = 2) +
            ggplot2::labs(x = "dG folding",
                y = "dG binding",
                color = "") +
            ggplot2::theme(legend.position = "none")
        if (is.factor(vd_plot$color_type)) {
            df_db <- df_db + ggplot2::scale_color_brewer(palette = "Set1")
        } else {
            df_db <- df_db + ggplot2::scale_color_gradient(low = col_orange, high = col_purple)
        }



    } else if (parlist[["no_folded_states"]] == 2) {
        dfA_dist <- ggplot2::ggplot(data = vd_plot)
        if (is.factor(vd_plot$color_type)) {
            dfA_dist <- dfA_dist +
                ggplot2::geom_density(
                    ggplot2::aes(x = fA_ddg + am[grep(paste0("^fA", datasets_ab[1], "_dgwt"), parameter), value],
                        color = color_type)) +
                ggplot2::scale_color_brewer(palette = "Set1")
        } else {
            dfA_dist <- dfA_dist +
                ggplot2::geom_density(
                    ggplot2::aes(x = fA_ddg + am[grep(paste0("^fA", datasets_ab[1], "_dgwt"), parameter), value]),
                        color = 'black')
        }
        dfA_dist <- dfA_dist +
            ggplot2::geom_vline(xintercept = am[grep(paste0("^fA", datasets_ab[1], "_dgwt"), parameter), value], linetype = 2) +
            ggplot2::theme(legend.position = c(0.75, 0.75)) +
            ggplot2::labs(x = "dG folding state A", color = "")

        dfB_dist <- ggplot2::ggplot(data = vd_plot)
        if (is.factor(vd_plot$color_type)) {
            dfB_dist <- dfB_dist +
                ggplot2::geom_density(
                    ggplot2::aes(x = fB_ddg + am[grep(paste0("^fB", datasets_ab[1], "_dgwt"), parameter), value],
                        color = color_type)) +
                ggplot2::scale_color_brewer(palette = "Set1")
        } else {
            dfB_dist <- dfB_dist +
                ggplot2::geom_density(
                    ggplot2::aes(x = fB_ddg + am[grep(paste0("^fB", datasets_ab[1], "_dgwt"), parameter), value]),
                        color = 'black')
        }
        dfB_dist <- dfB_dist +
            ggplot2::geom_vline(xintercept = am[grep(paste0("^fB", datasets_ab[1], "_dgwt"), parameter), value], linetype = 2) +
            ggplot2::theme(legend.position = "none") +
            ggplot2::coord_flip() +
            ggplot2::labs(x = "dG folding state B", color = "")

        dfA_ffitness <- ggplot2::ggplot(data = vd_plot,
                ggplot2::aes(x = fA_ddg + am[grep(paste0("^fA", datasets_ab[1], "_dgwt"), parameter), value])) +
            ggplot2::geom_density2d(ggplot2::aes_string(x = paste0('fA_ddg + am[grep(paste0("^fA", datasets_ab[1], "_dgwt"), parameter), value]'),
                y = paste0("f", datasets_ab[1], "_fitness")), color = 'black') +
            ggplot2::geom_point(ggplot2::aes_string(y = paste0("f", datasets_ab[1], "_fitness"), color = "color_type")) +
            ggplot2::geom_line(ggplot2::aes_string(y = paste0("f", datasets_ab[1], "_pred_fBasfA")), color = "black") +
            ggplot2::geom_line(ggplot2::aes_string(y = paste0("f", datasets_ab[1], "_pred_fBddg0")), color = "black", linetype = 3) +
            ggplot2::geom_hline(yintercept = am[grep(paste0("^f", datasets_ab[1], "_fit"), parameter), c(value)], linetype = 2) +
            ggplot2::geom_vline(xintercept = am[grep(paste0("^fA", datasets_ab[1], "_dgwt"), parameter), value], linetype = 2) +
            ggplot2::labs(x = "dG folding state A",
                y = paste0("folding fitness ", datasets_ab[1]),
                color = "") +
            ggplot2::theme(legend.position = "none")
        if (is.factor(vd_plot$color_type)) {
            dfA_ffitness <- dfA_ffitness + ggplot2::scale_color_brewer(palette = "Set1")
        } else {
            dfA_ffitness <- dfA_ffitness + ggplot2::scale_color_gradient(low = col_orange, high = col_purple)
        }

        dfB_ffitness <- ggplot2::ggplot(data = vd_plot,
                ggplot2::aes(x = fB_ddg + am[grep(paste0("^fB", datasets_ab[1], "_dgwt"), parameter), value])) +
            ggplot2::geom_density2d(ggplot2::aes_string(y = paste0("f", datasets_ab[1], "_fitness")), color = 'black') +
            ggplot2::geom_point(ggplot2::aes_string(y = paste0("f", datasets_ab[1], "_fitness"), color = "color_type")) +
            ggplot2::geom_line(
                ggplot2::aes_string(y = paste0("f", datasets_ab[1], "_pred_fAasfB")),color = "black") +
            ggplot2::geom_line(ggplot2::aes_string(y = paste0("f", datasets_ab[1], "_pred_fAddg0")), color = "black", linetype = 3) +
            ggplot2::geom_hline(yintercept = am[grep(paste0("^f", datasets_ab[1], "_fit"), parameter), c(value)], linetype = 2) +
            ggplot2::geom_vline(xintercept = am[grep(paste0("^fB", datasets_ab[1], "_dgwt"), parameter), value], linetype = 2) +
            ggplot2::labs(x = "dG folding state B",
                y = paste0("folding fitness ", datasets_ab[1]),
                color = "") +
            ggplot2::theme(legend.position = "none")
        if (is.factor(vd_plot$color_type)) {
            dfB_ffitness <- dfB_ffitness + ggplot2::scale_color_brewer(palette = "Set1")
        } else {
            dfB_ffitness <- dfB_ffitness + ggplot2::scale_color_gradient(low = col_orange, high = col_purple)
        }

        dfA_dfB <- ggplot2::ggplot(vd_plot,
              ggplot2::aes(x = fA_ddg + am[grep(paste0("^fA", datasets_ab[1], "_dgwt"), parameter), value],
                  y = fB_ddg + am[grep(paste0("^fB", datasets_ab[1], "_dgwt"), parameter), value], color = color_type)) +
            ggplot2::geom_density2d() +
            ggplot2::geom_point(alpha = 0.4) +
            ggplot2::geom_vline(xintercept = am[grep(paste0("^fA", datasets_ab[1], "_dgwt"), parameter), value], linetype = 2) +
            ggplot2::geom_hline(yintercept = am[grep(paste0("^fB", datasets_ab[1], "_dgwt"), parameter), value], linetype = 2) +
            ggplot2::geom_abline(color = "red",
                intercept = am[grep(paste0("^fB", datasets_ab[1], "_dgwt"), parameter), value] -
                            am[grep(paste0("^fA", datasets_ab[1], "_dgwt"), parameter), value]) +
            ggplot2::labs(x = "dG folding state A",
                y = "dG folding state B",
                color = "") +
            ggplot2::theme(legend.position = "none")
        if (is.factor(vd_plot$color_type)) {
            dfA_dfB <- dfA_dfB + ggplot2::scale_color_brewer(palette = "Set1")
        } else {
            dfA_dfB <- dfA_dfB + ggplot2::scale_color_gradient(low = col_orange, high = col_purple)
        }



        dfA_db <- ggplot2::ggplot(vd_plot,
              ggplot2::aes(x = fA_ddg + am[grep(paste0("^fA", datasets_ab[1], "_dgwt"), parameter), value],
                  y = b_ddg + am[grep(paste0("^b", datasets_ab[2], "_dgwt"), parameter), value], color = color_type)) +
            ggplot2::geom_density2d() +
            ggplot2::geom_point(alpha = 0.4) +
            ggplot2::geom_vline(xintercept = am[grep(paste0("^fA", datasets_ab[1], "_dgwt"), parameter), value],linetype = 2) +
            ggplot2::geom_hline(yintercept = am[grep(paste0("^b", datasets_ab[2], "_dgwt"), parameter), value],linetype = 2) +
            ggplot2::labs(x = "dG folding state A",
                y = "dG binding",
                color = "") +
            ggplot2::theme(legend.position = "none")
        if (is.factor(vd_plot$color_type)) {
            dfA_db <- dfA_db + ggplot2::scale_color_brewer(palette = "Set1")
        } else {
            dfA_db <- dfA_db + ggplot2::scale_color_gradient(low = col_orange, high = col_purple)
        }

        dfB_db <- ggplot2::ggplot(vd_plot,
              ggplot2::aes(x = fB_ddg + am[grep(paste0("^fB", datasets_ab[1], "_dgwt"), parameter), value],
                  y = b_ddg + am[grep(paste0("^b", datasets_ab[2], "_dgwt"), parameter), value], color = color_type)) +
            ggplot2::geom_density2d() +
            ggplot2::geom_point(alpha = 0.4) +
            ggplot2::geom_vline(xintercept = am[grep(paste0("^fB", datasets_ab[1], "_dgwt"), parameter), value],linetype = 2) +
            ggplot2::geom_hline(yintercept = am[grep(paste0("^b", datasets_ab[2], "_dgwt"), parameter), value],linetype = 2) +
            ggplot2::labs(x = "dG folding state B",
                y = "dG binding",
                color = "") +
            ggplot2::theme(legend.position = "none")
        if (is.factor(vd_plot$color_type)) {
            dfB_db <- dfB_db + ggplot2::scale_color_brewer(palette = "Set1")
        } else {
            dfB_db <- dfB_db + ggplot2::scale_color_gradient(low = col_orange, high = col_purple)
        }

        dfA_bfitness <- ggplot2::ggplot(data = vd_plot,
                ggplot2::aes(x = fA_ddg + am[grep(paste0("^bfA", datasets_ab[2], "_dgwt"), parameter), value])) +
            ggplot2::geom_density2d(ggplot2::aes_string(y = paste0("b", datasets_ab[2], "_fitness")), color = 'black') +
            ggplot2::geom_point(ggplot2::aes_string(y = paste0("b", datasets_ab[2], "_fitness"), color = "color_type")) +
            ggplot2::geom_line(ggplot2::aes_string(y = paste0("b", datasets_ab[2], "_pred_fBasfA")),color = "black") +
            ggplot2::geom_line(ggplot2::aes_string(y = paste0("b", datasets_ab[2], "_pred_fBddg0")), color = "black", linetype = 3) +
            ggplot2::geom_hline(yintercept = am[grep(paste0("^b", datasets_ab[2], "_fit"), parameter), c(value)],linetype = 2) +
            ggplot2::geom_vline(xintercept = am[grep(paste0("^bfA", datasets_ab[2], "_dgwt"), parameter), value],linetype = 2) +
            ggplot2::labs(x = "dG folding state A",
                y = paste0("binding fitness ", datasets_ab[2]),
                color = "") +
            ggplot2::theme(legend.position = "none")
        if (is.factor(vd_plot$color_type)) {
            dfA_bfitness <- dfA_bfitness + ggplot2::scale_color_brewer(palette = "Set1")
        } else {
            dfA_bfitness <- dfA_bfitness + ggplot2::scale_color_gradient(low = col_orange, high = col_purple)
        }

        dfB_bfitness <- ggplot2::ggplot(data = vd_plot,
                ggplot2::aes(x = fB_ddg + am[grep(paste0("^bfB", datasets_ab[2], "_dgwt"), parameter), value])) +
            ggplot2::geom_density2d(ggplot2::aes_string(y = paste0("b", datasets_ab[2], "_fitness")), color = 'black') +
            ggplot2::geom_point(ggplot2::aes_string(y = paste0("b", datasets_ab[2], "_fitness"), color = "color_type")) +
            ggplot2::geom_line(
                ggplot2::aes_string(y = paste0("b", datasets_ab[2], "_pred_fAasfB")),color = "black") +
            ggplot2::geom_line(ggplot2::aes_string(y = paste0("b", datasets_ab[2], "_pred_fAddg0")), color = "black", linetype = 3) +
            ggplot2::geom_hline(yintercept = am[grep(paste0("^b", datasets_ab[2], "_fit"), parameter), c(value)],linetype = 2) +
            ggplot2::geom_vline(xintercept = am[grep(paste0("^bfB", datasets_ab[2], "_dgwt"), parameter), value],linetype = 2) +
            ggplot2::labs(x = "dG folding state B",
                y = paste0("binding fitness ", datasets_ab[2]),
                color = "") +
            ggplot2::theme(legend.position = "none")
        if (is.factor(vd_plot$color_type)) {
            dfB_bfitness <- dfB_bfitness + ggplot2::scale_color_brewer(palette = "Set1")
        } else {
            dfB_bfitness <- dfB_bfitness + ggplot2::scale_color_gradient(low = col_orange, high = col_purple)
        }

    } else {
        print(paste0(parlist[["no_folded_states"]], " folded states currently not supported"))
    }

    # plot b_ddg distribution and relationship with b1_fitness
    db_dist <- ggplot2::ggplot()
    if (is.factor(vd_plot$color_type)) {
        db_dist <- db_dist +
            ggplot2::geom_density(data = vd_plot,
                ggplot2::aes(x = b_ddg + am[grep(paste0("^b", datasets_ab[2], "_dgwt"), parameter), value],
                    color = color_type)) +
            ggplot2::scale_color_brewer(palette = "Set1")
    } else {
        db_dist <- db_dist +
            ggplot2::geom_density(data = vd_plot,
                ggplot2::aes(x = b_ddg + am[grep(paste0("^b", datasets_ab[2], "_dgwt"), parameter), value]), color = 'black')
    }
    db_dist <- db_dist +
        ggplot2::geom_vline(xintercept = am[grep(paste0("^b", datasets_ab[2], "_dgwt"), parameter), value],linetype = 2) +
        # ggplot2::scale_x_continuous(breaks = seq(-6,6,1)) +
        ggplot2::coord_flip() +
        ggplot2::labs(x = "dG binding", color = "") +
        ggplot2::theme(legend.position = c(0.75, 0.25))


    db_bfitness <- ggplot2::ggplot(data = vd_plot,
            ggplot2::aes(x = b_ddg + am[grep(paste0("^b", datasets_ab[2], "_dgwt"), parameter), value])) +
        ggplot2::geom_density2d(ggplot2::aes_string(y = paste0("b", datasets_ab[2], "_fitness")), color = 'black') +
        ggplot2::geom_point(ggplot2::aes_string(y = paste0("b", datasets_ab[2], "_fitness"), color = "color_type")) +
        ggplot2::geom_line(ggplot2::aes_string(y = paste0("b", datasets_ab[2], "_pred_fddg0")), color = "black") +
        # ggplot2::scale_y_continuous(breaks = seq(-6, 1, 1)) +
        # ggplot2::scale_x_continuous(breaks = seq(-6,6,1))
        ggplot2::geom_hline(yintercept = am[grep(paste0("^b", datasets_ab[2], "_fit"), parameter), c(value)],linetype = 2) +
        ggplot2::geom_vline(xintercept = am[grep(paste0("^b", datasets_ab[2], "_dgwt"), parameter), value],linetype = 2) +
        ggplot2::labs(x = "dG binding",
            y = paste0("binding fitness ", datasets_ab[2]),
            color = "") +
        ggplot2::theme(legend.position = "none")
    if (is.factor(vd_plot$color_type)) {
        db_bfitness <- db_bfitness + ggplot2::scale_color_brewer(palette = "Set1")
    } else {
        db_bfitness <- db_bfitness + ggplot2::scale_color_gradient(low = col_orange, high = col_purple)
    }


    # plot
    if (parlist[["no_folded_states"]] == 1) {
        # plot folding and binding
        p <- gridExtra::grid.arrange(df_dist, df_ffitness, df_bfitness,
                                     df_db, db_dist, db_bfitness,
            nrow = 2)
        ggplot2::ggsave(p,
            file = file.path(dataset_folder, model_name,
                paste0("results/dG_relationship_fitness_f", datasets_ab[1], "b", datasets_ab[2], "_", color_type, ".pdf")),
            width = 10,
            height = 7)
    } else if (parlist[["no_folded_states"]] == 2) {
        # plot folding
        p <- gridExtra::grid.arrange(dfA_dist, dfA_ffitness, dfB_ffitness,
                                     dfA_dfB, dfB_dist,
            nrow = 2)
        ggplot2::ggsave(p,
            file = file.path(dataset_folder, model_name,
                paste0("results/dG_relationship_fitness_fA", datasets_ab[1], "_fB", datasets_ab[2], "_", color_type, ".pdf")),
            width = 10,
            height = 7)

        # plot binding
        p <- gridExtra::grid.arrange(db_dist, dfA_db, dfB_db,
                                     dfA_bfitness, dfB_bfitness, db_bfitness,
            nrow = 2)
        ggplot2::ggsave(p,
            file = file.path(dataset_folder, model_name,
                paste0("results/dG_relationship_b", datasets_ab[1], "_fA", datasets_ab[1], "_fB", datasets_ab[2], "_", color_type, ".pdf")),
            width = 10,
            height = 7)
    }




    ## plot relationship between different fitness measurements if effects where only driven by folding
    plot_list = list()
    pidx <- 1

    # abundance fitness datasets against one another
    if (varlist[["no_abd_datasets"]] > 1) {
        idx <- utils::combn(varlist[["no_abd_datasets"]],2)
        for (i in 1:ncol(idx)) {
            plot_list[[pidx]] <- ggplot2::ggplot(data = vd[
                        !is.na(get(paste0("f", idx[1, i], "_fitness"))) &
                        !is.na(get(paste0("f", idx[2, i], "_fitness"))) &
                        !is.na(color_type)]) +
                ggplot2::geom_point(ggplot2::aes_string(
                    x = paste0("f", idx[1, i], "_fitness"),
                    y = paste0("f", idx[2, i], "_fitness"),
                    color = "color_type")) +
                ggplot2::labs(color = "")
            if (parlist[["no_folded_states"]] == 1) {
                plot_list[[pidx]] <- plot_list[[pidx]] +
                    ggplot2::geom_line(ggplot2::aes_string(
                        x = paste0("f", idx[1, i], "_pred"),
                        y = paste0("f", idx[2, i], "_pred")), color = "black")
            } else {
                plot_list[[pidx]] <- plot_list[[pidx]] +
                    ggplot2::geom_point(ggplot2::aes_string(
                        x = paste0("f", idx[1, i], "_pred"),
                        y = paste0("f", idx[2, i], "_pred")), color = "black", alpha = 0.3) +
                    ggplot2::geom_smooth(ggplot2::aes_string(
                        x = paste0("f", idx[1, i], "_pred"),
                        y = paste0("f", idx[2, i], "_pred")), color = "black")
            }
            if (is.factor(vd_plot$color_type)) {
                plot_list[[pidx]] <- plot_list[[pidx]] +
                    ggplot2::scale_color_brewer(palette = "Set1")
            } else {
                plot_list[[pidx]] <- plot_list[[pidx]] +
                ggplot2::scale_color_gradient(low = col_orange, high = col_purple)
            }
            if (pidx == 1) {
                plot_list[[pidx]] <- plot_list[[pidx]] +
                    ggplot2::theme(legend.position = c(0.75, 0.25))
            } else {
                plot_list[[pidx]] <- plot_list[[pidx]] +
                    ggplot2::theme(legend.position = "none")
            }
            pidx <- pidx + 1
        }
    }

    # abundance versus binding fitness datasets
    for (fi in 1:varlist[["no_abd_datasets"]]) {
        for (bi in 1:varlist[["no_bind_datasets"]]) {
            plot_list[[pidx]] <- ggplot2::ggplot(vd[
                    !is.na(get(paste0("b", bi, "_fitness"))) &
                    !is.na(get(paste0("f", fi, "_fitness"))) &
                    !is.na(color_type)]) +
                ggplot2::geom_point(ggplot2::aes_string(
                    x = paste0("b", bi, "_fitness"),
                    y = paste0("f", fi, "_fitness"),
                    color = "color_type")) +
                ggplot2::theme(legend.position = "none")
            if (parlist[["no_folded_states"]] == 1) {
                plot_list[[pidx]] <- plot_list[[pidx]] +
                    ggplot2::geom_line(ggplot2::aes_string(
                        x = paste0("b", bi, "_pred_bddg0"),
                        y = paste0("f", fi, "_pred")), color = "black")
            } else {
                plot_list[[pidx]] <- plot_list[[pidx]] +
                    ggplot2::geom_line(ggplot2::aes_string(
                        x = paste0("b", bi, "_pred_fAasfB"),
                        y = paste0("f", fi, "_pred_fAasfB")), color = "black") +
                    ggplot2::geom_line(ggplot2::aes_string(
                        x = paste0("b", bi, "_pred_fBddg0"),
                        y = paste0("f", fi, "_pred_fBddg0")), color = "black", lty = 2) +
                    ggplot2::geom_line(ggplot2::aes_string(
                        x = paste0("b", bi, "_pred_fAddg0"),
                        y = paste0("f", fi, "_pred_fAddg0")), color = "black", lty = 3)

            }
            if (is.factor(vd_plot$color_type)) {
                plot_list[[pidx]] <- plot_list[[pidx]] + ggplot2::scale_color_brewer(palette = "Set1")
            } else {
                plot_list[[pidx]] <- plot_list[[pidx]] + ggplot2::scale_color_gradient(low = col_orange, high = col_purple)
            }

            if (fi == datasets_ab[1] & bi == datasets_ab[2]) {
                plot_ind <- pidx
            }
            pidx <- pidx + 1
        }
    }

    # binding fitness datasets against one another
    if (varlist[["no_bind_datasets"]] > 1) {
        idx <- utils::combn(varlist[["no_bind_datasets"]],2)

        for (i in 1:ncol(idx)) {
            plot_list[[pidx]] <- ggplot2::ggplot(vd[
                    !is.na(get(paste0("b", idx[1, i], "_fitness"))) &
                    !is.na(get(paste0("b", idx[2, i], "_fitness"))) &
                    !is.na(color_type)]) +
                ggplot2::geom_point(ggplot2::aes_string(
                    x = paste0("b", idx[1, i], "_fitness"),
                    y = paste0("b", idx[2, i], "_fitness"),
                    color = "color_type")) +
                ggplot2::theme(legend.position = "none") +
                ggplot2::geom_point(ggplot2::aes_string(
                    x = paste0("b", idx[1, i], "_pred"),
                    y = paste0("b", idx[2, i], "_pred")), color = "black", alpha = 0.3) +
                ggplot2::geom_smooth(ggplot2::aes_string(
                    x = paste0("b", idx[1, i], "_pred"),
                    y = paste0("b", idx[2, i], "_pred")), color = "black")
            if (is.factor(vd_plot$color_type)) {
                plot_list[[pidx]] <- plot_list[[pidx]] + ggplot2::scale_color_brewer(palette = "Set1")
            } else {
                plot_list[[pidx]] <- plot_list[[pidx]] + ggplot2::scale_color_gradient(low = col_orange, high = col_purple)
            }
            pidx <- pidx + 1
        }
    }

    # plot comparison between desired datasets
    ggplot2::ggsave(plot_list[[plot_ind]],
        file = file.path(dataset_folder, model_name, paste0("results/fitness_scatter_dg_relationship_fA", datasets_ab[1], "_fB", datasets_ab[2], "_", color_type, ".pdf")),
        width = 4,
        height = 3.5)

    # plot all comparisions
    if (varlist[["no_abd_datasets"]] > 1 & varlist[["no_bind_datasets"]] > 1) {
        p <- gridExtra::grid.arrange(grobs = plot_list,
            nrow = ceiling(sqrt(length(plot_list))),
            ncol = ceiling(sqrt(length(plot_list))))
        ggplot2::ggsave(p,
            file = file.path(dataset_folder, model_name, paste0("results/fitness_scatter_dg_relationship_all_", color_type, ".pdf")),
            width = 4 * ceiling(sqrt(length(plot_list))),
            height = 3.5 * ceiling(sqrt(length(plot_list))))
    }

}
