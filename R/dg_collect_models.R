#' Reads and collects individual dg models from tmp folder
#'
#' @param dataset_folder absolute path to the dataset folder, is created if non-existent
#' @param model_name name of the model that should be computed on the dataset
#' @param model_averaging c("mean" [default], "median"), how to average estimated parameters over different train/test sets
#' @param which_test_set for 0, the default, reads individual models from training/test sets, for a number greater 0, only reads models from this training/test set
#' @param stage c("model" [default], "bootstrap", "pll"), which stage to collect models from, for bootstrap and pll stage, the functions adds their results to the dg_model.txt table
#'
#' @return returns a data.table with fitted values for all model parameters
#' @import data.table
#'
#' @export
#'

dg_collect_models <- function(
    dataset_folder,
    model_name,
    model_averaging = "mean",
    which_test_set = 0,
    stage = "model"
){

    objective <- test_set <- varlist <- parlist <- x <- y <- value <- variable <- parameter <- convergence <- type <- dataset <- dgwt_type <- fit_type <- iteration <- NULL

    ggplot2::theme_set(ggplot2::theme_bw(base_size = 8))

    ## load dataset and parameters
    # load varlist
    load(file = file.path(dataset_folder, "data/fitness_dataset.RData"))
    # load parlist
    load(file = file.path(dataset_folder, model_name, "parameter_list.RData"))

    if (stage == "model") {

        ## create results folder
        dir.create(file.path(dataset_folder, model_name, "results"))

        ## collect all models from tmp folder
        model_files <- list.files(file.path(dataset_folder, model_name, "tmp"))
        if (which_test_set == 0) {
            model_files <- model_files[grep("^dg_model", model_files)]
        } else {
            model_files <- model_files[grep(paste0("^dg_model_testset", which_test_set, "_"), model_files)]
        }
        for (i in 1:length(model_files)) {
          X = fread(file.path(dataset_folder, model_name, "tmp", model_files[i]))
          if (i == 1){
            dt_models <- X
          } else {
            if (nrow(X) > 0) {
                dt_models <- rbind(dt_models, X)
            } else {
                print(paste0(model_files[i], " empty!"))
            }
          }
        }

        print(paste0("collected ", nrow(dt_models), " models"))
        print(paste0(dt_models[convergence == 0,.N], "/", nrow(dt_models), " models converged"))

        ## compare global parameters across all models
        global_pars <- dt_models[,
            .SD,
            .SDcols = names(dt_models)[!grepl("ddg", names(dt_models))]]
        global_pars[, objective := log10(objective)]

        p <- GGally::ggpairs(global_pars[objective < min(c(min(objective) + 0.2, stats::quantile(objective, 0.9)))],
            columns = grep("^[obf]",names(global_pars)),
            ggplot2::aes(colour = factor(test_set)))
        ggplot2::ggsave(p,
            file = file.path(dataset_folder, model_name, "results", "global_par_distribution_allmodels.pdf"),
            width = 15,
            height = 15)



        ## extract best model per train/test set
        setorder(dt_models, objective)
        best_models <- dt_models[, .SD[1, ], test_set]


        ## compare global parameters across best models
        if (nrow(best_models) > 1) {

            global_pars <- reshape2::melt(best_models[,
                .SD,
                .SDcols = names(best_models)[grepl("^[fbt]", names(best_models))]], id.vars = "test_set")

            global_pars[grep("dgwt", variable), type := "dgwt"]
            global_pars[grep("fit", variable), type := "fitness"]
            global_pars[, type := factor(type, levels = c("dgwt", "fitness"))]
            levels(global_pars$type) = c("dG of wild-type state", "fitness parameters of DMS dataset")

            global_pars[, dataset := strsplit(as.character(variable), "_")[[1]][1], variable]
            global_pars[grep("^bf", dataset), dataset := gsub("[fAB]", "", dataset)]
            global_pars[grep("^f", dataset), dataset := gsub("[AB]", "", dataset)]

            p <- ggplot2::ggplot(global_pars, ggplot2::aes(x = variable, y = value, group = variable, color = dataset)) +
                ggplot2::geom_boxplot() +
                ggplot2::geom_point() +
                ggplot2::facet_wrap(type ~ ., scales = "free") +
                ggplot2::expand_limits(y = 0) +
                ggplot2::labs(x = "", y = "estimates") +
                ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1))

            ggplot2::ggsave(p,
                    file = file.path(dataset_folder, model_name, "results/global_par_distribution_bestmodels.pdf"),
                    width = 10,
                    height = 4)
        }



        ## predictive performance per train/test set
        X <- dg__fitness_from_model(
            dg_model0 = best_models,
            variant_data0 = copy(varlist[["variant_data"]]),
            varlist = varlist,
            parlist = parlist,
            calc_performance = TRUE,
            per_testset = TRUE
        )
        variant_data <- X[["variant_data"]]
        prediction_performance <- X[["prediction_performance"]]


        ## plot predited versus measured fitness
        plot_list = list()
        for (ds in 1:varlist[["no_abd_datasets"]]) {
            plot_list[[ds]] <- ggplot2::ggplot(
                    variant_data[, list(
                            x = unlist(.SD[, 1]),
                            y = unlist(.SD[, 2]),
                            .SD[,3]),
                        .SDcols = c(paste0("f", ds, "_pred"),
                            paste0("f", ds, "_fitness"),
                            "test_set")],
                    ggplot2::aes(x,y)) +
                ggplot2::geom_point() +
                ggplot2::geom_abline(color = "red") +
                ggplot2::geom_smooth() +
                ggplot2::scale_x_continuous(breaks = seq(0, 1, 0.1)) +
                ggplot2::scale_y_continuous(breaks = seq(0, 1, 0.1)) +
                ggpubr::stat_cor(ggplot2::aes(color = factor(test_set))) +
                ggplot2::labs(x = "predicted fitness",
                    y = "measured fitness",
                    title = paste0("abundance dataset ", ds),
                    color = "train/test set")
        }
        for (ds in 1:varlist[["no_bind_datasets"]]) {
            plot_list[[varlist[["no_abd_datasets"]] + ds]] <-
                ggplot2::ggplot(
                    variant_data[, list(
                                x = unlist(.SD[, 1]),
                                y = unlist(.SD[, 2]),
                                .SD[,3]),
                            .SDcols = c(paste0("b", ds, "_pred"),
                                paste0("b", ds, "_fitness"),
                                "test_set")],
                        ggplot2::aes(x,y)) +
                ggplot2::geom_point() +
                ggplot2::geom_abline(color = "red") +
                ggplot2::geom_smooth() +
                ggplot2::scale_x_continuous(breaks = seq(0, 1, 0.1)) +
                ggplot2::scale_y_continuous(breaks = seq(0, 1, 0.1)) +
                ggpubr::stat_cor(ggplot2::aes(color = factor(test_set))) +
                ggplot2::labs(
                    x = "predicted fitness",
                    y = "measured fitness",
                    title = paste0("binding dataset ", ds),
                    color = "train/test set")
        }
        # p <-
        ggplot2::ggsave(gridExtra::grid.arrange(grobs = plot_list,
                          nrow = ceiling(sqrt(length(plot_list))),
                          ncol = ceiling(sqrt(length(plot_list)))),
            file = file.path(dataset_folder, model_name, "results/prediction_performance_fitness.pdf"),
            width = 6 * ceiling(sqrt(length(plot_list))),
            height = 5 * ceiling(sqrt(length(plot_list))))


        ## plot Pearson's R across train/test sets
        pp_melt <- reshape2::melt(prediction_performance, id.vars = "test_set")
        p <- ggplot2::ggplot(pp_melt, ggplot2::aes(x = test_set, y = value, fill = variable)) +
            ggplot2::geom_bar(stat = "identity", position="dodge") +
            ggplot2::scale_x_continuous(breaks = sort(unique(best_models[, test_set]))) +
            ggplot2::labs(x = "train/test set",
                    y = "Pearson's R",
                    title = paste0("average R: ", pp_melt[, list(value = mean(value)), variable][, paste0(variable, " = ", round(value, 3), collapse = ", ")]))
        ggplot2::ggsave(p,
                file = file.path(dataset_folder, model_name, "results/prediction_performance_fitness_R.pdf"),
                width = 5,
                height = 3)


        ## average over parameters from different train/test sets
        best_models_melt <- reshape2::melt(best_models[, .SD, .SDcols = !grepl("^[toci]", names(best_models))])
        parameter_vec <- unique(best_models_melt$variable)
        if (model_averaging == "median") {
            avg_model <- best_models_melt[value != 0, list(value = stats::median(value)), list(parameter = variable)]
        } else if (model_averaging == "mean") {
            avg_model <- best_models_melt[value != 0, list(value = mean(value)), list(parameter = variable)]
        } else {
            print("model_averaging parameter not 'median' or 'mean'")
        }
        data.table::setkey(avg_model, parameter)
        avg_model <- avg_model[list(parameter_vec)]
        avg_model[is.na(value), value := 0]

        ## predict fitness with average model parameters
        X <- dg__fitness_from_model(
            dg_model0 = avg_model,
            variant_data0 = copy(varlist[["variant_data"]]),
            varlist = varlist,
            parlist = parlist,
            calc_performance = TRUE
        )
        # output Pearson R
        pp <- melt(X[["prediction_performance"]][, .SD, .SDcols = grep("^[fb]", names(X[["prediction_performance"]]))])
        cat("prediction performance of average model: \n", pp[, paste0(variable, " = ", round(value,3), "\n")])

        # predict all fitness values given average parameters
        variant_data <- X[["variant_data"]]

        ## write models/parameter to RData file
        model_results = list(best_models = best_models,
                avg_model = avg_model,
                variant_data = variant_data)

        save(model_results,
          file = file.path(dataset_folder, model_name, "model_results.RData"),
          quote = F,
          row.names = F)

        ## plot parameter/fitness relationships


    } else if (stage == "bootstrap") {

        ## load initial models
        load(file = file.path(dataset_folder, model_name, "model_results.RData"))


        ## collect bootstrapped models
        model_files <- list.files(file.path(dataset_folder, model_name, "tmp"))
        model_files <- model_files[grep("^dg_bootstrap", model_files)]
        for (i in 1:length(model_files)) {
          X = fread(file.path(dataset_folder, model_name, "tmp", model_files[i]))
          if (i == 1){
            dt_boot <- X
          } else {
            dt_boot <- rbind(dt_boot, X)
          }
        }

        print(paste0("collected ", nrow(dt_boot), " bootstrapped models (", dt_boot[convergence == 0,.N], " converged)"))


        ## compare objective of bootstrapped models to avg models
        # for this, run the avg model quickly through the optimizer
        # get parameters from average dG model
        start_par <- model_results[["avg_model"]][, value]
        names(start_par) <- model_results[["avg_model"]][, parameter]
        varlist[["mutxvar"]] <- Matrix::t(varlist[["varxmut"]])
        # call optimzer
        dt_avg <- dg__model_optim(
            start_par = start_par,
            parlist = parlist,
            varlist = varlist,
            maxit = 1
        )
        # convert output to data.table
        dt_avg <- data.table::data.table(
            t(dt_avg[["par"]]),
            dt_avg[["value"]],
            dt_avg[["convergence"]],
            0
        )
        names(dt_avg) <- c(
            parlist[["dt_par"]][, parameter],
            "objective",
            "convergence",
            "iteration"
        )


        ## calculate mean and sd over bootstraps
        boot_mean <- dt_boot[convergence == 0 & iteration != 0, lapply(.SD,mean), .SDcols = !grepl("^[toci]", names(dt_boot))]
        boot_sd <- dt_boot[convergence == 0 & iteration != 0, lapply(.SD,stats::sd), .SDcols = !grepl("^[toci]", names(dt_boot))]
        # long table format
        avg_boot_model <- merge(data.table(parameter = names(boot_mean),
                                            boot_mean = boot_mean[, unlist(.SD)]),
                                data.table(parameter = names(boot_sd),
                                            boot_sd = boot_sd[, unlist(.SD)]), all = T)

        ## add to avg_model table
        model_results[["avg_model"]] <- merge(model_results[["avg_model"]][,
                                                    .SD,
                                                    .SDcols = c("parameter",
                                                                "value",
                                                                grep("pll", names(model_results[["avg_model"]])))],
                                              avg_boot_model,
                                              by = "parameter",
                                              all = T)


        ## plot ddg parameter values
        ddg <- model_results[["avg_model"]][grepl("ddg", parameter) & value != 0]
        ddg[, type := paste0(strsplit(parameter, "_")[[1]][2:3], collapse = "_"), parameter]
        if (parlist[["no_folded_states"]] == 1) {
            ddg[, type := factor(type, levels = c("f_ddg", "b_ddg"))]
            levels(ddg$type) = c("ddG of folding", "ddG of binding")
        } else {
            ddg[, type := factor(type, levels = c("fA_ddg", "fB_ddg", "b_ddg"))]
            levels(ddg$type) = c("ddG of folding state A", "ddG of folding state B", "ddG of binding")
        }

        # plot original versus bootstrapped (mean) values
        p1 <- ggplot2::ggplot(ddg, ggplot2::aes(value, boot_mean)) +
            ggplot2::geom_point(alpha = 0.1) +
            ggplot2::facet_wrap(type ~ .) +
            ggplot2::expand_limits(y = 0) +
            ggplot2::geom_hline(yintercept = 0, linetype = 2) +
            ggplot2::geom_vline(xintercept = 0, linetype = 2) +
            ggplot2::geom_abline(linetype = 3, color = "red") +
            ggplot2::labs(x = "initial estimate", y = "bootstrapped mean")

        # plot bootstrapped mean versus sd
        p2 <- ggplot2::ggplot(ddg, ggplot2::aes(boot_mean, boot_sd)) +
            ggplot2::geom_point(alpha = 0.1) +
            ggplot2::facet_wrap(type ~ .) +
            ggplot2::expand_limits(y = 0) +
            ggplot2::geom_hline(yintercept = 0, linetype = 2) +
            ggplot2::geom_abline(linetype = 3, color = "red") +
            ggplot2::labs(x = "bootstrapped mean", y = "bootstrapped sd")

        # check how difference between original and bootstrapped (mean) values relate to bootstrapped sd
        p3 <- ggplot2::ggplot(ddg, ggplot2::aes(value - boot_mean, boot_sd)) +
            ggplot2::geom_point(alpha = 0.1) +
            ggplot2::facet_wrap(type ~ .) +
            ggplot2::expand_limits(y = 0) +
            ggplot2::geom_vline(xintercept = 0, linetype = 2) +
            ggplot2::geom_abline(linetype = 3, color = "red", slope = c(-1, 1)) +
            ggplot2::labs(x = "diff (initial - bootstrapped) estimate", y = "bootstrapped sd")

        # p <-

        ggplot2::ggsave(gridExtra::grid.arrange(grobs = list(p1, p2, p3),
                  nrow = 3,
                  ncol = 1),
                file = file.path(dataset_folder, model_name, "results/bootstrap_ddg_pars.pdf"),
                width = 4 * ddg[, length(unique(type))] ,
                height = 12)


        ## plot global parameter values
        global_pars <- model_results[["avg_model"]][!grepl("ddg", parameter)]

        # add in log10(objective)
        global_pars <- rbind(global_pars,
            data.table(parameter = "objective",
                value = log10(dt_avg$objective),
                boot_mean = mean(log10(dt_boot$objective)),
                boot_sd = stats::sd(log10(dt_boot$objective))))

        global_pars[grep("dgwt", parameter), type := "dgwt"]
        global_pars[grep("fit", parameter), type := "fitness"]
        global_pars[parameter == "objective", type := "objective"]

        # global_pars[, type := factor(type, levels = c("dgwt", "fitness","objective"))]
        # levels(global_pars$type) = c("dG of wild-type state", "fitness parameters of DMS dataset", "log10(objective func)")

        global_pars[, dataset := strsplit(parameter, "_")[[1]][1], parameter]
        global_pars[grep("^bf", dataset), dataset := gsub("[fAB]", "", dataset)]
        global_pars[grep("^f", dataset), dataset := gsub("[AB]", "", dataset)]

        global_pars[type == "dgwt" & grepl("^b?f", parameter), dgwt_type := "folding"]
        global_pars[type == "dgwt" & grepl("^b[^f]", parameter), dgwt_type := "binding"]
        p1 <- ggplot2::ggplot(global_pars[type == "dgwt"],
                ggplot2::aes(value,
                    boot_mean,
                    color = dataset,
                    shape = dgwt_type)) +
            ggplot2::geom_point() +
            ggplot2::geom_pointrange(ggplot2::aes(ymin = boot_mean - boot_sd, ymax = boot_mean + boot_sd)) +
            ggplot2::expand_limits(x = 0, y = 0) +
            ggplot2::geom_abline(linetype = 2) +
            ggplot2::labs(title = "wild-type dG",
                x = "initial estimate",
                y = "bootstrapped estimate (+-sd)")

        global_pars[type == "fitness" & grepl("fit0", parameter), fit_type := "background"]
        global_pars[type == "fitness" & grepl("fitwt", parameter), fit_type := "wild-type"]
        p2 <- ggplot2::ggplot(global_pars[type == "fitness"],
                ggplot2::aes(value,
                    boot_mean,
                    color = dataset,
                    shape = fit_type)) +
            ggplot2::geom_point() +
            ggplot2::geom_pointrange(ggplot2::aes(ymin = boot_mean - boot_sd, ymax = boot_mean + boot_sd)) +
            ggplot2::expand_limits(x = c(0,1), y = c(0,1)) +
            ggplot2::scale_x_continuous(breaks = seq(0,1,0.25)) +
            ggplot2::scale_y_continuous(breaks = seq(0,1,0.25)) +
            ggplot2::geom_abline(linetype = 2) +
            ggplot2::labs(title = "fitness parameters DMS dataset",
                x = "initial estimate",
                y = "bootstrapped estimate (+-sd)")

        p3 <- ggplot2::ggplot(global_pars[type == "objective"],
                ggplot2::aes(value,
                    boot_mean)) +
            ggplot2::geom_point() +
            ggplot2::geom_pointrange(ggplot2::aes(ymin = boot_mean - boot_sd, ymax = boot_mean + boot_sd)) +
            ggplot2::scale_x_continuous(
                limits = global_pars[type == "objective", c(value - 0.5, value + 0.5)],
                breaks = seq(0,10,0.25)) +
            ggplot2::scale_y_continuous(
                limits = global_pars[type == "objective", c(boot_mean - 0.5, boot_mean + 0.5)],
                breaks = seq(0,10,0.25)) +
            ggplot2::geom_abline(linetype = 2) +
            ggplot2::labs(title = "value objective function",
                x = "initial estimate",
                y = "bootstrapped estimate (+-sd)")

        # p <-
        ggplot2::ggsave(gridExtra::grid.arrange(grobs = list(p1, p2, p3),
                  nrow = 1,
                  ncol = 3),
                file = file.path(dataset_folder, model_name, "results/bootstrap_global_pars.pdf"),
                width = 13,
                height = 4)


        ## write to file (append to previous table)
        save(model_results,
          file = file.path(dataset_folder, model_name, "model_results.RData"),
          quote = F,
          row.names = F)


    } else if (stage == "pll") {
        ## collect pll runs


        ## plot some basic stuff


        ## write to file (append to previous table?!)



    } else {print("stage parameter not defined properly")}
}
