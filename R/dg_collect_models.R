#' Reads and collects individual dg models from tmp folder
#'
#' @param dataset_folder absolute path to the dataset folder, is created if non-existent
#' @param model_name name of the model that should be computed on the dataset
#' @param model_averaging c("median" [default], "mean"), how to average estimated parameters over different train/test sets
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
    model_averaging = "median",
    which_test_set = 0,
    stage = "model"
){

    ggplot2::theme_set(ggplot2::theme_bw(base_size = 8))

    if (stage == "model") {

        ## create results folder
        dir.create(file.path(dataset_folder, model_name, "results"))


        ## load dataset and parameters
        # load varlist
        load(file = file.path(dataset_folder, "data/fitness_dataset.RData"))
        # load parlist
        load(file = file.path(dataset_folder, model_name, "parameter_list.RData"))



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
            dt_models <- rbind(dt_models, X)
          }
        }


        ## compare global parameters across all models
        global_pars <- dt_models[,
            .SD,
            .SDcols = names(dt_models)[!grepl("ddg", names(dt_models))]]
        global_pars[, objective := log10(objective)]

        p <- GGally::ggpairs(global_pars[objective < (min(objective) + 0.5)],
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
        global_pars <- best_models[,
            .SD,
            .SDcols = names(best_models)[!grepl("ddg", names(best_models))]]
        global_pars[, objective := log10(objective)]

        p <- GGally::ggpairs(global_pars,
            columns = grep("^[obf]",names(global_pars)))
        ggplot2::ggsave(p,
            file = file.path(dataset_folder, model_name, "results", "global_par_distribution_bestmodels.pdf"),
            width = 15,
            height = 15)


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
                    variant_data[, .(
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
                    variant_data[, .(
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
        p <- gridExtra::grid.arrange(grobs = plot_list,
            nrow = ceiling(sqrt(length(plot_list))),
            ncol = ceiling(sqrt(length(plot_list))))
        ggplot2::ggsave(p,
            file = file.path(dataset_folder, model_name, "results", paste0("prediction_performance_fitness.pdf")),
            width = 6 * ceiling(sqrt(length(plot_list))),
            height = 5 * ceiling(sqrt(length(plot_list))))


        ## plot Pearson's R across train/test sets
        pp_melt <- reshape2::melt(prediction_performance, id.vars = "test_set")
        p <- ggplot2::ggplot(pp_melt, ggplot2::aes(x = test_set, y = value, fill = variable)) +
            ggplot2::geom_bar(stat = "identity", position="dodge") +
            ggplot2::scale_x_continuous(breaks = sort(unique(best_models[, test_set]))) +
            ggplot2::labs(x = "train/test set",
                    y = "Pearson's R",
                    title = paste0("average R: ", pp_melt[, .(value = mean(value)), variable][, paste0(variable, " = ", round(value, 3), collapse = ", ")]))
        ggplot2::ggsave(p,
                file = file.path(dataset_folder, model_name, "results", paste0("prediction_performance_fitness_R.pdf")),
                width = 5,
                height = 3)


        ## average over parameters from different train/test sets
        if (model_averaging == "median") {
            avg_model <- best_models[, lapply(.SD,stats::median), .SDcols = !grepl("^[toci]", names(best_models))]
        } else if (model_averaging == "mean") {
            avg_model <- best_models[, lapply(.SD, mean), .SDcols = !grepl("^[toci]", names(best_models))]
        } else {
            print("model_averaging parameter not 'median' or 'mean'")
        }

        # long table format
        avg_model <- data.table(parameter = names(avg_model), value = avg_model[, unlist(.SD)])


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
            dt_models <- X
          } else {
            dt_models <- rbind(dt_models, X)
          }
        }


        ## calculate mean and sd over bootstraps
        boot_mean <- dt_models[, lapply(.SD,mean), .SDcols = !grepl("^[toci]", names(dt_models))]
        boot_sd <- dt_models[, lapply(.SD,stats::sd), .SDcols = !grepl("^[toci]", names(dt_models))]
        # long table format
        avg_boot_model <- merge(data.table(parameter = names(boot_mean),
                                            boot_mean = boot_mean[, unlist(.SD)]),
                                data.table(parameter = names(boot_sd),
                                            boot_sd = boot_sd[, unlist(.SD)]), all = T)

        ## add to avg_model table
        model_results[["avg_model"]] <- merge(model_results[["avg_model"]],
                                              avg_boot_model,
                                              by = "parameter",
                                              all = T)


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
