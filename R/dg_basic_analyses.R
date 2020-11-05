#' Produces some basic analyses plots
#'
#' @param dataset_folder absolute path to the dataset folder, is created if non-existent
#' @param model_name name of the model that should be computed on the dataset
#' @param color_type variable with which to color variants; can be Nmut or any column_name from structural_properties file; any other than the default "" will exclude variants with multiple mutations and the wild-type variant
#' @param max_variants maximual number of variants to plot, default = 3e3, will be downsampled if # variants > max_variants
#' @param datasets_ab index of datasets to use (in case multiple abundance or binding datasets were supplied)
#' @param stage specifies the stage ("model" or "bootstrap") from which to plot results, default = "" uses bootstrapped results if they are detected, initial model results otherwise
#' @param fdr_thres FDR cutoff to call ddg value significantly different from zero (if plotting bootstrapped parameter estimates), default = 0.1
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
    datasets_ab = c(1, 1),
    stage = "",
    fdr_thres = 1e-2
){

  varlist <- parlist <- f_ddg <- parameter <- value <- fA_ddg <- fB_ddg <-
    b_ddg <- aa_subs <- Pos <- b_ddg_fdr <- b_ddg_sd <- boot_mean <- boot_sd <-
    fA_ddg_fdr <- fA_ddg_sd <- fB_ddg_fdr <- fB_ddg_sd <- f_ddg_fdr <-
    f_ddg_sd <- RSA_unbound <- ddg <- ddg_sd <- ddg_type <- ddg_weight <-
    scHAmin_ligand <- structural_property <- structural_property_value <- y <-
    NULL

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

  ## stage-specific plotting
  if (stage == "") {
    if (any(names(model_results[["avg_model"]]) == "boot_mean")) {
      stage <- "bootstrap"
    } else {
      stage <- "model"
    }
  }

  print(paste0("plotting results for ", ifelse(stage == "bootstrap", "bootstrapped", "initial"), " model estimates"))

  # use bootstrapped mean to plot results
  if (stage == "bootstrap") {
    model_results[["avg_model"]][, value := boot_mean]
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
      if (stage == "bootstrap") {
        model_results[["variant_data"]][, f_ddg_sd := as.numeric(varlist[["varxmut"]] %*% model_results[["avg_model"]][grep("f_ddg", parameter), boot_sd])]
      }
  } else {
      model_results[["variant_data"]][, fA_ddg := as.numeric(varlist[["varxmut"]] %*% model_results[["avg_model"]][grep("fA_ddg", parameter), value])]
      model_results[["variant_data"]][, fB_ddg := as.numeric(varlist[["varxmut"]] %*% model_results[["avg_model"]][grep("fB_ddg", parameter), value])]
      if (stage == "bootstrap") {
        model_results[["variant_data"]][, fA_ddg_sd := as.numeric(varlist[["varxmut"]] %*% model_results[["avg_model"]][grep("fA_ddg", parameter), boot_sd])]
        model_results[["variant_data"]][, fB_ddg_sd := as.numeric(varlist[["varxmut"]] %*% model_results[["avg_model"]][grep("fB_ddg", parameter), boot_sd])]
      }
  }
  model_results[["variant_data"]][, b_ddg := as.numeric(varlist[["varxmut"]] %*% model_results[["avg_model"]][grep("b_ddg", parameter), value])]
  if (stage == "bootstrap") {
    model_results[["variant_data"]][, b_ddg_sd := as.numeric(varlist[["varxmut"]] %*% model_results[["avg_model"]][grep("b_ddg", parameter), boot_sd])]
  }

  ## define helper pars and include structural properties
  vd <- copy(model_results[["variant_data"]])

  # note position for singles
  vd[!grepl("_", aa_subs) & grepl("[0-9]", aa_subs), Pos := as.integer(paste0(strsplit(aa_subs,"")[[1]][1:(nchar(aa_subs)-1)], collapse = "")), aa_subs]

  if (color_type == "") {
      vd[, color_type := factor("all vars")]
      color_type = "allvars"
  } else if (color_type == "Nmut") {
      vd[grepl("[0-9]", aa_subs), color_type := factor(length(grep("_", aa_subs)) + 1), aa_subs]
  } else if (exists("structural_properties") == TRUE) {

      if (length(structural_properties[, unique(unlist(.SD)), .SDcols = color_type]) < 10) { #assume it's a factor
          vd <- merge(vd,
              structural_properties[, list(Pos, color_type = factor(unlist(.SD[,1]))), .SDcols = color_type],
              by = "Pos")
      } else {
          vd <- merge(vd,
              structural_properties[, list(Pos, color_type = unlist(.SD[,1])), .SDcols = color_type],
              by = "Pos")
      }
  } else {
      print("error: no structural properties file provided to color variants")
  }
  am <- model_results[["avg_model"]]


  ##############################################################################
  ## plot inter-relationship between estimated ddg values and fitness values
  ##############################################################################

  # plot f_ddg distribution and relationship with f_fitness, b_ddg and b_fitness
  vd_plot <- vd[!is.na(get(paste0("f", parlist[["str_abd"]][datasets_ab[1]], "_fitness"))) &
                !is.na(get(paste0("b", parlist[["str_bind"]][datasets_ab[2]], "_fitness"))) &
                !is.na(color_type)][sample(.N, min(c(5e3, .N)))]

  if (parlist[["no_folded_states"]] == 1) {
    if (parlist[["fix_f_dgwt"]] == FALSE) {
      f_dgwt_helper <- am[grep(paste0("^f", parlist[["str_abd"]][datasets_ab[1]], "_dgwt"), parameter), value]
      bf_dgwt_helper <- am[grep(paste0("^bf", parlist[["str_bind"]][datasets_ab[2]], "_dgwt"), parameter), value]
    } else {
      f_dgwt_helper <- bf_dgwt_helper <- am[grep("^f_dgwt", parameter), value]
    }
    if (parlist[["fix_b_dgwt"]] == FALSE) {
      b_dgwt_helper <- am[grep(paste0("^b", parlist[["str_bind"]][datasets_ab[2]], "_dgwt"), parameter), value]
    } else {
      b_dgwt_helper <- am[grep("^b_dgwt", parameter), value]
    }

    df_dist <- ggplot2::ggplot(data = vd_plot)
    if (is.factor(vd_plot$color_type)) {
        df_dist <- df_dist +
            ggplot2::geom_density(ggplot2::aes(x = f_ddg + f_dgwt_helper,
                                        color = color_type)) +
            ggplot2::scale_color_brewer(palette = "Set1")
    } else {
        df_dist <- df_dist +
            ggplot2::geom_density(ggplot2::aes(x = f_ddg + f_dgwt_helper), color = 'black')
    }
   df_dist <- df_dist +
        ggplot2::geom_vline(xintercept = f_dgwt_helper,linetype = 2) +
        ggplot2::theme(legend.position = c(0.75, 0.75)) +
        ggplot2::labs(x = "dG folding", color = "")

    df_ffitness <- ggplot2::ggplot(data = vd_plot,
            ggplot2::aes(x = f_ddg + f_dgwt_helper)) +
        ggplot2::geom_point(ggplot2::aes_string(y = paste0("f", parlist[["str_abd"]][datasets_ab[1]], "_fitness"), color = "color_type")) +
        ggplot2::geom_density2d(ggplot2::aes_string(y = paste0("f", parlist[["str_abd"]][datasets_ab[1]], "_fitness")), color = 'black', alpha = 0.5) +
        ggplot2::geom_line(ggplot2::aes_string(y = paste0("f", parlist[["str_abd"]][datasets_ab[1]], "_pred")), color = "black") +
        ggplot2::geom_hline(yintercept = am[grep(paste0("^f", parlist[["str_abd"]][datasets_ab[1]], "_fit"), parameter), c(value)],linetype = 2) +
        ggplot2::geom_vline(xintercept = f_dgwt_helper,linetype = 2) +
        ggplot2::labs(x = "dG folding",
            y = paste0("folding fitness ", parlist[["str_abd"]][datasets_ab[1]]),
            color = "") +
        ggplot2::theme(legend.position = "none")
    if (is.factor(vd_plot$color_type)) {
        df_ffitness <- df_ffitness + ggplot2::scale_color_brewer(palette = "Set1")
    } else {
        df_ffitness <- df_ffitness + ggplot2::scale_color_gradient(low = col_orange, high = col_purple)
    }


    df_bfitness <- ggplot2::ggplot(data = vd_plot,
            ggplot2::aes(x = f_ddg + bf_dgwt_helper)) +
        ggplot2::geom_point(ggplot2::aes_string(y = paste0("b", parlist[["str_bind"]][datasets_ab[2]], "_fitness"), color = "color_type")) +
        ggplot2::geom_density2d(ggplot2::aes_string(y = paste0("b", parlist[["str_bind"]][datasets_ab[2]], "_fitness")), color = 'black', alpha = 0.5) +
        ggplot2::geom_line(ggplot2::aes_string(y = paste0("b", parlist[["str_bind"]][datasets_ab[2]], "_pred_bddg0")), color = "black") +
        ggplot2::geom_hline(yintercept = am[grep(paste0("^b", parlist[["str_bind"]][datasets_ab[2]], "_fit"), parameter), c(value)],linetype = 2) +
        ggplot2::geom_vline(xintercept = bf_dgwt_helper,linetype = 2) +
        ggplot2::labs(x = "dG folding",
            y = paste0("binding fitness ", parlist[["str_bind"]][datasets_ab[2]]),
            color = "") +
        ggplot2::theme(legend.position = "none")
    if (is.factor(vd_plot$color_type)) {
        df_bfitness <- df_bfitness + ggplot2::scale_color_brewer(palette = "Set1")
    } else {
        df_bfitness <- df_bfitness + ggplot2::scale_color_gradient(low = col_orange, high = col_purple)
    }

    df_db <- ggplot2::ggplot(vd_plot[!is.na(f_ddg) & !is.na(b_ddg)],
          ggplot2::aes(x = f_ddg + f_dgwt_helper,
              y = b_ddg + b_dgwt_helper, color = color_type)) +
        ggplot2::geom_point(alpha = 0.4) +
        ggplot2::geom_density2d() +
        ggplot2::geom_vline(xintercept = f_dgwt_helper,linetype = 2) +
        ggplot2::geom_hline(yintercept = b_dgwt_helper,linetype = 2) +
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
    if (parlist[["fix_f_dgwt"]] == FALSE) {
      fA_dgwt_helper <- am[grep(paste0("^fA", parlist[["str_abd"]][datasets_ab[1]], "_dgwt"), parameter), value]
      fB_dgwt_helper <- am[grep(paste0("^fB", parlist[["str_abd"]][datasets_ab[1]], "_dgwt"), parameter), value]
      bfA_dgwt_helper <- am[grep(paste0("^bfA", parlist[["str_bind"]][datasets_ab[2]], "_dgwt"), parameter), value]
      bfB_dgwt_helper <- am[grep(paste0("^bfB", parlist[["str_bind"]][datasets_ab[2]], "_dgwt"), parameter), value]
    } else {
      fA_dgwt_helper <- bfA_dgwt_helper <- am[grep("^fA_dgwt", parameter), value]
      fB_dgwt_helper <- bfB_dgwt_helper <- am[grep("^fB_dgwt", parameter), value]
    }
    if (parlist[["fix_b_dgwt"]] == FALSE) {
      b_dgwt_helper <- am[grep(paste0("^b", parlist[["str_bind"]][datasets_ab[2]], "_dgwt"), parameter), value]
    } else {
      b_dgwt_helper <- am[grep("^b_dgwt", parameter), value]
    }

    dfA_dist <- ggplot2::ggplot(data = vd_plot)
    if (is.factor(vd_plot$color_type)) {
        dfA_dist <- dfA_dist +
            ggplot2::geom_density(
                ggplot2::aes(x = fA_ddg + fA_dgwt_helper,
                    color = color_type)) +
            ggplot2::scale_color_brewer(palette = "Set1")
    } else {
        dfA_dist <- dfA_dist +
            ggplot2::geom_density(
                ggplot2::aes(x = fA_ddg + fA_dgwt_helper),
                    color = 'black')
    }
    dfA_dist <- dfA_dist +
        ggplot2::geom_vline(xintercept = fA_dgwt_helper, linetype = 2) +
        ggplot2::theme(legend.position = c(0.75, 0.75)) +
        ggplot2::labs(x = "dG folding state A", color = "")

    dfB_dist <- ggplot2::ggplot(data = vd_plot)
    if (is.factor(vd_plot$color_type)) {
        dfB_dist <- dfB_dist +
            ggplot2::geom_density(
                ggplot2::aes(x = fB_ddg + fB_dgwt_helper,
                    color = color_type)) +
            ggplot2::scale_color_brewer(palette = "Set1")
    } else {
        dfB_dist <- dfB_dist +
            ggplot2::geom_density(
                ggplot2::aes(x = fB_ddg + fB_dgwt_helper),
                    color = 'black')
    }
    dfB_dist <- dfB_dist +
        ggplot2::geom_vline(xintercept = fB_dgwt_helper, linetype = 2) +
        ggplot2::theme(legend.position = "none") +
        ggplot2::coord_flip() +
        ggplot2::labs(x = "dG folding state B", color = "")

    dfA_ffitness <- ggplot2::ggplot(data = vd_plot,
            ggplot2::aes(x = fA_ddg + fA_dgwt_helper)) +
        ggplot2::geom_point(ggplot2::aes_string(y = paste0("f", parlist[["str_abd"]][datasets_ab[1]], "_fitness"), color = "color_type")) +
        ggplot2::geom_density2d(ggplot2::aes_string(y = paste0("f", parlist[["str_abd"]][datasets_ab[1]], "_fitness")), color = 'black', alpha = 0.5) +
        ggplot2::geom_line(ggplot2::aes_string(y = paste0("f", parlist[["str_abd"]][datasets_ab[1]], "_pred_fBasfA")), color = "black") +
        ggplot2::geom_line(ggplot2::aes_string(y = paste0("f", parlist[["str_abd"]][datasets_ab[1]], "_pred_fBddg0")), color = "black", linetype = 3) +
        ggplot2::geom_hline(yintercept = am[grep(paste0("^f", parlist[["str_abd"]][datasets_ab[1]], "_fit"), parameter), c(value)], linetype = 2) +
        ggplot2::geom_vline(xintercept = fA_dgwt_helper, linetype = 2) +
        ggplot2::labs(x = "dG folding state A",
            y = paste0("folding fitness ", parlist[["str_abd"]][datasets_ab[1]]),
            color = "") +
        ggplot2::theme(legend.position = "none")
    if (is.factor(vd_plot$color_type)) {
        dfA_ffitness <- dfA_ffitness + ggplot2::scale_color_brewer(palette = "Set1")
    } else {
        dfA_ffitness <- dfA_ffitness + ggplot2::scale_color_gradient(low = col_orange, high = col_purple)
    }

    dfB_ffitness <- ggplot2::ggplot(data = vd_plot,
            ggplot2::aes(x = fB_ddg + fB_dgwt_helper)) +
        ggplot2::geom_point(ggplot2::aes_string(y = paste0("f", parlist[["str_abd"]][datasets_ab[1]], "_fitness"), color = "color_type")) +
        ggplot2::geom_density2d(ggplot2::aes_string(y = paste0("f", parlist[["str_abd"]][datasets_ab[1]], "_fitness")), color = 'black', alpha = 0.5) +
        ggplot2::geom_line(
            ggplot2::aes_string(y = paste0("f", parlist[["str_abd"]][datasets_ab[1]], "_pred_fAasfB")),color = "black") +
        ggplot2::geom_line(ggplot2::aes_string(y = paste0("f", parlist[["str_abd"]][datasets_ab[1]], "_pred_fAddg0")), color = "black", linetype = 3) +
        ggplot2::geom_hline(yintercept = am[grep(paste0("^f", parlist[["str_abd"]][datasets_ab[1]], "_fit"), parameter), c(value)], linetype = 2) +
        ggplot2::geom_vline(xintercept = fB_dgwt_helper, linetype = 2) +
        ggplot2::labs(x = "dG folding state B",
            y = paste0("folding fitness ", parlist[["str_abd"]][datasets_ab[1]]),
            color = "") +
        ggplot2::theme(legend.position = "none")
    if (is.factor(vd_plot$color_type)) {
        dfB_ffitness <- dfB_ffitness + ggplot2::scale_color_brewer(palette = "Set1")
    } else {
        dfB_ffitness <- dfB_ffitness + ggplot2::scale_color_gradient(low = col_orange, high = col_purple)
    }

    dfA_dfB <- ggplot2::ggplot(vd_plot,
          ggplot2::aes(x = fA_ddg + fA_dgwt_helper,
              y = fB_ddg + fB_dgwt_helper, color = color_type)) +
        ggplot2::geom_point(alpha = 0.4) +
        ggplot2::geom_density2d() +
        ggplot2::geom_vline(xintercept = fA_dgwt_helper, linetype = 2) +
        ggplot2::geom_hline(yintercept = fB_dgwt_helper, linetype = 2) +
        ggplot2::geom_abline(color = "red",
            intercept = fB_dgwt_helper -
                        fA_dgwt_helper) +
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
          ggplot2::aes(x = fA_ddg + fA_dgwt_helper,
              y = b_ddg + b_dgwt_helper, color = color_type)) +
        ggplot2::geom_point(alpha = 0.4) +
        ggplot2::geom_density2d() +
        ggplot2::geom_vline(xintercept = fA_dgwt_helper,linetype = 2) +
        ggplot2::geom_hline(yintercept = b_dgwt_helper,linetype = 2) +
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
          ggplot2::aes(x = fB_ddg + fB_dgwt_helper,
              y = b_ddg + b_dgwt_helper, color = color_type)) +
        ggplot2::geom_point(alpha = 0.4) +
        ggplot2::geom_density2d() +
        ggplot2::geom_vline(xintercept = fB_dgwt_helper,linetype = 2) +
        ggplot2::geom_hline(yintercept = b_dgwt_helper,linetype = 2) +
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
            ggplot2::aes(x = fA_ddg + bfA_dgwt_helper)) +
        ggplot2::geom_point(ggplot2::aes_string(y = paste0("b", parlist[["str_bind"]][datasets_ab[2]], "_fitness"), color = "color_type")) +
        ggplot2::geom_density2d(ggplot2::aes_string(y = paste0("b", parlist[["str_bind"]][datasets_ab[2]], "_fitness")), color = 'black', alpha = 0.5) +
        ggplot2::geom_line(ggplot2::aes_string(y = paste0("b", parlist[["str_bind"]][datasets_ab[2]], "_pred_fBasfA")),color = "black") +
        ggplot2::geom_line(ggplot2::aes_string(y = paste0("b", parlist[["str_bind"]][datasets_ab[2]], "_pred_fBddg0")), color = "black", linetype = 3) +
        ggplot2::geom_hline(yintercept = am[grep(paste0("^b", parlist[["str_bind"]][datasets_ab[2]], "_fit"), parameter), c(value)],linetype = 2) +
        ggplot2::geom_vline(xintercept = bfA_dgwt_helper,linetype = 2) +
        ggplot2::labs(x = "dG folding state A",
            y = paste0("binding fitness ", parlist[["str_bind"]][datasets_ab[2]]),
            color = "") +
        ggplot2::theme(legend.position = "none")
    if (is.factor(vd_plot$color_type)) {
        dfA_bfitness <- dfA_bfitness + ggplot2::scale_color_brewer(palette = "Set1")
    } else {
        dfA_bfitness <- dfA_bfitness + ggplot2::scale_color_gradient(low = col_orange, high = col_purple)
    }

    dfB_bfitness <- ggplot2::ggplot(data = vd_plot,
            ggplot2::aes(x = fB_ddg + bfB_dgwt_helper)) +
        ggplot2::geom_point(ggplot2::aes_string(y = paste0("b", parlist[["str_bind"]][datasets_ab[2]], "_fitness"), color = "color_type")) +
        ggplot2::geom_density2d(ggplot2::aes_string(y = paste0("b", parlist[["str_bind"]][datasets_ab[2]], "_fitness")), color = 'black', alpha = 0.5) +
        ggplot2::geom_line(
            ggplot2::aes_string(y = paste0("b", parlist[["str_bind"]][datasets_ab[2]], "_pred_fAasfB")),color = "black") +
        ggplot2::geom_line(ggplot2::aes_string(y = paste0("b", parlist[["str_bind"]][datasets_ab[2]], "_pred_fAddg0")), color = "black", linetype = 3) +
        ggplot2::geom_hline(yintercept = am[grep(paste0("^b", parlist[["str_bind"]][datasets_ab[2]], "_fit"), parameter), c(value)],linetype = 2) +
        ggplot2::geom_vline(xintercept = bfB_dgwt_helper,linetype = 2) +
        ggplot2::labs(x = "dG folding state B",
            y = paste0("binding fitness ", parlist[["str_bind"]][datasets_ab[2]]),
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
              ggplot2::aes(x = b_ddg + b_dgwt_helper,
                  color = color_type)) +
          ggplot2::scale_color_brewer(palette = "Set1")
  } else {
      db_dist <- db_dist +
          ggplot2::geom_density(data = vd_plot,
              ggplot2::aes(x = b_ddg + b_dgwt_helper), color = 'black')
  }
  db_dist <- db_dist +
      ggplot2::geom_vline(xintercept = b_dgwt_helper,linetype = 2) +
      # ggplot2::scale_x_continuous(breaks = seq(-6,6,1)) +
      ggplot2::coord_flip() +
      ggplot2::labs(x = "dG binding", color = "") +
      ggplot2::theme(legend.position = c(0.75, 0.25))


  db_bfitness <- ggplot2::ggplot(data = vd_plot,
          ggplot2::aes(x = b_ddg + b_dgwt_helper)) +
      ggplot2::geom_point(ggplot2::aes_string(y = paste0("b", parlist[["str_bind"]][datasets_ab[2]], "_fitness"), color = "color_type")) +
      ggplot2::geom_density2d(ggplot2::aes_string(y = paste0("b", parlist[["str_bind"]][datasets_ab[2]], "_fitness")), color = 'black', alpha = 0.5) +
      ggplot2::geom_line(ggplot2::aes_string(y = paste0("b", parlist[["str_bind"]][datasets_ab[2]], "_pred_fddg0")), color = "black") +
      ggplot2::geom_hline(yintercept = am[grep(paste0("^b", parlist[["str_bind"]][datasets_ab[2]], "_fit"), parameter), c(value)],linetype = 2) +
      ggplot2::geom_vline(xintercept = b_dgwt_helper,linetype = 2) +
      ggplot2::labs(x = "dG binding",
          y = paste0("binding fitness ", parlist[["str_bind"]][datasets_ab[2]]),
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
      ggplot2::ggsave(gridExtra::grid.arrange(df_dist, df_ffitness, df_bfitness,
                        df_db, db_dist, db_bfitness,
                        nrow = 2),
          file = file.path(dataset_folder, model_name,
              paste0("results/dG_fitness_f",
                parlist[["str_abd"]][datasets_ab[1]], "b", parlist[["str_bind"]][datasets_ab[2]], "_", color_type,
                ifelse(stage == "bootstrap", "_boot",""), ".pdf")),
          width = 10,
          height = 7)
  } else if (parlist[["no_folded_states"]] == 2) {
      # plot folding
      ggplot2::ggsave(gridExtra::grid.arrange(dfA_dist, dfA_ffitness, dfB_ffitness,
                        dfA_dfB, dfB_dist,
                        nrow = 2),
          file = file.path(dataset_folder, model_name,
              paste0("results/dG_folding_fitness_fA",
                parlist[["str_abd"]][datasets_ab[1]], "_fB", parlist[["str_abd"]][datasets_ab[1]], "_", color_type,
                ifelse(stage == "bootstrap", "_boot",""),  ".pdf")),
          width = 10,
          height = 7)

      # plot binding
      ggplot2::ggsave(gridExtra::grid.arrange(db_dist, dfA_db, dfB_db,
                        dfA_bfitness, dfB_bfitness, db_bfitness,
                        nrow = 2),
          file = file.path(dataset_folder, model_name,
              paste0("results/dG_binding_fitness_b",
                parlist[["str_bind"]][datasets_ab[2]], "_fA", parlist[["str_abd"]][datasets_ab[1]], "_fB", parlist[["str_abd"]][datasets_ab[1]],
                "_", color_type,
                ifelse(stage == "bootstrap", "_boot",""), ".pdf")),
          width = 10,
          height = 7)
  }



  ##############################################################################
  ## evaluate uncertainty of ddg parameters
  ##############################################################################

  if (stage == "bootstrap") {

    # select variants that match a color_type and have at least one measured folding and binding value each
    vd_plot <- vd[!is.na(color_type), x := rowSums(!is.na(.SD)) > 0, .SDcols = grep("f[0-9]*_fitness", names(vd))][
      x == TRUE, y := rowSums(!is.na(.SD)) > 0, .SDcols = grep("b[0-9]*_fitness", names(vd))][y == TRUE]

    # evaluate if parameter estimates are significantly different from zero
    if (parlist[["no_folded_states"]] == 1){
      vd_plot[, f_ddg_fdr := stats::p.adjust(2*stats::pnorm(-abs(f_ddg / f_ddg_sd)),
        method = "fdr")]
    } else {
      vd_plot[, fA_ddg_fdr := stats::p.adjust(2*stats::pnorm(-abs(fA_ddg / fA_ddg_sd)),
        method = "fdr")]
      vd_plot[, fB_ddg_fdr := stats::p.adjust(2*stats::pnorm(-abs(fB_ddg / fB_ddg_sd)),
        method = "fdr")]
    }
    vd_plot[, b_ddg_fdr := stats::p.adjust(2*stats::pnorm(-abs(b_ddg / b_ddg_sd)),
      method = "fdr")]

    if (parlist[["no_folded_states"]] == 1){
      df_db_dfsig <- ggplot2::ggplot(vd_plot,
            ggplot2::aes(x = f_ddg + f_dgwt_helper,
                y = b_ddg + b_dgwt_helper)) +
          ggplot2::geom_density2d(color = 'black', alpha = 0.25) +
          ggplot2::geom_point(ggplot2::aes(color = color_type,shape = f_ddg_fdr < fdr_thres)) +
          ggplot2::geom_vline(xintercept = f_dgwt_helper,linetype = 2) +
          # ggplot2::geom_hline(yintercept = b_dgwt_helper,linetype = 2) +
          ggplot2::scale_shape_manual(values = c(1, 19)) +
          ggplot2::labs(x = "dG folding",
              y = "dG binding",
              color = color_type,
              shape = paste0("FDR < ", fdr_thres),
              title = "sig. folding ddG")
      if (is.factor(vd_plot$color_type)) {
          df_db_dfsig <- df_db_dfsig + ggplot2::scale_color_brewer(palette = "Set1")
      } else {
          df_db_dfsig <- df_db_dfsig + ggplot2::scale_color_gradient(low = col_orange, high = col_purple)
      }

      df_db_dbsig <- ggplot2::ggplot(vd_plot,
            ggplot2::aes(x = f_ddg + f_dgwt_helper,
                y = b_ddg + b_dgwt_helper)) +
          ggplot2::geom_density2d(color = 'black', alpha = 0.25) +
          ggplot2::geom_point(ggplot2::aes(color = color_type,shape = b_ddg_fdr < fdr_thres)) +
          # ggplot2::geom_vline(xintercept = f_dgwt_helper,linetype = 2) +
          ggplot2::geom_hline(yintercept = b_dgwt_helper,linetype = 2) +
          ggplot2::scale_shape_manual(values = c(1, 19)) +
          ggplot2::labs(x = "dG folding",
              y = "dG binding",
              color = color_type,
              shape = paste0("FDR < ", fdr_thres),
              title = "sig. binding ddG")
      if (is.factor(vd_plot$color_type)) {
          df_db_dbsig <- df_db_dbsig + ggplot2::scale_color_brewer(palette = "Set1")
      } else {
          df_db_dbsig <- df_db_dbsig + ggplot2::scale_color_gradient(low = col_orange, high = col_purple)
      }

      df_db_sig <- ggplot2::ggplot(vd_plot,
            ggplot2::aes(x = f_ddg + f_dgwt_helper,
                y = b_ddg + b_dgwt_helper)) +
          ggplot2::geom_density2d(color = 'black', alpha = 0.25) +
          ggplot2::geom_point(ggplot2::aes(color = color_type,shape = f_ddg_fdr < fdr_thres & b_ddg_fdr < fdr_thres)) +
          ggplot2::geom_vline(xintercept = f_dgwt_helper,linetype = 2) +
          ggplot2::geom_hline(yintercept = b_dgwt_helper,linetype = 2) +
          ggplot2::scale_shape_manual(values = c(1, 19)) +
          ggplot2::labs(x = "dG folding",
              y = "dG binding",
              color = color_type,
              shape = paste0("FDR < ", fdr_thres),
              title = "sig. folding & binding ddGs")
      if (is.factor(vd_plot$color_type)) {
          df_db_sig <- df_db_sig + ggplot2::scale_color_brewer(palette = "Set1")
      } else {
          df_db_sig <- df_db_sig + ggplot2::scale_color_gradient(low = col_orange, high = col_purple)
      }
      # p <-
      ggplot2::ggsave(gridExtra::grid.arrange(df_db_dfsig, df_db_dbsig, df_db_sig,
                        nrow = 3,
                        ncol = 1),
          file = file.path(dataset_folder, model_name,
              paste0("results/bootstrap_dG_uncertainty_", color_type, ".pdf")),
          width = 4.5,
          height = 9)



    } else if (parlist[["no_folded_states"]] == 2) {

      # fA vs fB
      dfA_dfB_dfAsig <- ggplot2::ggplot(vd_plot,
            ggplot2::aes(x = fA_ddg + fA_dgwt_helper,
                y = fB_ddg + fB_dgwt_helper)) +
          ggplot2::geom_density2d(color = 'black', alpha = 0.25) +
          ggplot2::geom_point(ggplot2::aes(color = color_type, shape = fA_ddg_fdr < fdr_thres)) +
          ggplot2::geom_vline(xintercept = fA_dgwt_helper,linetype = 2) +
          # ggplot2::geom_hline(yintercept = afB_dgwt_helper,linetype = 4) +
          ggplot2::scale_shape_manual(values = c(1, 19)) +
          ggplot2::labs(x = "dG folding state A",
              y = "dG folding state B",
              color = color_type,
              shape = paste0("FDR < ", fdr_thres),
              title = "sig. state A folding ddG")
      if (is.factor(vd_plot$color_type)) {
          dfA_dfB_dfAsig <- dfA_dfB_dfAsig + ggplot2::scale_color_brewer(palette = "Set1")
      } else {
          dfA_dfB_dfAsig <- dfA_dfB_dfAsig + ggplot2::scale_color_gradient(low = col_orange, high = col_purple)
      }

      dfA_dfB_dfBsig <- ggplot2::ggplot(vd_plot,
            ggplot2::aes(x = fA_ddg + fA_dgwt_helper,
                y = fB_ddg + fB_dgwt_helper)) +
          ggplot2::geom_density2d(color = 'black', alpha = 0.25) +
          ggplot2::geom_point(ggplot2::aes(color = color_type,shape = fB_ddg_fdr < fdr_thres)) +
          # ggplot2::geom_vline(xintercept = fA_dgwt_helper,linetype = 4) +
          ggplot2::geom_hline(yintercept = fB_dgwt_helper,linetype = 2) +
          ggplot2::scale_shape_manual(values = c(1, 19)) +
          ggplot2::labs(x = "dG folding state A",
              y = "dG folding state B",
              color = color_type,
              shape = paste0("FDR < ", fdr_thres),
              title = "sig. state B folding ddG")
      if (is.factor(vd_plot$color_type)) {
          dfA_dfB_dfBsig <- dfA_dfB_dfBsig + ggplot2::scale_color_brewer(palette = "Set1")
      } else {
          dfA_dfB_dfBsig <- dfA_dfB_dfBsig + ggplot2::scale_color_gradient(low = col_orange, high = col_purple)
      }

      dfA_dfB_sig <- ggplot2::ggplot(vd_plot,
            ggplot2::aes(x = fA_ddg + fA_dgwt_helper,
                y = fB_ddg + fB_dgwt_helper)) +
          ggplot2::geom_density2d(color = 'black', alpha = 0.25) +
          ggplot2::geom_point(ggplot2::aes(color = color_type,shape = fA_ddg_fdr < fdr_thres & fB_ddg_fdr < fdr_thres)) +
          ggplot2::geom_vline(xintercept = fA_dgwt_helper,linetype = 2) +
          ggplot2::geom_hline(yintercept = fB_dgwt_helper,linetype = 2) +
          ggplot2::scale_shape_manual(values = c(1, 19)) +
          ggplot2::labs(x = "dG folding state A",
              y = "dG folding state B",
              color = color_type,
              shape = paste0("FDR < ", fdr_thres),
              title = "sig. folding ddGs both states")
      if (is.factor(vd_plot$color_type)) {
          dfA_dfB_sig <- dfA_dfB_sig + ggplot2::scale_color_brewer(palette = "Set1")
      } else {
          dfA_dfB_sig <- dfA_dfB_sig + ggplot2::scale_color_gradient(low = col_orange, high = col_purple)
      }


      # fA vs b
      dfA_db_dfsig <- ggplot2::ggplot(vd_plot,
            ggplot2::aes(x = fA_ddg + fA_dgwt_helper,
                y = b_ddg + b_dgwt_helper)) +
          ggplot2::geom_density2d(color = 'black', alpha = 0.25) +
          ggplot2::geom_point(ggplot2::aes(color = color_type,shape = fA_ddg_fdr < fdr_thres)) +
          ggplot2::geom_vline(xintercept = fA_dgwt_helper,linetype = 2) +
          # ggplot2::geom_hline(yintercept = b_dgwt_helper,linetype = 4) +
          ggplot2::scale_shape_manual(values = c(1, 19)) +
          ggplot2::labs(x = "dG folding state A",
              y = "dG binding",
              color = color_type,
              shape = paste0("FDR < ", fdr_thres),
              title = "sig. state A folding ddG")
      if (is.factor(vd_plot$color_type)) {
          dfA_db_dfsig <- dfA_db_dfsig + ggplot2::scale_color_brewer(palette = "Set1")
      } else {
          dfA_db_dfsig <- dfA_db_dfsig + ggplot2::scale_color_gradient(low = col_orange, high = col_purple)
      }

      dfA_db_dbsig <- ggplot2::ggplot(vd_plot,
            ggplot2::aes(x = fA_ddg + fA_dgwt_helper,
                y = b_ddg + b_dgwt_helper)) +
          ggplot2::geom_density2d(color = 'black', alpha = 0.25) +
          ggplot2::geom_point(ggplot2::aes(color = color_type,shape = b_ddg_fdr < fdr_thres)) +
          # ggplot2::geom_vline(xintercept = fA_dgwt_helper,linetype = 4) +
          ggplot2::geom_hline(yintercept = b_dgwt_helper,linetype = 2) +
          ggplot2::scale_shape_manual(values = c(1, 19)) +
          ggplot2::labs(x = "dG folding state A",
              y = "dG binding",
              color = color_type,
              shape = paste0("FDR < ", fdr_thres),
              title = "sig. binding ddG")
      if (is.factor(vd_plot$color_type)) {
          dfA_db_dbsig <- dfA_db_dbsig + ggplot2::scale_color_brewer(palette = "Set1")
      } else {
          dfA_db_dbsig <- dfA_db_dbsig + ggplot2::scale_color_gradient(low = col_orange, high = col_purple)
      }

      dfA_db_sig <- ggplot2::ggplot(vd_plot,
            ggplot2::aes(x = fA_ddg + fA_dgwt_helper,
                y = b_ddg + b_dgwt_helper)) +
          ggplot2::geom_density2d(color = 'black', alpha = 0.25) +
          ggplot2::geom_point(ggplot2::aes(color = color_type,shape = fA_ddg_fdr < fdr_thres & b_ddg_fdr < fdr_thres)) +
          ggplot2::geom_vline(xintercept = fA_dgwt_helper,linetype = 2) +
          ggplot2::geom_hline(yintercept = b_dgwt_helper,linetype = 2) +
          ggplot2::scale_shape_manual(values = c(1, 19)) +
          ggplot2::labs(x = "dG folding state A",
              y = "dG binding",
              color = color_type,
              shape = paste0("FDR < ", fdr_thres),
              title = "sig. state A folding & binding ddGs")
      if (is.factor(vd_plot$color_type)) {
          dfA_db_sig <- dfA_db_sig + ggplot2::scale_color_brewer(palette = "Set1")
      } else {
          dfA_db_sig <- dfA_db_sig + ggplot2::scale_color_gradient(low = col_orange, high = col_purple)
      }

      # fB vs b
      dfB_db_dfsig <- ggplot2::ggplot(vd_plot,
            ggplot2::aes(x = fB_ddg + fB_dgwt_helper,
                y = b_ddg + b_dgwt_helper)) +
          ggplot2::geom_density2d(color = 'black', alpha = 0.25) +
          ggplot2::geom_point(ggplot2::aes(color = color_type,shape = fB_ddg_fdr < fdr_thres)) +
          ggplot2::geom_vline(xintercept = fB_dgwt_helper,linetype = 2) +
          # ggplot2::geom_hline(yintercept = b_dgwt_helper,linetype = 4) +
          ggplot2::scale_shape_manual(values = c(1, 19)) +
          ggplot2::labs(x = "dG folding state B",
              y = "dG binding",
              color = color_type,
              shape = paste0("FDR < ", fdr_thres),
              title = "sig. state B folding ddG")
      if (is.factor(vd_plot$color_type)) {
          dfB_db_dfsig <- dfB_db_dfsig + ggplot2::scale_color_brewer(palette = "Set1")
      } else {
          dfB_db_dfsig <- dfB_db_dfsig + ggplot2::scale_color_gradient(low = col_orange, high = col_purple)
      }

      dfB_db_dbsig <- ggplot2::ggplot(vd_plot,
            ggplot2::aes(x = fB_ddg + fB_dgwt_helper,
                y = b_ddg + b_dgwt_helper)) +
          ggplot2::geom_density2d(color = 'black', alpha = 0.25) +
          ggplot2::geom_point(ggplot2::aes(color = color_type,shape = b_ddg_fdr < fdr_thres)) +
          # ggplot2::geom_vline(xintercept = fB_dgwt_helper,linetype = 4) +
          ggplot2::geom_hline(yintercept = b_dgwt_helper,linetype = 2) +
          ggplot2::scale_shape_manual(values = c(1, 19)) +
          ggplot2::labs(x = "dG folding state B",
              y = "dG binding",
              color = color_type,
              shape = paste0("FDR < ", fdr_thres),
              title = "sig. binding ddG")
      if (is.factor(vd_plot$color_type)) {
          dfB_db_dbsig <- dfB_db_dbsig + ggplot2::scale_color_brewer(palette = "Set1")
      } else {
          dfB_db_dbsig <- dfB_db_dbsig + ggplot2::scale_color_gradient(low = col_orange, high = col_purple)
      }

      dfB_db_sig <- ggplot2::ggplot(vd_plot,
            ggplot2::aes(x = fB_ddg + fB_dgwt_helper,
                y = b_ddg + b_dgwt_helper)) +
          ggplot2::geom_density2d(color = 'black', alpha = 0.25) +
          ggplot2::geom_point(ggplot2::aes(color = color_type,shape = fB_ddg_fdr < fdr_thres & b_ddg_fdr < fdr_thres)) +
          ggplot2::geom_vline(xintercept = fB_dgwt_helper,linetype = 2) +
          ggplot2::geom_hline(yintercept = b_dgwt_helper,linetype = 2) +
          ggplot2::scale_shape_manual(values = c(1, 19)) +
          ggplot2::labs(x = "dG folding state B",
              y = "dG binding",
              color = color_type,
              shape = paste0("FDR < ", fdr_thres),
              title = "sig. state B folding & binding ddGs")
      if (is.factor(vd_plot$color_type)) {
          dfB_db_sig <- dfB_db_sig + ggplot2::scale_color_brewer(palette = "Set1")
      } else {
          dfB_db_sig <- dfB_db_sig + ggplot2::scale_color_gradient(low = col_orange, high = col_purple)
      }
      # p <-
      ggplot2::ggsave(gridExtra::grid.arrange(
                          dfA_dfB_dfAsig, dfA_db_dfsig, dfB_db_dfsig,
                          dfA_dfB_dfBsig, dfA_db_dbsig, dfB_db_dbsig,
                          dfA_dfB_sig, dfA_db_sig, dfB_db_sig,
                        nrow = 3,
                        ncol = 3),
          file = file.path(dataset_folder, model_name,
              paste0("results/bootstrap_dG_uncertainty_", color_type, ".pdf")),
          width = 3*4.5,
          height = 9)
    }
  }





  ##############################################################################
  ## plot relationship between different fitness measurements if effects where only driven by folding
  ##############################################################################
  plot_list = list()
  pidx <- 1

  # abundance fitness datasets against one another
  if (varlist[["no_abd_datasets"]] > 1) {
      idx0 <- utils::combn(varlist[["no_abd_datasets"]],2)
      idx <- rbind(parlist[["str_abd"]][idx0[1, ]],
        parlist[["str_abd"]][idx0[2, ]])
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
  for (fi in parlist[["str_abd"]]) {
      for (bi in parlist[["str_bind"]]) {
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

          if (fi == parlist[["str_abd"]][datasets_ab[1]] &
            bi == parlist[["str_bind"]][datasets_ab[2]]) {
              plot_ind <- pidx
          }
          pidx <- pidx + 1
      }
  }

  # binding fitness datasets against one another
  if (varlist[["no_bind_datasets"]] > 1) {
      idx0 <- utils::combn(varlist[["no_bind_datasets"]],2)
      idx <- rbind(parlist[["str_bind"]][idx0[1, ]],
        parlist[["str_bind"]][idx0[2, ]])
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
      file = file.path(dataset_folder, model_name,
        paste0("results/fitness_scatter_dg_relationship_fA",
          parlist[["str_abd"]][datasets_ab[1]], "_fB", parlist[["str_abd"]][datasets_ab[1]], "_", color_type,
                ifelse(stage == "bootstrap", "_boot",""), ".pdf")),
      width = 4,
      height = 3.5)

  # plot all comparisions
  if (varlist[["no_abd_datasets"]] > 1 | varlist[["no_bind_datasets"]] > 1 & length(plot_list) < 20) {
      # p <-
      ggplot2::ggsave(gridExtra::grid.arrange(grobs = plot_list,
                        nrow = ceiling(length(plot_list)/min(3, ceiling(sqrt(length(plot_list))))),
                        ncol = min(3, ceiling(sqrt(length(plot_list))))),
          file = file.path(dataset_folder, model_name,
            paste0("results/fitness_scatter_dg_relationship_all_", color_type,
                ifelse(stage == "bootstrap", "_boot",""), ".pdf")),
          width = 4 * min(3, ceiling(sqrt(length(plot_list)))),
          height = 3.5 * ceiling(length(plot_list)/min(3, ceiling(sqrt(length(plot_list))))),
          limitsize = FALSE)
  }



  ##################################################################
  ## plot relationship between dG and distance to ligand and RSA ###
  ##################################################################
  if (file.exists(file.path(dataset_folder, "data/structural_properties.txt"))) {

    # select single mutants that match a color_type and have at least one measured folding and binding value each
    vd_plot <- vd[!is.na(Pos) & !is.na(color_type), x := rowSums(!is.na(.SD)) > 0, .SDcols = grep("f[0-9]*_fitness", names(vd))][
      x == TRUE, y := rowSums(!is.na(.SD)) > 0, .SDcols = grep("b[0-9]*_fitness", names(vd))][y == TRUE]

    dg_list <- grep("^[bf][AB]?_ddg$", names(vd_plot), value = T)
    for (dg_idx in dg_list) {
      if (stage == "model") {
      x <- vd_plot[, list(color_type = unique(color_type),
                        ddg_type = dg_idx,
                        ddg = mean(unlist(.SD[, 1])),
                        ddg_sd = stats::sd(unlist(.SD[, 1]))),
          Pos, .SDcols = dg_idx]
      } else if (stage == "bootstrap"){
        vd_plot[, paste0(dg_idx, "_sd_reg") := .SD +
              stats::quantile(vd_plot[, unlist(.SD), .SDcols = paste0(dg_idx, "_sd")], 0.25),
              Pos,
              .SDcols = paste0(dg_idx, "_sd")]
        x <- vd_plot[, list(color_type = unique(color_type),
                          ddg_type = dg_idx,
                          ddg = sum(unlist(.SD[, 1])/unlist(.SD[, 2])^2)/sum(1/unlist(.SD[, 2])^2),
                          ddg_sd = stats::sd(unlist(.SD[, 1])),
                          ddg_weight = sum(1/unlist(.SD[, 2])^2)),
            Pos, .SDcols = c(dg_idx, paste0(dg_idx, "_sd_reg"))]
        x[, ddg_weight := ddg_weight / max(ddg_weight)]
      }
      if (dg_idx == dg_list[1]) {
        vd_pos <- x
        names(x) <- gsub("ddg", dg_idx, names(x))
        vd_write <- x
      } else {
        vd_pos <- rbind(vd_pos, x)
        names(x) <- gsub("ddg", dg_idx, names(x))
        x <- x[,.SD, .SDcols = !grepl("color_type", names(x))]
        vd_write <- merge(vd_write, x)
      }
    }

    vd_pos[, ddg_type := factor(ddg_type, levels = dg_list)]

    vd_pos <- merge(
      vd_pos,
      structural_properties[, list(Pos, scHAmin_ligand, RSA_unbound)],
        by = "Pos")
    vd_write <- merge(
      vd_write,
      structural_properties[, list(Pos, scHAmin_ligand, RSA_unbound)],
        by = "Pos")

    vd_pos1 <- melt(vd_pos,
      id.vars = c("Pos", "color_type", grep("ddg", names(vd_pos), value = T)),
      measure.vars = c("scHAmin_ligand", "RSA_unbound"))
    setnames(vd_pos1, "variable", "structural_property")
    vd_pos1[, structural_property := factor(structural_property, levels = c("RSA_unbound", "scHAmin_ligand"))]
    levels(vd_pos1$structural_property) <- c("rel. surface accessibility", "distance to ligand")
    setnames(vd_pos1, "value", "structural_property_value")


    if (stage == "bootstrap") {
      p <- ggplot2::ggplot(vd_pos1,
          ggplot2::aes(structural_property_value, ddg)) +
          ggplot2::geom_pointrange(ggplot2::aes(
            color = color_type,
            ymin = ddg - ddg_sd,
            ymax = ddg + ddg_sd,
            alpha = ddg_weight)) +
          ggplot2::geom_smooth(ggplot2::aes(weight = ddg_weight), span = 0.5, color = "black")
    } else {
      p <- ggplot2::ggplot(vd_pos1,
        ggplot2::aes(structural_property_value, ddg)) +
        ggplot2::geom_pointrange(ggplot2::aes(
          color = color_type,
          ymin = ddg - ddg_sd,
          ymax = ddg + ddg_sd)) +
        ggplot2::geom_smooth(span = 0.5, color = "black")
    }
    p <- p + ggplot2::geom_hline(yintercept = 0, linetype = 3) +
      ggplot2::facet_grid(ddg_type ~ structural_property, scales = "free") +
      ggplot2::labs(y = "ddG values [+- SD over values per position]",
        x = "structural value (%RSA or distance[A])",
        color = color_type,
        alpha = "rel. ddG certainty")
    if (is.factor(vd_plot$color_type)) {
        p <- p + ggplot2::scale_color_brewer(palette = "Set1")
    } else {
        p <- p + ggplot2::scale_color_gradient(low = col_orange, high = col_purple)
    }

    ggplot2::ggsave(p,
      file = file.path(dataset_folder, model_name,
          paste0("results/dG_structural_properties_", color_type,
                  ifelse(stage == "bootstrap", "_boot",""), ".pdf")),
      width = 8,
      height = (parlist[["no_folded_states"]] + 1) * 3.5)

    # write results to file
    utils::write.table(vd_write,
      file = file.path(dataset_folder, model_name,
          paste0("results/dG_structural_properties",
                  ifelse(stage == "bootstrap", "_boot",""), ".txt")),
      quote = F, row.names = F, sep = "\t")
  }
}
