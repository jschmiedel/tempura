#' Predict fitness from computed dG model
#'
#' @param dg_model0 dG model(s) for which fitness is to be predicted, either long or wide format
#' @param variant_data0 variant_data data.table
#' @param varlist varlist
#' @param parlist parlist
#' @param calc_performance logical, if Pearson's R should be calculated for predicted versus measured fitness, default = FALSE
#' @param per_testset logical, if TRUE, predicts fitness only for test set variants (dg_model0 needs to contain test_set variable and is expected to be in wide format), default = FALSE
#' @param overwrite logical, overwrite predicted fitness values if already present in variant_data0?, default = TRUE
#' @param predfit_ddgvar logical, predict fitness values for a set of variations of ddg values (setting one to zero, holding parameters equal etc, used for plotting analysis), default = FALSE
#'
#' @return returns a list with a variant_data data.table with predited fitness values, and a prediction_performance data.table if calc_performance == TRUE
#' @import data.table
#'
#' @export
#'

dg__fitness_from_model <- function(
    dg_model0,
    variant_data0,
    varlist,
    parlist,
    calc_performance = FALSE,
    per_testset = FALSE,
    overwrite = TRUE,
    predfit_ddgvar = FALSE
) {

  test_set <- parameter <- value <- aa_subs <- NULL

  if (per_testset == TRUE) {
    idx = sort(dg_model0[, test_set])
  } else {idx = 0}

  if (calc_performance == TRUE) {
    prediction_performance = data.table(test_set = idx)
  }

  for (i in idx) {
    if (per_testset == TRUE) {
      variant_data <- variant_data0[test_set == i]
      varxmut <- varlist[["varxmut"]][variant_data0[, which(test_set == i)], ]
      # cast long
      dg_model <- data.table(parameter = names(dg_model0), value = dg_model0[test_set == i, unlist(.SD)])
    } else {
      variant_data <- variant_data0
      varxmut <- varlist[["varxmut"]]
      if (nrow(dg_model0) < ncol(dg_model0)) { #wide format -> cast long
          dg_model <- data.table(parameter = names(dg_model0), value = dg_model0[, unlist(.SD)])
      } else {
          dg_model <- dg_model0
      }
    }


    ## collect folding ddgs & predict fitness in folding datasets
    if (parlist[["no_folded_states"]] == 1) {
        f_ddg_var <- list(varxmut %*% dg_model[grep("f_ddg", parameter), value])
        f_ddg_var0 <- list(rep(0, nrow(varxmut)))
    } else {
        f_ddg_var <- list(
            varxmut %*% dg_model[grep("fA_ddg", parameter), value],
            varxmut %*% dg_model[grep("fB_ddg", parameter), value]
        )
        f_ddg_var0 <- list(rep(0, nrow(varxmut)), rep(0, nrow(varxmut)))
    }

    # extract folding dgwt values
    f_dgwt = list()
    for (ds in 1:varlist[["no_abd_datasets"]]) {
      if (parlist[["fix_dgwt"]] == FALSE) {
        f_dgwt[[ds]] <- dg_model[grep(paste0("^f[AB]?", ds, "_dgwt"), parameter), value]
      } else {
        f_dgwt[[ds]] <- dg_model[grep("^f[AB]?_dgwt", parameter), value]
      }
    }

    # predict folding fitness
    for (ds in 1:varlist[["no_abd_datasets"]]) {

      # predict folding fitness
      if (overwrite == TRUE | length(grep(paste0("^f", ds, "_pred$"), names(variant_data))) == 0) {
        variant_data[, paste0("f", ds, "_pred") := as.numeric(
                convert_dg2foldingfitness(
                  f_ddg_var = f_ddg_var,
                  f_dgwt = f_dgwt[[ds]],
                  f_fitwt = dg_model[grep(paste0("f", ds, "_fitwt"), parameter), value],
                  f_fit0 = dg_model[grep(paste0("f", ds, "_fit0"), parameter), value],
                  fitness_scale = varlist[["fitness_scale"]],
                  no_folded_states = parlist[["no_folded_states"]]
                )
            )
        ]
      }

      if (predfit_ddgvar == TRUE & parlist[["no_folded_states"]] == 2) {

        # predict folding fitness if folding ddg values of state B are set equal to state A
        if (overwrite == TRUE | length(grep(paste0("^f", ds, "_pred_fBasfA$"), names(variant_data))) == 0) {
          f_ddg_var_helper <- list(
            varxmut %*% dg_model[grep("fA_ddg", parameter), value],
            varxmut %*% dg_model[grep("fA_ddg", parameter), value]
          )

          variant_data[, paste0("f", ds, "_pred_fBasfA") := as.numeric(
              convert_dg2foldingfitness(
                f_ddg_var = f_ddg_var_helper,
                f_dgwt = f_dgwt[[ds]],
                f_fitwt = dg_model[grep(paste0("f", ds, "_fitwt"), parameter), value],
                f_fit0 = dg_model[grep(paste0("f", ds, "_fit0"), parameter), value],
                fitness_scale = varlist[["fitness_scale"]],
                no_folded_states = parlist[["no_folded_states"]]
              )
            )
          ]
        }

        # predict folding fitness if folding ddg values of state A are set equal to state B
        if (overwrite == TRUE | length(grep(paste0("^f", ds, "_pred_fAasfB$"), names(variant_data))) == 0) {
          f_ddg_var_helper <- list(
            varxmut %*% dg_model[grep("fB_ddg", parameter), value],
            varxmut %*% dg_model[grep("fB_ddg", parameter), value]
          )

          variant_data[, paste0("f", ds, "_pred_fAasfB") := as.numeric(
              convert_dg2foldingfitness(
                f_ddg_var = f_ddg_var_helper,
                f_dgwt = f_dgwt[[ds]],
                f_fitwt = dg_model[grep(paste0("f", ds, "_fitwt"), parameter), value],
                f_fit0 = dg_model[grep(paste0("f", ds, "_fit0"), parameter), value],
                fitness_scale = varlist[["fitness_scale"]],
                no_folded_states = parlist[["no_folded_states"]]
              )
            )
          ]
        }

        # predict folding fitness if only folding ddg values of state A are non-zero
        if (overwrite == TRUE | length(grep(paste0("^f", ds, "_pred_fBddg0$"), names(variant_data))) == 0) {
          f_ddg_var_helper <- list(
            varxmut %*% dg_model[grep("fA_ddg", parameter), value],
            rep(0, nrow(varxmut))
          )

          variant_data[, paste0("f", ds, "_pred_fBddg0") := as.numeric(
              convert_dg2foldingfitness(
                f_ddg_var = f_ddg_var_helper,
                f_dgwt = f_dgwt[[ds]],
                f_fitwt = dg_model[grep(paste0("f", ds, "_fitwt"), parameter), value],
                f_fit0 = dg_model[grep(paste0("f", ds, "_fit0"), parameter), value],
                fitness_scale = varlist[["fitness_scale"]],
                no_folded_states = parlist[["no_folded_states"]]
              )
            )
          ]
        }

        # predict folding fitness if only folding ddg values of state A are non-zero
        if (overwrite == TRUE | length(grep(paste0("^f", ds, "_pred_fAddg0$"), names(variant_data))) == 0) {
          f_ddg_var_helper <- list(
            rep(0, nrow(varxmut)),
            varxmut %*% dg_model[grep("fB_ddg", parameter), value]
          )

          variant_data[, paste0("f", ds, "_pred_fAddg0") := as.numeric(
              convert_dg2foldingfitness(
                f_ddg_var = f_ddg_var_helper,
                f_dgwt = f_dgwt[[ds]],
                f_fitwt = dg_model[grep(paste0("f", ds, "_fitwt"), parameter), value],
                f_fit0 = dg_model[grep(paste0("f", ds, "_fit0"), parameter), value],
                fitness_scale = varlist[["fitness_scale"]],
                no_folded_states = parlist[["no_folded_states"]]
              )
            )
          ]
        }
      }

      # compute Pearson's R between predictions and measurements
      if (calc_performance == TRUE) {
        prediction_performance[test_set == i,
          paste0("f", ds) := variant_data[,
              stats::cor(.SD, use = "na.or.complete")[2,1],
              .SDcols = grep(paste0("^f", ds, "_[pf]"), names(variant_data))]]
      }
    }


    ## collect binding ddgs & predict fitness in binding datasets
    b_ddg_var <- varxmut %*% dg_model[grep("b_ddg", parameter), value]
    b_ddg_var0 <- rep(0, nrow(varxmut))


    # extract folding dgwt values
    b_dgwt = c()
    bf_dgwt = list()
    for (ds in 1:varlist[["no_bind_datasets"]]) {
      if (parlist[["fix_dgwt"]] == FALSE) {
        b_dgwt[ds] <- dg_model[grep(paste0("b", ds, "_dgwt"), parameter), value]
        bf_dgwt[[ds]] <- dg_model[grep(paste0("^f[AB]?", ds, "_dgwt"), parameter), value]
      } else {
        b_dgwt[ds] <- dg_model[grep("b_dgwt", parameter), value]
        bf_dgwt[[ds]] <- dg_model[grep("^f[AB]?_dgwt", parameter), value]
      }
    }


    for (ds in 1:varlist[["no_bind_datasets"]]) {

      # predict binding fitness
      if (overwrite == TRUE | length(grep(paste0("^b", ds, "_pred$"), names(variant_data))) == 0) {
        variant_data[, paste0("b", ds, "_pred") := as.numeric(
                convert_dg2bindingfitness(
                  b_ddg_var = b_ddg_var,
                  f_ddg_var = f_ddg_var,
                  b_dgwt = b_dgwt[ds],
                  f_dgwt = bf_dgwt[[ds]],
                  b_fitwt = dg_model[grep(paste0("b", ds, "_fitwt"), parameter), value],
                  b_fit0 = dg_model[grep(paste0("b", ds, "_fit0"), parameter), value],
                  fitness_scale = varlist[["fitness_scale"]],
                  no_folded_states = parlist[["no_folded_states"]]
                )
            )
        ]
      }

      if (predfit_ddgvar == TRUE) {

        # predict binding fitness when only folding ddg values are non-zero
        if (overwrite == TRUE | length(grep(paste0("^b", ds, "_pred_bddg0$"), names(variant_data))) == 0) {
          variant_data[, paste0("b", ds, "_pred_bddg0") := as.numeric(
              convert_dg2bindingfitness(
                b_ddg_var = b_ddg_var0,
                f_ddg_var = f_ddg_var,
                b_dgwt = b_dgwt[ds],
                f_dgwt = bf_dgwt[[ds]],
                b_fitwt = dg_model[grep(paste0("b", ds, "_fitwt"), parameter), value],
                b_fit0 = dg_model[grep(paste0("b", ds, "_fit0"), parameter), value],
                fitness_scale = varlist[["fitness_scale"]],
                no_folded_states = parlist[["no_folded_states"]]
              )
            )
          ]
        }

        # predict binding fitness when only binding ddg values are non-zero
        if (overwrite == TRUE | length(grep(paste0("^b", ds, "_pred_fddg0$"), names(variant_data))) == 0) {
          variant_data[, paste0("b", ds, "_pred_fddg0") := as.numeric(
              convert_dg2bindingfitness(
                b_ddg_var = b_ddg_var,
                f_ddg_var = f_ddg_var0,
                b_dgwt = b_dgwt[ds],
                f_dgwt = bf_dgwt[[ds]],
                b_fitwt = dg_model[grep(paste0("b", ds, "_fitwt"), parameter), value],
                b_fit0 = dg_model[grep(paste0("b", ds, "_fit0"), parameter), value],
                fitness_scale = varlist[["fitness_scale"]],
                no_folded_states = parlist[["no_folded_states"]]
              )
            )
          ]
        }

        if (parlist[["no_folded_states"]] == 2) {

          # predict folding fitness if folding ddg values of state B are set equal to state A
          if (overwrite == TRUE | length(grep(paste0("^b", ds, "_pred_fBasfA$"), names(variant_data))) == 0) {
            f_ddg_var_helper <- list(
              varxmut %*% dg_model[grep("fA_ddg", parameter), value],
              varxmut %*% dg_model[grep("fA_ddg", parameter), value]
            )

            variant_data[, paste0("b", ds, "_pred_fBasfA") := as.numeric(
                convert_dg2bindingfitness(
                  b_ddg_var = b_ddg_var0,
                  f_ddg_var = f_ddg_var_helper,
                  b_dgwt = b_dgwt[ds],
                  f_dgwt = bf_dgwt[[ds]],
                  b_fitwt = dg_model[grep(paste0("b", ds, "_fitwt"), parameter), value],
                  b_fit0 = dg_model[grep(paste0("b", ds, "_fit0"), parameter), value],
                  fitness_scale = varlist[["fitness_scale"]],
                  no_folded_states = parlist[["no_folded_states"]]
                )
              )
            ]
          }

          # predict folding fitness if folding ddg values of state A are set equal to state B
          if (overwrite == TRUE | length(grep(paste0("^b", ds, "_pred_fAasfB$"), names(variant_data))) == 0) {
            f_ddg_var_helper <- list(
              varxmut %*% dg_model[grep("fB_ddg", parameter), value],
              varxmut %*% dg_model[grep("fB_ddg", parameter), value]
            )

            variant_data[, paste0("b", ds, "_pred_fAasfB") := as.numeric(
                convert_dg2bindingfitness(
                  b_ddg_var = b_ddg_var0,
                  f_ddg_var = f_ddg_var_helper,
                  b_dgwt = b_dgwt[ds],
                  f_dgwt = bf_dgwt[[ds]],
                  b_fitwt = dg_model[grep(paste0("b", ds, "_fitwt"), parameter), value],
                  b_fit0 = dg_model[grep(paste0("b", ds, "_fit0"), parameter), value],
                  fitness_scale = varlist[["fitness_scale"]],
                  no_folded_states = parlist[["no_folded_states"]]
                )
              )
            ]
          }

          # predict folding fitness if only folding ddg values of state A are non-zero
          if (overwrite == TRUE | length(grep(paste0("^b", ds, "_pred_fBddg0$"), names(variant_data))) == 0) {
            f_ddg_var_helper <- list(
              varxmut %*% dg_model[grep("fA_ddg", parameter), value],
              rep(0, nrow(varxmut))
            )

            variant_data[, paste0("b", ds, "_pred_fBddg0") := as.numeric(
                convert_dg2bindingfitness(
                  b_ddg_var = b_ddg_var0,
                  f_ddg_var = f_ddg_var_helper,
                  b_dgwt = b_dgwt[ds],
                  f_dgwt = bf_dgwt[[ds]],
                  b_fitwt = dg_model[grep(paste0("b", ds, "_fitwt"), parameter), value],
                  b_fit0 = dg_model[grep(paste0("b", ds, "_fit0"), parameter), value],
                  fitness_scale = varlist[["fitness_scale"]],
                  no_folded_states = parlist[["no_folded_states"]]
                )
              )
            ]
          }

          # predict folding fitness if only folding ddg values of state A are non-zero
          if (overwrite == TRUE | length(grep(paste0("^b", ds, "_pred_fAddg0$"), names(variant_data))) == 0) {
            f_ddg_var_helper <- list(
              rep(0, nrow(varxmut)),
              varxmut %*% dg_model[grep("fB_ddg", parameter), value]
            )

            variant_data[, paste0("b", ds, "_pred_fAddg0") := as.numeric(
                convert_dg2bindingfitness(
                  b_ddg_var = b_ddg_var0,
                  f_ddg_var = f_ddg_var_helper,
                  b_dgwt = b_dgwt[ds],
                  f_dgwt = bf_dgwt[[ds]],
                  b_fitwt = dg_model[grep(paste0("b", ds, "_fitwt"), parameter), value],
                  b_fit0 = dg_model[grep(paste0("b", ds, "_fit0"), parameter), value],
                  fitness_scale = varlist[["fitness_scale"]],
                  no_folded_states = parlist[["no_folded_states"]]
                )
              )
            ]
          }
        }

      }



      # compute Pearson's R between predictions and measurements
      if (calc_performance == TRUE) {
        prediction_performance[test_set == i,
            paste0("b", ds) := variant_data[,
                stats::cor(.SD, use = "na.or.complete")[2,1],
                .SDcols = grep(paste0("^b", ds, "_[pf]"), names(variant_data))]]
      }
    }

    # collect variant_data subsets
    if (per_testset == TRUE) {
      if (i == idx[1]) {
        all_data = variant_data
      } else {
        all_data = rbind(all_data, variant_data)
      }
    }
  }

  # return variant data and prediction performance
  if (per_testset == TRUE) {
    data.table::setkey(all_data, aa_subs)
    return_list = list(variant_data = all_data)
  } else {
    return_list = list(variant_data = variant_data)
  }
  if (calc_performance == TRUE) {
    return_list = c(return_list, list(prediction_performance = prediction_performance))
  }
  return(return_list)
}
