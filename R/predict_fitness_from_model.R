#' Predict fitness from dG model
#'
#' @param dg_model0 dG model(s) for which fitness is to be predicted
#' @param variant_data0 variant_data data.table
#' @param varlist varlist
#' @param parlist parlist
#' @param calc_performance logical, if Pearson's R should be calculated for predicted versus measured fitness, default = FALSE
#' @param per_testset logical, if TRUE, predicts fitness only for test set variants (dg_model0 needs to contain test_set variable), default = FALSE
#'
#' @return returns a list with a variant_data data.table with predited fitness values, and a prediction_performance data.table if calc_performance == TRUE
#' @import data.table
#'
#' @export
#'

predict_fitness_from_model <- function(
    dg_model0,
    variant_data0,
    varlist,
    parlist,
    calc_performance = FALSE,
    per_testset = FALSE
) {

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
      dg_model <- dg_model0[test_set == i]
    } else {
      variant_data <- variant_data0
      varxmut <- varlist[["varxmut"]]
      dg_model <- dg_model0
    }


    ## collect folding ddgs & predict fitness in folding datasets
    if (parlist[["no_folded_states"]] == 1) {
        # collect
        f_ddg <- data.table(variant = sapply(X = grep("f_ddg", names(dg_model)),
                                   FUN = function(X) {strsplit(names(dg_model)[X], "_")[[1]][1]}),
                  f_ddg = dg_model[, unlist(.SD), .SDcols = grep("f_ddg", names(dg_model))])
        # calculate ddG per variant
        f_ddg_var <- list(varxmut %*% f_ddg[, f_ddg])
    } else {
        fA_ddg <- data.table(variant = sapply(X = grep("fA_ddg", names(dg_model)),
                                   FUN = function(X) {strsplit(names(dg_model)[X], "_")[[1]][1]}),
                  fA_ddg = dg_model[, unlist(.SD), .SDcols = grep("fA_ddg", names(dg_model))])
        fB_ddg <- data.table(variant = sapply(X = grep("fB_ddg", names(dg_model)),
                                   FUN = function(X) {strsplit(names(dg_model)[X], "_")[[1]][1]}),
                  fB_ddg = dg_model[, unlist(.SD), .SDcols = grep("fB_ddg", names(dg_model))])

        f_ddg_var <- list(
            varxmut %*% fA_ddg[, fA_ddg],
            varxmut %*% fB_ddg[, fB_ddg]
        )
    }
    # predict folding fitness
    for (ds in 1:varlist[["no_abd_datasets"]]) {
        variant_data[, paste0("f", ds, "_pred") := as.numeric(
                convert_dg2foldingfitness(
                  f_ddg = f_ddg_var,
                  f_dgwt = dg_model[, unlist(.SD), .SDcols = grep(paste0("^f[AB]?", ds, "_dgwt"), names(dg_model))],
                  f_fitwt = dg_model[, unlist(.SD), .SDcols = paste0("f", ds, "_fitwt")],
                  f_fit0 = dg_model[, unlist(.SD), .SDcols = paste0("f", ds, "_fit0")],
                  fitness_scale = varlist[["fitness_scale"]],
                  no_folded_states = parlist[["no_folded_states"]]
                )
            )
        ]

        # compute Pearson's R
        if (calc_performance == TRUE) {
          prediction_performance[test_set == i,
            paste0("f", ds) := variant_data[,
                stats::cor(.SD, use = "na.or.complete")[2,1],
                .SDcols = grep(paste0("^f", ds, "_[pf]"), names(variant_data))]]
        }
    }


    ## collect binding ddgs & predict fitness in binding datasets
    b_ddg <- data.table(variant = sapply(X = grep("b_ddg", names(dg_model)),
                                   FUN = function(X) {strsplit(names(dg_model)[X], "_")[[1]][1]}),
                  b_ddg = dg_model[, unlist(.SD), .SDcols = grep("b_ddg", names(dg_model))])
    b_ddg_var <- varxmut %*% b_ddg[, b_ddg]
    for (ds in 1:varlist[["no_bind_datasets"]]) {
        variant_data[, paste0("b", ds, "_pred") := as.numeric(
                convert_dg2bindingfitness(
                  b_ddg = b_ddg_var,
                  f_ddg = f_ddg_var,
                  b_dgwt = dg_model[, unlist(.SD), .SDcols = paste0("b", ds, "_dgwt")],
                  f_dgwt = dg_model[, unlist(.SD), .SDcols = grep(paste0("^bf[AB]?", ds, "_dgwt"), names(dg_model))],
                  b_fitwt = dg_model[, unlist(.SD), .SDcols = paste0("b", ds, "_fitwt")],
                  b_fit0 = dg_model[, unlist(.SD), .SDcols = paste0("b", ds, "_fit0")],
                  fitness_scale = varlist[["fitness_scale"]],
                  no_folded_states = parlist[["no_folded_states"]]
                )
            )
        ]

        if (calc_performance == TRUE) {
          # compute Pearson's R
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
