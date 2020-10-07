#!/usr/bin/env Rscript

# Rscript to run the computation-heavy model and parameter uncertainty estimation stages of tempura via the command line

###########################
### CHECK DEPENDENCIES
###########################

#R packages
required_packages <- c(
  "tempura", "optparse")
missing_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(missing_packages)!=0){
  stop(paste0("Required R packages not installed. Please install the following packages: ", paste(missing_packages, sep = ", ")), call. = FALSE)
}


# read in varlist file deposited in the copied folder structure
option_list <- list(
  optparse::make_option(
    opt_str = c("-s", "--stage"),
    default = "model",
    type = "character",
    help = "modeling stage, 'model' for initial dG model estimation, 'bootstrap' for bootstrapping of dG parameters, 'pll' for profile-likelihood estimates of dG parameters"
  ),
  optparse::make_option(
    opt_str = c("-d", "--dataset_folder"),
    type = "character",
    help = "absolute path to dataset folder"
  ),
  optparse::make_option(
    opt_str = c("-m", "--model_name"),
    type = "character",
    help = "name of dg model to be computed"
  ),
  optparse::make_option(
    opt_str = c("-t", "--which_test_set"),
    default = 0,
    type = "integer"
  ),
  optparse::make_option(
    opt_str = c("-i", "--iteration"),
    default = 1,
    type = "integer"
  ),
  optparse::make_option(
    opt_str = c("--maxit"),
    default = 1e4,
    type = "integer"
  )
)

opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

library(tempura)

if (opt$stage == "model") {
  dg_estimate(
    dataset_folder = opt$dataset_folder,
    model_name = opt$model_name,
    which_test_set = opt$which_test_set,
    iteration = opt$iteration,
    return_model = FALSE,
    maxit = opt$maxit
  )
}

