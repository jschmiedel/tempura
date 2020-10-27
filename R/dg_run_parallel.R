#' A wrapper for parallel exectuion of dg_estimate() or dg_bootstrap()
#'
#'
#' @param dataset_folder absolute path to the dataset folder, is created if non-existent
#' @param model_name name of the model that should be computed on the dataset
#' @param stage "model" or "bootstrap"
#' @param which_test_set integer, (for stage = "model" only): test_set to exclude from model training, default = 0, i.e. will be chosen depending on iteration parameter
#' @param iterations a vector of integers to iterate over in dg_estimate(iteration = iterations[i])
#' @param Ncores integer, how many cores to use for parallelization
#'
#' @return writes the model parameters as .txt file to $dataset_folder/$model_name/tmp/dg_model_$testset_$iteration
#' @import data.table
#' @import Matrix
#'
#' @export
#'
dg_run_parallel <- function(
  dataset_folder,
  model_name,
  which_test_set = 0,
  iterations,
  Ncores = 8,
  stage = "model"
) {

  i <- NULL

  doMC::registerDoMC(cores = Ncores)

  `%dopar%` <- foreach::`%dopar%`

  if (stage == "model") {
    foreach::foreach(i = 1:length(iterations)) %dopar% {
      dg_estimate(
          dataset_folder = dataset_folder,
          model_name = model_name,
          iteration = iterations[i],
          which_test_set = which_test_set
      )
    }
  } else if (stage == "bootstrap") {
    foreach::foreach(i = 1:length(iterations)) %dopar% {
      dg_bootstrap(
          dataset_folder = dataset_folder,
          model_name = model_name,
          iteration = iterations[i]
      )
    }
  }



}
