## test tempura package functionality

setwd("tempura")
require(devtools)

load_all()

dataset_folder = "../doubledeepPCA/dg_models/RBD"
model_name = "three_state"


abd = grep("SortSeq", list.files(file.path(dataset_folder, "data")), value = T)
bind = grep("TiteSeq", list.files(file.path(dataset_folder, "data")), value = T)
dg_prepare_datasets(dataset_folder = dataset_folder,
    abundancepca_files = paste0(dataset_folder, "/data/", abd),
    bindingpca_files = paste0(dataset_folder, "/data/", bind),
    fitness_scale = "log"
)

dg_prepare_model(
	dataset_folder = dataset_folder,
	model_name = model_name,
    fix_f_dgwt = TRUE,
    fix_b_dgwt = FALSE
)

dg_estimate_parallel(
    dataset_folder = dataset_folder,
    model_name = model_name,
    iterations = 1:20,
    Ncores = 8
)

for (it in seq(1,5)) {
    print (it)
    dg_estimate(
        dataset_folder = dataset_folder,
        model_name = model_name,
        iteration = it,
        which_test_set = 1,
        maxit = 100,
        trace_optim = TRUE
    )
}

dg_collect_models(
    dataset_folder = dataset_folder,
    model_name = model_name
)



for (it in seq(2,5)) {
    print (it)
    dg_bootstrap(
        dataset_folder = dataset_folder,
        model_name = model_name,
        iteration = it,
        maxit = 20
    )
}

dg_collect_models(
    dataset_folder = dataset_folder,
    model_name = model_name,
    stage = "bootstrap"
)

dg_basic_analyses(
    dataset_folder = dataset_folder,
    model_name = model_name,
    color_type = "type",
    datasets_ab = c(1,1)
)
