## test tempura package functionality

setwd("tempura")
require(devtools)

load_all()

dataset_folder = "../doubledeepPCA/dg_models/SH3"
model_name = "three_state_test"
model_name = "four_state"

dataset_folder = "/nfs/users/blehner/jschmiedel/doubledeepPCA/dg_models/SH3_regtest"
model_name = "three_lm1"

dg_prepare_datasets(dataset_folder = dataset_folder,
	abundancepca_files = c(paste0(dataset_folder, "/data/01a-GRB2_epPCR_stabilityPCA_aa012_thresholded.RData"),
		paste0(dataset_folder, "/data/01c-GRB2_NM2_stabilityPCA_aa012_thresholded.RData")),
	bindingpca_files = paste0(dataset_folder, "/data/01b-GRB2_epPCR_bindingPCA_aa012_thresholded.RData"),
	wt_seq = "TYVQALFDFDPQEDGELGFRRGDFIHVMDNSDPNWWKGACHGQTGMFPRNYVTPVN*"
)

dg_prepare_model(
	dataset_folder = dataset_folder,
	model_name = model_name,
    lambda = 0.1,
    fix_dgwt = TRUE
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
