#' Prepare DiMSum-processed doubledeepPCA datasets for free energy estimation
#'
#'
#' @param dataset_folder absolute path to the dataset folder, is created if non-existent
#' @param abundancepca_files absolute path to the DiMSum-processed abundancePCA .RData files containing aa_seq (or aa_subs), fitness and sigma cols; multiple files will be treated as independent measurements
#' @param bindingpca_files absolute path to the DiMSum-processed bindingPCA .RData files containing aa_seq (or aa_subs), fitness and sigma cols; multiple files will be used as independent measurements
#' @param wt_seq a character string of the wild-type/reference aa sequence used to extract aa_subs from aa_seq; if wt_seq is empty, assumes aa_subs column is present
#' @param train_test_split c(A, B), Ax test split of variants, where each split has size B; default c(10, 0.1), i.e. 10x cross validation
#' @param fitness_scale either "lin" or "log"; for "lin", assumes fitness = [0,Inf], with wild-type fitness 1; for "log", assumes fitness = [-Inf, Inf], with wild-type fitness = 0, default = "lin"
#'
#' @return writes a .RData file to $dataset_folder/data/fitness_dataset.RData containing the varlist list with all necessary variables to compute dg models for the dataset
#' @import data.table
#' @import Matrix
#' @export
dg_prepare_datasets = function(
	dataset_folder,
	abundancepca_files,
	bindingpca_files,
	wt_seq = "",
	train_test_split = c(10, 0.1),
	fitness_scale = "lin"
) {

	tts <- test_set <- aa_seq <- aa_subs <- mutation <- NULL

	## create dataset folders, if not already present
	dir.create(file.path(dataset_folder))
	dir.create(file.path(dataset_folder, "data"))

	## load and merge datasets
	all_files = list()
	no_abd_datasets <- length(abundancepca_files)
	no_bind_datasets <- length(bindingpca_files)

  if (wt_seq != "") {
    common_col = "aa_seq"
  } else {
    common_col = "aa_subs"
  }

	# load abundancePCA files
	for (i in seq_along(abundancepca_files)) {
		load(abundancepca_files[i])
		names(all_variants)[grep("^fitness$", names(all_variants))] <- paste0("f", i, "_fitness")
		names(all_variants)[grep("^sigma$", names(all_variants))] <- paste0("f", i, "_sigma")
		all_files[[i]] <- all_variants[, #only allow missense variants
			.SD,
			.SDcols = c(common_col,
				paste0("f", i, "_fitness"),
				paste0("f", i, "_sigma"))]
	}
	# load bindingPCA files
	for (i in seq_along(bindingpca_files)) {
		load(bindingpca_files[i])
		names(all_variants)[grep("^fitness$", names(all_variants))] <- paste0("b", i, "_fitness")
		names(all_variants)[grep("^sigma$", names(all_variants))] <- paste0("b", i, "_sigma")
		all_files[[i + no_abd_datasets]] <- all_variants[,  #only allow missense variants
			.SD,
			.SDcols = c(common_col,
				paste0("b", i, "_fitness"),
				paste0("b", i, "_sigma"))]
	}
	# merge files
	variant_data <- all_files[[1]]
	for (i in 2:length(all_files)) {
		variant_data <- merge(variant_data, all_files[[i]], by = common_col, all = T)
	}


	## define train/test split of data
	variant_data <- variant_data[sample(.N)]
	variant_data[, tts := (1:.N)/.N]
	split_vec <- seq(train_test_split[2], train_test_split[2]*train_test_split[1], train_test_split[2])
	# increment by 1 to have test_sets run from 1 to train_test_split[1]
	variant_data[, test_set := findInterval(tts, split_vec, rightmost.closed = T) + 1]
	# if train_test_split[2]*train_test_split[1] < 1, set fraction of variants not assigned a test_set to NA
	variant_data[test_set > train_test_split[1], test_set := NA]
	variant_data[, tts := NULL]


  if (wt_seq != "") {
    ## extract aa mutations from aa_seq
    wt_seq_split <- strsplit(wt_seq, "")[[1]]
    variant_data[, aa_subs := paste0(stringr::str_pad(which(strsplit(aa_seq, "")[[1]] != wt_seq_split), ceiling(log10(length(wt_seq_split))), "left", "0"),
        strsplit(aa_seq, "")[[1]][which(strsplit(aa_seq, "")[[1]] != wt_seq_split)],
        collapse = "_"),
      aa_seq]
  }
	data.table::setkey(variant_data, aa_subs)


	## extract list of all mutations
	mut_list <- data.table::data.table(mutation = sort(unique(unlist(lapply(variant_data[, aa_subs], strsplit, "_")))))


	## create 'variant X mutations' matrix; this is a slow implementation
	varxmut <- Matrix::Matrix(0, nrow = variant_data[, .N], ncol = nrow(mut_list), sparse = T)
	for (i in 1:nrow(variant_data)) {
	  x <- variant_data[i, strsplit(aa_subs, "_")[[1]]]
	  if (length(x) > 0){
	    varxmut[i, mut_list[, mutation] %in% x] <- 1
	  }
	}


    ## collect all variables
    varlist <- list(variant_data = variant_data,
				mut_list = mut_list,
				varxmut = varxmut,
				no_abd_datasets = no_abd_datasets,
				no_bind_datasets = no_bind_datasets,
				fitness_scale = fitness_scale)


    ## save varlist to .RData file
    save(varlist, file = file.path(dataset_folder, "data/fitness_dataset.RData"))

}
