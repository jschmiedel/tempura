# tempura
An R package to fit thermodynamic models to doubledeepPCA data

Example workflow:

1) dg_prepare_datasets()\
Merges and prepares DiMSum-processed .RData files of folding/abundance and binding measurements. \
These can also be custom .RData files including a variant_data data.table with aa_seq, fitness and sigma columns. \
Creates training/test set splits (default: 10 equally-sized sets)

2) dg_prepare_model()\
Creates necessary parameter sets to fit 3 (unfolded, folded, folded-bound) or 4 (unfolded, folded, folded-binding competent, folded-bound) state thermodynamic models on a prepared dataset.

3) dg_estimate()\
Fits the thermodynamic model. \
Estimates ddG values of folding and binding for each mutation present in the dataset. \
Default is to use 10x cross validation. \
Around 50 models per training/test set are advisable. Use tempura/tempura_cmdln.R to create batchjobs.

4) dg_collect_models()\
Collects models estimated by dg_estimate(). \
Provides plots comparing models and displaying predictive performance on test set variants. \
Computes an average best model across training/test sets.

5) dg_bootstrap()\
Bootstraps the average best model to estimate parameter uncertainty.

6) dg_collect_models(stage = "bootstrap")\
Collects bootstrapped models and provides diagnostic plots on parameter uncertainty.

7) dg_basic_analyses()\
Provides basic analyses plots on inter-relationship of ddG parameters and fitness measurements.

For an example workflow of 3- and 4- state models of the SH3 domain see: \
/nfs/users/blehner/jschmiedel/doubledeepPCA/dg_models/SH3/workflow_SH3.R
