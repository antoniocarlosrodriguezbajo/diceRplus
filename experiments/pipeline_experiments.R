library(diceRplus)
library(mclust)
library(R.matlab)
library(scales)

algorithms = c("km", "gmm", "sc", "cmeans", "pam")


########################################################
########################################################
# Meat
########################################################
########################################################

data(Meat)

data <- Meat

mean(data$x)
sd(data$x)

mean_abs_correlation(data$x)
# 0.8426

# UFS
results_UFS <- get_best_top_features(data$x)

result <- run_pipeline(data = data$x,
                       dataset_name = "Meat",
                       top_features = results_UFS$top_features,
                       nk = 5,
                       algorithms = algorithms,
                       UFS_method = "Inf-FS2020",
                       alpha = results_UFS$alpha_best,
                       n_reps=1,
                       seed=100)

exp_ids <- c(result$best_experiment_UFS_id,
             result$best_experiment_NO_UFS_id,
             result$best_experiment_UFS_RPClu_id)

em <- calculate_external_metrics_from_ids(exp_ids, as.integer(data$y))
em$ari

# "Best Clustering UFS: cmeans, Best Consensus: CSPA, Experiment id: 265"
# "Best Clustering NO UFS: sc, Best Consensus: LCA, Experiment id: 287"
# [1] "RPClu with UFS"
# [1] 298
# > em$ari
# [1] 0.3298860 0.2425831 0.6897873

# Run DICE ensemble clustering with all algorithms
exp_id <- run_clustering_ensemble_experiments(data_all = data$x,
                                              dataset  = "Meat",
                                              top_features = NULL,
                                              nk = 5,
                                              algorithms = algorithms,
                                              cons.funs = "LCA",
                                              UFS_method = "NO UFS",
                                              seed=100)

print(exp_id)
em <- calculate_external_metrics_from_ids(exp_id, as.integer(data$y))
em$ari
#
# > print(exp_id)
# [1] 646
# > em <- calculate_external_metrics_from_ids(exp_id, as.integer(data$y))
# > em$ari
# [1] 0.1835855

########################################################
########################################################
# lymphoma
########################################################
########################################################

data(lymphoma)

data <- lymphoma

data$y <- data$y+1

mean(data$x)
sd(data$x)

mean_abs_correlation(data$x)
# 0.2108

# UFS
results_UFS <- get_best_top_features(data$x)

result <- run_pipeline(data = data$x,
                      dataset_name = "Lymphoma",
                      top_features = results_UFS$top_features,
                      nk = 3,
                      algorithms = algorithms,
                      UFS_method = "Inf-FS2020",
                      alpha = results_UFS$alpha_best,
                      n_reps=1,
                      seed=100)

exp_ids <- c(result$best_experiment_UFS_id,
             result$best_experiment_NO_UFS_id,
             result$best_experiment_UFS_RPClu_id)

em <- calculate_external_metrics_from_ids(exp_ids, as.integer(data$y))
em$ari

# "Best Clustering UFS: sc, Best Consensus: LCE, Experiment id: 312"
# "Best Clustering NO UFS: cmeans, Best Consensus: CSPA, Experiment id: 341"
# [1] "RPClu with UFS"
# [1] 349
# > em$ari
# [1] 0.7884032 0.8518607 1.0000000


# Run DICE ensemble clustering with all algorithms
exp_id <- run_clustering_ensemble_experiments(data_all = data$x,
                                              dataset  = "Lymphoma",
                                              top_features = NULL,
                                              nk = 3,
                                              algorithms = algorithms,
                                              cons.funs = "CSPA",
                                              UFS_method = "NO UFS",
                                              seed=100)

print(exp_id)
em <- calculate_external_metrics_from_ids(exp_id, as.integer(data$y))
em$ari
#
# > print(exp_id)
# [1] 647
# > em <- calculate_external_metrics_from_ids(exp_id, as.integer(data$y))
# > em$ari
# [1] 0.9471277

########################################################
########################################################
# warpPIE10P
########################################################
########################################################

data(warpPIE10P)

data <- warpPIE10P

mean(data$x)
sd(data$x)
min(data$x)
max(data$x)

data$x <- apply(data$x, 2, rescale)

mean(data$x)
sd(data$x)
min(data$x)
max(data$x)

mean_abs_correlation(data$x)
# 0.4103

# UFS
results_UFS <- get_best_top_features(data$x)

result <- run_pipeline(data = data$x,
             dataset_name = "warpPIE10P",
             top_features = results_UFS$top_features,
             nk = 10,
             algorithms = algorithms,
             UFS_method = "Inf-FS2020",
             alpha = results_UFS$alpha_best,
             n_reps=1,
             seed=100)

exp_ids <- c(result$best_experiment_UFS_id,
             result$best_experiment_NO_UFS_id,
             result$best_experiment_UFS_RPClu_id)
em <- calculate_external_metrics_from_ids(exp_ids, as.integer(data$y))
em$ari

# "Best Clustering UFS: cmeans, Best Consensus: kmodes, Experiment id: 365"
# "Best Clustering NO UFS: sc, Best Consensus: majority, Experiment id: 386"
# [1] "RPClu with UFS"
# [1] 400
# > em$ari
# [1] 0.1597852 0.1071034 0.5379900


# Run DICE ensemble clustering with all algorithms
exp_id <- run_clustering_ensemble_experiments(data_all = data$x,
                                              dataset  = "warpPIE10P",
                                              top_features = NULL,
                                              nk = 10,
                                              algorithms = algorithms,
                                              cons.funs = "majority",
                                              UFS_method = "NO UFS",
                                              seed=100)

print(exp_id)
em <- calculate_external_metrics_from_ids(exp_id, as.integer(data$y))
em$ari

# > print(exp_id)
# [1] 648
# > em <- calculate_external_metrics_from_ids(exp_id, as.integer(data$y))
# > em$ari
# [1] 0.1139348


########################################################
########################################################
# ALLAML
########################################################
########################################################

data(ALLAML)

data <- ALLAML

mean(data$x)
sd(data$x)
min(data$x)
max(data$x)

mean_abs_correlation(data$x)
# 0.1643

# UFS
results_UFS <- get_best_top_features(data$x)

result <- run_pipeline(data = data$x,
                       dataset_name = "ALLAML",
                       top_features = results_UFS$top_features,
                       nk = 2,
                       algorithms = setdiff(algorithms, "pam"),
                       UFS_method = "Inf-FS2020",
                       alpha = results_UFS$alpha_best,
                       n_reps=1,
                       seed=100)

exp_ids <- c(result$best_experiment_UFS_id,
             result$best_experiment_NO_UFS_id,
             result$best_experiment_UFS_RPClu_id)
em <- calculate_external_metrics_from_ids(exp_ids, as.integer(data$y))
em$ari

# "Best Clustering UFS: cmeans, Best Consensus: LCE, Experiment id: 419"
# "Best Clustering NO UFS: cmeans, Best Consensus: majority, Experiment id: 437"
# [1] "RPClu with UFS"
# [1] 441
#
# em$ari
# [1]  0.02432598  0.23640287 -0.01293935

# Run DICE ensemble clustering with all algorithms
exp_id <- run_clustering_ensemble_experiments(data_all = data$x,
                                              dataset  = "ALLAML",
                                              top_features = NULL,
                                              nk = 2,
                                              algorithms = setdiff(algorithms, "pam"),
                                              cons.funs = "majority",
                                              UFS_method = "NO UFS",
                                              seed=100)

print(exp_id)
em <- calculate_external_metrics_from_ids(exp_id, as.integer(data$y))
em$ari

########################################################
########################################################
# Lung Cancer
########################################################
########################################################

data(lung_cancer)

data <- lung_cancer

mean(data$x)
sd(data$x)
min(data$x)
max(data$x)

mean_abs_correlation(data$x)
# 0.1467

# UFS
results_UFS <- get_best_top_features(data$x)

result <- run_pipeline(data = data$x,
                       dataset_name = "Lung Cancer",
                       top_features = results_UFS$top_features,
                       nk = 5,
                       algorithms = algorithms,
                       UFS_method = "Inf-FS2020",
                       alpha = results_UFS$alpha_best,
                       n_reps=1,
                       seed=100)

exp_ids <- c(result$best_experiment_UFS_id,
             result$best_experiment_NO_UFS_id,
             result$best_experiment_UFS_RPClu_id)
em <- calculate_external_metrics_from_ids(exp_ids, as.integer(data$y))
em$ari

# "Best Clustering NO UFS: pam, Best Consensus: CSPA, Experiment id: 489"
# Best Clustering NO UFS:pam, Best Consensus:CSPA, Experiment id:489"
# [1] "RPClu with UFS"
# [1] 492
# > em$ari
# [1] 0.2092839 0.4814490 0.4013991

# Run DICE ensemble clustering with all algorithms
exp_id <- run_clustering_ensemble_experiments(data_all = data$x,
                                              dataset  = "Lung Cancer",
                                              top_features = NULL,
                                              nk = 5,
                                              algorithms = algorithms,
                                              cons.funs = "CSPA",
                                              UFS_method = "NO UFS",
                                              seed=100)

print(exp_id)
em <- calculate_external_metrics_from_ids(exp_id, as.integer(data$y))
em$ari

# > print(exp_id)
# [1] 650
# > em <- calculate_external_metrics_from_ids(exp_id, as.integer(data$y))
# > em$ari
# [1] 0.3053927

########################################################
########################################################
# Prostate GE
########################################################
########################################################

data(prostate_ge)

data <- prostate_ge

mean(data$x)
sd(data$x)
min(data$x)
max(data$x)

mean_abs_correlation(data$x)
#  0.4118

# UFS
results_UFS <- get_best_top_features(data$x)

result <- run_pipeline(data = data$x,
                       dataset_name = "Prostate GE",
                       top_features = results_UFS$top_features,
                       nk = 2,
                       algorithms = algorithms,
                       UFS_method = "Inf-FS2020",
                       alpha = results_UFS$alpha_best,
                       n_reps=1,
                       seed=100)

exp_ids <- c(result$best_experiment_UFS_id,
             result$best_experiment_NO_UFS_id,
             result$best_experiment_UFS_RPClu_id)
em <- calculate_external_metrics_from_ids(exp_ids, as.integer(data$y))
em$ari

# "Best Clustering UFS: gmm, Best Consensus: CSPA, Experiment id: 500"
# "Best Clustering NO UFS: sc, Best Consensus: kmodes, Experiment id: 528"
# [1] "RPClu with UFS"
# [1] 543
#
# > em$ari
# [1] 0.0007846037 0.0305681440 0.0157535120


# Run DICE ensemble clustering with all algorithms
exp_id <- run_clustering_ensemble_experiments(data_all = data$x,
                                              dataset  = "Prostate GE",
                                              top_features = NULL,
                                              nk = 2,
                                              algorithms = algorithms,
                                              cons.funs = "kmodes",
                                              UFS_method = "NO UFS",
                                              seed=100)

print(exp_id)
em <- calculate_external_metrics_from_ids(exp_id, as.integer(data$y))
em$ari


########################################################
########################################################
# TOX_171
########################################################
########################################################

data(TOX_171)

data <- TOX_171

mean(data$x)
sd(data$x)
min(data$x)
max(data$x)

data$x <- scale(data$x)

mean(data$x)
sd(data$x)
min(data$x)
max(data$x)

mean_abs_correlation(data$x)
#  0.1453

# UFS
results_UFS <- get_best_top_features(data$x)

result <- run_pipeline(data = data$x,
                       dataset_name = "TOX_171",
                       top_features = results_UFS$top_features,
                       nk = 4,
                       algorithms = algorithms,
                       UFS_method = "Inf-FS2020",
                       alpha = results_UFS$alpha_best,
                       n_reps=1,
                       seed=100)

exp_ids <- c(result$best_experiment_UFS_id,
             result$best_experiment_NO_UFS_id,
             result$best_experiment_UFS_RPClu_id)
em <- calculate_external_metrics_from_ids(exp_ids, as.integer(data$y))
em$ari

# [1] "Best Clustering UFS: gmm, Best Consensus: kmodes, Experiment id: 549"
# [1] "Best Clustering NO UFS: km, Best Consensus: LCE, Experiment id: 572"
# 1] "RPClu with UFS"
# [1] 594
# > em$ari
# [1] 0.15973709 0.23570909 0.08408656


# Run DICE ensemble clustering with all algorithms
exp_id <- run_clustering_ensemble_experiments(data_all = data$x,
                                              dataset  = "TOX_171",
                                              top_features = NULL,
                                              nk = 4,
                                              algorithms = algorithms,
                                              cons.funs = "kmodes",
                                              UFS_method = "NO UFS",
                                              seed=100)

print(exp_id)
em <- calculate_external_metrics_from_ids(exp_id, as.integer(data$y))
em$ari

# > print(exp_id)
# [1] 653
# > em <- calculate_external_metrics_from_ids(exp_id, as.integer(data$y))
# > em$ari
# [1] 0.4051207

########################################################
########################################################
# GLIOMA
########################################################
########################################################

data(GLIOMA)

data <- GLIOMA

mean(data$x)
sd(data$x)
min(data$x)
max(data$x)


mean_abs_correlation(data$x)
# [1] 0.3441

# UFS
results_UFS <- get_best_top_features(data$x)

result <- run_pipeline(data = data$x,
                       dataset_name = "GLIOMA",
                       top_features = results_UFS$top_features,
                       nk = 4,
                       algorithms = algorithms,
                       UFS_method = "Inf-FS2020",
                       alpha = results_UFS$alpha_best,
                       n_reps=1,
                       seed=100)

exp_ids <- c(result$best_experiment_UFS_id,
             result$best_experiment_NO_UFS_id,
             result$best_experiment_UFS_RPClu_id)
em <- calculate_external_metrics_from_ids(exp_ids, as.integer(data$y))
em$ari

# [1] "Best Clustering UFS: cmeans, Best Consensus: kmodes, Experiment id: 610"
# [1] "Best Clustering NO UFS: pam, Best Consensus: kmodes, Experiment id: 640"
# [1] "RPClu with UFS"
# [1] 645
# [1] 0.13115664 0.36979329 0.05237299

# Run DICE ensemble clustering with all algorithms
exp_id <- run_clustering_ensemble_experiments(data_all = data$x,
                                              dataset  = "GLIOMA",
                                              top_features = NULL,
                                              nk = 4,
                                              algorithms = algorithms,
                                              cons.funs = "kmodes",
                                              UFS_method = "NO UFS",
                                              seed=100)

print(exp_id)
em <- calculate_external_metrics_from_ids(exp_id, as.integer(data$y))
em$ari
# > print(exp_id)
# [1] 652
# > em <- calculate_external_metrics_from_ids(exp_id, as.integer(data$y))
# > em$ari
# [1] 0.4051207
