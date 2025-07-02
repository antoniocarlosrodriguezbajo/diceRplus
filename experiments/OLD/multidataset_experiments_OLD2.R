library(diceRplus)
library(mclust)

UFS_Methods <- list(
  "Inf-FS2020" = TRUE
)

algorithms = c("nmf", "hc", "diana", "km", "pam", "ap", "sc",
               "gmm", "block", "som", "cmeans", "hdbscan")

calculate_external_metrics <- function (E, k, ref_labels) {
  # List of consensus methods and their specific arguments
  consensus_methods <- list(
    kmodes   = function(E) k_modes(E, is.relabelled = FALSE, seed = 1),
    majority = function(E) majority_voting(E, is.relabelled = FALSE),
    lce      = function(E) LCE(E, k, sim.mat = "cts"),
    lca      = function(E) LCA(E, is.relabelled = FALSE, seed = 1)
  )

  # Apply each consensus method
  consensus_labels <- lapply(consensus_methods, function(f) f(E))

  # Relabel only the methods that need it (adjust as needed)
  aligned_labels <- list(
    kmodes   = relabel_class(consensus_labels$kmodes, ref_labels),
    majority = relabel_class(consensus_labels$majority, ref_labels),
    lce      = relabel_class(consensus_labels$lce, ref_labels),
    lca      = relabel_class(consensus_labels$lca, ref_labels)
  )

  # Evaluate results
  eval_results <- lapply(aligned_labels, function(labels) list(
    confmat = ev_confmat(labels, ref_labels),
    nmi     = ev_nmi(labels, ref_labels),
    ari     = adjustedRandIndex(labels, ref_labels)
  ))

  return(eval_results)
}

# Implementación simplificada del algoritmo CGUFS
cgufs_selection <- function(data, k1=5, k2=5,quantile=0.75) {
  # Paso 1: Agrupación de muestras
  sample_clusters <- kmeans(data, centers=k1)$cluster

  # Paso 2: Cálculo de importancia de características
  feature_importance <- apply(data, 2, function(x) {
    summary(aov(x ~ as.factor(sample_clusters)))[[1]][1,4]
  })

  # Paso 3: Selección adaptativa
  selected_features <- which(feature_importance > quantile(feature_importance, quantile))
  return(selected_features)
}

prepare_consensus_evaluation <- function(data,cluster_labels, alg_name="RPGMMClu") {
  num_labels <- length(cluster_labels)

  # Create structure for consensus_evaluate
  cc_data <- array(cluster_labels, dim = c(num_labels, 1, 1, 1))
  # row_names <- rownames(data)
  row_names <- if (!is.null(rownames(data))) rownames(data) else seq_len(nrow(data))
  dimnames(cc_data) <- list(
    row_names,  # Primer nivel de nombres: nombres de las filas de Meat$x
    "R1",       # Repetition
    alg_name, # Clustering algorithm
    "5"         # Number of clusters
  )
  # Structure
  return(cc_data)
}

############################
# Meat
############################
data(Meat)

data <- Meat

data$x <- scale(data$x)
mean(data$x)
sd(data$x)

# Run just once
UFS_results_file <- "experiments/UFS_Meat_SC_Inf-FS2020.RData"
UFS_results <- runUFS(data$x, UFS_Methods)
save(UFS_results, file = UFS_results_file)

load(UFS_results_file)

method = "Inf-FS2020"
top_features <- UFS_results$Results[[method]]$Result[[3]]

sprintf("The %s method selected %d features.", method, length(top_features))

B <- 100
execution_time <- system.time(RPClu_results <- RPClu_parallel(data$x[, top_features],
                                                                 clust_fun = "gmm",
                                                                 g=5,
                                                                 B=B,
                                                                 verb=TRUE))["elapsed"]

E <- as.matrix(RPClu_results$clusterings)

k=5
ref_labels <- as.integer(data$y)

external_metrics <- calculate_external_metrics(E,k,ref_labels)

external_metrics

# > external_metrics
# $kmodes
# $kmodes$confmat
# # A tibble: 13 × 3
# .metric              .estimator .estimate
# <chr>                <chr>          <dbl>
#   1 accuracy             multiclass     0.823
# 2 kap                  multiclass     0.775
# 3 sens                 macro          0.843
# 4 spec                 macro          0.954
# 5 ppv                  macro          0.847
# 6 npv                  macro          0.954
# 7 mcc                  multiclass     0.776
# 8 j_index              macro          0.797
# 9 bal_accuracy         macro          0.898
# 10 detection_prevalence macro          0.2
# 11 precision            macro          0.847
# 12 recall               macro          0.843
# 13 f_meas               macro          0.843
#
# $kmodes$nmi
# [1] 0.7794098
#
# $kmodes$ari
# [1] 0.666964
#
#
# $majority
# $majority$confmat
# # A tibble: 13 × 3
# .metric              .estimator .estimate
# <chr>                <chr>          <dbl>
#   1 accuracy             multiclass     0.814
# 2 kap                  multiclass     0.764
# 3 sens                 macro          0.836
# 4 spec                 macro          0.951
# 5 ppv                  macro          0.872
# 6 npv                  macro          0.958
# 7 mcc                  multiclass     0.786
# 8 j_index              macro          0.787
# 9 bal_accuracy         macro          0.894
# 10 detection_prevalence macro          0.2
# 11 precision            macro          0.872
# 12 recall               macro          0.836
# 13 f_meas               macro          0.821
#
# $majority$nmi
# [1] 0.8084113
#
# $majority$ari
# [1] 0.6827528
#
#
# $lce
# $lce$confmat
# # A tibble: 13 × 3
# .metric              .estimator .estimate
# <chr>                <chr>          <dbl>
#   1 accuracy             multiclass     0.745
# 2 kap                  multiclass     0.676
# 3 sens                 macro          0.778
# 4 spec                 macro          0.933
# 5 ppv                  macro          0.684
# 6 npv                  macro          0.948
# 7 mcc                  multiclass     0.730
# 8 j_index              macro          0.711
# 9 bal_accuracy         macro          0.855
# 10 detection_prevalence macro          0.2
# 11 precision            macro          0.684
# 12 recall               macro          0.778
# 13 f_meas               macro          0.713
#
# $lce$nmi
# [1] 0.8507881
#
# $lce$ari
# [1] 0.6882715
#
#
# $lca
# $lca$confmat
# # A tibble: 13 × 3
# .metric              .estimator .estimate
# <chr>                <chr>          <dbl>
#   1 accuracy             multiclass     0.805
# 2 kap                  multiclass     0.753
# 3 sens                 macro          0.829
# 4 spec                 macro          0.949
# 5 ppv                  macro          0.837
# 6 npv                  macro          0.951
# 7 mcc                  multiclass     0.759
# 8 j_index              macro          0.778
# 9 bal_accuracy         macro          0.889
# 10 detection_prevalence macro          0.2
# 11 precision            macro          0.837
# 12 recall               macro          0.829
# 13 f_meas               macro          0.825
#
# $lca$nmi
# [1] 0.7791981
#
# $lca$ari
# [1] 0.6592031
#
# "The Inf-FS2020 method selected 91 features."

############################
# ALLAML  TBC
############################
data(ALLAML)

data <- ALLAML

mean(data$x)
sd(data$x)

############################
# leukemia TBC
############################


############################
# lung_cancer
############################
data(lung_cancer)

data <- lung_cancer

mean(data$x)
sd(data$x)

# Run just once
UFS_results_file <- "experiments/UFS_lung_cancer_Inf-FS2020.RData"
# UFS_results <- runUFS(data$x, UFS_Methods)
# save(UFS_results, file = UFS_results_file)

load(UFS_results_file)

method = "Inf-FS2020"
top_features <- UFS_results$Results[[method]]$Result[[3]]
sprintf("The %s method selected %d features.", method, length(top_features))

g <- 5
B <- 100
execution_time <- system.time(RPClu_results <- RPClu_parallel(data$x[, top_features],
                                                              clust_fun = "gmm",
                                                              g=g,
                                                              B=B,
                                                              verb=TRUE))["elapsed"]

E <- as.matrix(RPClu_results$clusterings)

external_metrics <- calculate_external_metrics(E,g, as.integer(data$y))

external_metrics

# > external_metrics
# $kmodes
# $kmodes$confmat
# # A tibble: 13 × 3
# .metric              .estimator .estimate
# <chr>                <chr>          <dbl>
#   1 accuracy             multiclass     0.690
# 2 kap                  multiclass     0.497
# 3 sens                 macro          0.634
# 4 spec                 macro          0.912
# 5 ppv                  macro          0.646
# 6 npv                  macro          0.888
# 7 mcc                  multiclass     0.527
# 8 j_index              macro          0.546
# 9 bal_accuracy         macro          0.773
# 10 detection_prevalence macro          0.2
# 11 precision            macro          0.646
# 12 recall               macro          0.634
# 13 f_meas               macro          0.627
#
# $kmodes$nmi
# [1] 0.4651067
#
# $kmodes$ari
# [1] 0.3685128
#
#
# $majority
# $majority$confmat
# # A tibble: 13 × 3
# .metric              .estimator .estimate
# <chr>                <chr>          <dbl>
#   1 accuracy             multiclass     0.704
# 2 kap                  multiclass     0.499
# 3 sens                 macro          0.621
# 4 spec                 macro          0.911
# 5 ppv                  macro          0.696
# 6 npv                  macro          0.890
# 7 mcc                  multiclass     0.524
# 8 j_index              macro          0.532
# 9 bal_accuracy         macro          0.766
# 10 detection_prevalence macro          0.2
# 11 precision            macro          0.696
# 12 recall               macro          0.621
# 13 f_meas               macro          0.634
#
# $majority$nmi
# [1] 0.4598495
#
# $majority$ari
# [1] 0.4061217
#
#
# $lce
# $lce$confmat
# # A tibble: 13 × 3
# .metric              .estimator .estimate
# <chr>                <chr>          <dbl>
#   1 accuracy             multiclass     0.700
# 2 kap                  multiclass     0.483
# 3 sens                 macro          0.577
# 4 spec                 macro          0.912
# 5 ppv                  macro          0.588
# 6 npv                  macro          0.890
# 7 mcc                  multiclass     0.504
# 8 j_index              macro          0.489
# 9 bal_accuracy         macro          0.744
# 10 detection_prevalence macro          0.2
# 11 precision            macro          0.588
# 12 recall               macro          0.577
# 13 f_meas               macro          0.560
#
# $lce$nmi
# [1] 0.4692405
#
# $lce$ari
# [1] 0.4093116
#
#
# $lca
# $lca$confmat
# # A tibble: 13 × 3
# .metric              .estimator .estimate
# <chr>                <chr>          <dbl>
#   1 accuracy             multiclass     0.611
# 2 kap                  multiclass     0.429
# 3 sens                 macro          0.657
# 4 spec                 macro          0.902
# 5 ppv                  macro          0.579
# 6 npv                  macro          0.877
# 7 mcc                  multiclass     0.480
# 8 j_index              macro          0.559
# 9 bal_accuracy         macro          0.779
# 10 detection_prevalence macro          0.2
# 11 precision            macro          0.579
# 12 recall               macro          0.657
# 13 f_meas               macro          0.579
#
# $lca$nmi
# [1] 0.4035596
#
# $lca$ari
# [1] 0.2510463
#
# [1] "The Inf-FS2020 method selected 583 features."

############################
# lymphoma TBC
############################
data(lymphoma)

data <- lymphoma

mean(data$x)
sd(data$x)

############################
# prostate_ge
############################

data(prostate_ge)

data <- prostate_ge

mean(data$x)
sd(data$x)

# Run just once
UFS_results_file <- "experiments/UFS_prostate_ge_Inf-FS2020.RData"
# UFS_results <- runUFS(data$x, UFS_Methods)
# save(UFS_results, file = UFS_results_file)

load(UFS_results_file)

method = "Inf-FS2020"
top_features <- UFS_results$Results[[method]]$Result[[3]]
sprintf("The %s method selected %d features.", method, length(top_features))

g <- 2
B <- 100
execution_time <- system.time(RPClu_results <- RPClu_parallel(data$x[, top_features],
                                                              clust_fun = "km",
                                                              g=g,
                                                              B=B,
                                                              verb=TRUE))["elapsed"]

E <- as.matrix(RPClu_results$clusterings)

external_metrics <- calculate_external_metrics(E,g, as.integer(data$y))

external_metrics

# > external_metrics
# $kmodes
# $kmodes$confmat
# # A tibble: 13 × 3
# .metric              .estimator .estimate
# <chr>                <chr>          <dbl>
#   1 accuracy             binary         0.716
# 2 kap                  binary         0.429
# 3 sens                 binary         0.62
# 4 spec                 binary         0.808
# 5 ppv                  binary         0.756
# 6 npv                  binary         0.689
# 7 mcc                  binary         0.436
# 8 j_index              binary         0.428
# 9 bal_accuracy         binary         0.714
# 10 detection_prevalence binary         0.402
# 11 precision            binary         0.756
# 12 recall               binary         0.62
# 13 f_meas               binary         0.681
#
# $kmodes$nmi
# [1] 0.1444463
#
# $kmodes$ari
# [1] 0.1782498
#
#
# $majority
# $majority$confmat
# # A tibble: 13 × 3
# .metric              .estimator .estimate
# <chr>                <chr>          <dbl>
#   1 accuracy             binary         0.569
# 2 kap                  binary         0.127
# 3 sens                 binary         0.26
# 4 spec                 binary         0.865
# 5 ppv                  binary         0.65
# 6 npv                  binary         0.549
# 7 mcc                  binary         0.158
# 8 j_index              binary         0.125
# 9 bal_accuracy         binary         0.563
# 10 detection_prevalence binary         0.196
# 11 precision            binary         0.65
# 12 recall               binary         0.26
# 13 f_meas               binary         0.371
#
# $majority$nmi
# [1] 0.02151663
#
# $majority$ari
# [1] 0.01253643
#
#
# $lce
# $lce$confmat
# # A tibble: 13 × 3
# .metric              .estimator .estimate
# <chr>                <chr>          <dbl>
#   1 accuracy             binary        0.569
# 2 kap                  binary        0.123
# 3 sens                 binary        0.14
# 4 spec                 binary        0.981
# 5 ppv                  binary        0.875
# 6 npv                  binary        0.543
# 7 mcc                  binary        0.225
# 8 j_index              binary        0.121
# 9 bal_accuracy         binary        0.560
# 10 detection_prevalence binary        0.0784
# 11 precision            binary        0.875
# 12 recall               binary        0.14
# 13 f_meas               binary        0.241
#
# $lce$nmi
# [1] 0.06406633
#
# $lce$ari
# [1] 0.01575351
#
#
# $lca
# $lca$confmat
# # A tibble: 13 × 3
# .metric              .estimator .estimate
# <chr>                <chr>          <dbl>
#   1 accuracy             binary         0.755
# 2 kap                  binary         0.509
# 3 sens                 binary         0.72
# 4 spec                 binary         0.788
# 5 ppv                  binary         0.766
# 6 npv                  binary         0.745
# 7 mcc                  binary         0.510
# 8 j_index              binary         0.508
# 9 bal_accuracy         binary         0.754
# 10 detection_prevalence binary         0.461
# 11 precision            binary         0.766
# 12 recall               binary         0.72
# 13 f_meas               binary         0.742
#
# $lca$nmi
# [1] 0.1971804
#
# $lca$ari
# [1] 0.2525461
#

############################
# COIL20
############################
data(COIL20)

data <- COIL20

data$x <- scale(data$x)

mean(data$x)
sd(data$x)

# Run just once
UFS_results_file <- "experiments/UFS_COIL20_Inf-FS2020.RData"
UFS_results <- runUFS(data$x, UFS_Methods)
save(UFS_results, file = UFS_results_file)

load(UFS_results_file)

method = "Inf-FS2020"
top_features <- UFS_results$Results[[method]]$Result[[3]]
sprintf("The %s method selected %d features.", method, length(top_features))

g <- 20
B <- 100
execution_time <- system.time(RPClu_results <- RPClu_parallel(data$x[, top_features],
                                                              clust_fun = "km",
                                                              g=g,
                                                              B=B,
                                                              verb=TRUE))["elapsed"]

E <- as.matrix(RPClu_results$clusterings)

external_metrics <- calculate_external_metrics(E,g, as.integer(data$y))

external_metrics

# [1] "The Inf-FS2020 method selected 175 features."
#
# $kmodes
# $kmodes$confmat
# # A tibble: 13 × 3
# .metric              .estimator .estimate
# <chr>                <chr>          <dbl>
#   1 accuracy             multiclass     0.569
# 2 kap                  multiclass     0.546
# 3 sens                 macro          0.569
# 4 spec                 macro          0.977
# 5 ppv                  macro          0.553
# 6 npv                  macro          0.978
# 7 mcc                  multiclass     0.550
# 8 j_index              macro          0.546
# 9 bal_accuracy         macro          0.773
# 10 detection_prevalence macro          0.05
# 11 precision            macro          0.553
# 12 recall               macro          0.569
# 13 f_meas               macro          0.537
#
# $kmodes$nmi
# [1] 0.6955416
#
# $kmodes$ari
# [1] 0.4702302
#
#
# $majority
# $majority$confmat
# # A tibble: 13 × 3
# .metric              .estimator .estimate
# <chr>                <chr>          <dbl>
#   1 accuracy             multiclass     0.580
# 2 kap                  multiclass     0.558
# 3 sens                 macro          0.580
# 4 spec                 macro          0.978
# 5 ppv                  macro          0.583
# 6 npv                  macro          0.978
# 7 mcc                  multiclass     0.562
# 8 j_index              macro          0.558
# 9 bal_accuracy         macro          0.779
# 10 detection_prevalence macro          0.05
# 11 precision            macro          0.583
# 12 recall               macro          0.580
# 13 f_meas               macro          0.549
#
# $majority$nmi
# [1] 0.7190488
#
# $majority$ari
# [1] 0.4917776
############################
# warpAR10P
############################

data(warpAR10P)

data <- warpAR10P

data$x <- scale(data$x)

mean(data$x)
sd(data$x)

# Run just once
UFS_results_file <- "experiments/UFS_warpAR10P_Inf-FS2020.RData"
UFS_results <- runUFS(data$x, UFS_Methods)
save(UFS_results, file = UFS_results_file)

load(UFS_results_file)

method = "Inf-FS2020"
top_features <- UFS_results$Results[[method]]$Result[[3]]
sprintf("The %s method selected %d features.", method, length(top_features))

g <- 10
B <- 100
execution_time <- system.time(RPClu_results <- RPClu_parallel(data$x[, top_features],
                                                              clust_fun = "km",
                                                              g=g,
                                                              B=B,
                                                              verb=TRUE))["elapsed"]

E <- as.matrix(RPClu_results$clusterings)

external_metrics <- calculate_external_metrics(E,g, as.integer(data$y))

external_metrics


############################
# warpPIE10P
############################

data(warpPIE10P)

data <- warpPIE10P

data$x <- apply(data$x, 2, rescale)

mean(data$x)
sd(data$x)
min(data$x)
max(data$x)


# Run just once
UFS_results_file <- "experiments/UFS_warpPIE10P_Inf-FS2020.RData"
# UFS_results <- runUFS(data$x, UFS_Methods)
# save(UFS_results, file = UFS_results_file)
#
load(UFS_results_file)

g <- 10
B <- 20

top_features <- UFS_results$Results[[method]]$Result[[3]]
sprintf("The %s method selected %d features.", method, length(top_features))

execution_time <- system.time(RPClu_results <- RPClu_parallel(data$x[, top_features],
                                                              clust_fun = "gmm",
                                                              g=g,
                                                              B=B,
                                                              verb=TRUE))["elapsed"]

E <- as.matrix(RPClu_results$clusterings)

external_metrics <- calculate_external_metrics(E,g, as.integer(data$y))

print(external_metrics)

# [1] "The Inf-FS2020 method selected 382 features."
#
#
# $kmodes
# $kmodes$confmat
# # A tibble: 13 × 3
# .metric              .estimator .estimate
# <chr>                <chr>          <dbl>
#   1 accuracy             multiclass     0.514
# 2 kap                  multiclass     0.460
# 3 sens                 macro          0.514
# 4 spec                 macro          0.946
# 5 ppv                  macro          0.562
# 6 npv                  macro          0.947
# 7 mcc                  multiclass     0.468
# 8 j_index              macro          0.460
# 9 bal_accuracy         macro          0.730
# 10 detection_prevalence macro          0.1
# 11 precision            macro          0.562
# 12 recall               macro          0.514
# 13 f_meas               macro          0.510
#
# $kmodes$nmi
# [1] 0.6865815
#
# $kmodes$ari
# [1] 0.4246193
#
#
# $majority
# $majority$confmat
# # A tibble: 13 × 3
# .metric              .estimator .estimate
# <chr>                <chr>          <dbl>
#   1 accuracy             multiclass     0.633
# 2 kap                  multiclass     0.593
# 3 sens                 macro          0.633
# 4 spec                 macro          0.959
# 5 ppv                  macro          0.679
# 6 npv                  macro          0.960
# 7 mcc                  multiclass     0.599
# 8 j_index              macro          0.593
# 9 bal_accuracy         macro          0.796
# 10 detection_prevalence macro          0.1
# 11 precision            macro          0.679
# 12 recall               macro          0.633
# 13 f_meas               macro          0.633
#
# $majority$nmi
# [1] 0.7351729
#
# $majority$ari
# [1] 0.5161567
#
#
# $lce
# $lce$confmat
# # A tibble: 13 × 3
# .metric              .estimator .estimate
# <chr>                <chr>          <dbl>
#   1 accuracy             multiclass     0.714
# 2 kap                  multiclass     0.683
# 3 sens                 macro          0.714
# 4 spec                 macro          0.968
# 5 ppv                  macro          0.796
# 6 npv                  macro          0.968
# 7 mcc                  multiclass     0.688
# 8 j_index              macro          0.683
# 9 bal_accuracy         macro          0.841
# 10 detection_prevalence macro          0.1
# 11 precision            macro          0.796
# 12 recall               macro          0.714
# 13 f_meas               macro          0.735
#
# $lce$nmi
# [1] 0.7803953
#
# $lce$ari
# [1] 0.5842485
#
#
# $lca
# $lca$confmat
# # A tibble: 13 × 3
# .metric              .estimator .estimate
# <chr>                <chr>          <dbl>
#   1 accuracy             multiclass     0.6
# 2 kap                  multiclass     0.556
# 3 sens                 macro          0.6
# 4 spec                 macro          0.956
# 5 ppv                  macro          0.668
# 6 npv                  macro          0.956
# 7 mcc                  multiclass     0.560
# 8 j_index              macro          0.556
# 9 bal_accuracy         macro          0.778
# 10 detection_prevalence macro          0.1
# 11 precision            macro          0.668
# 12 recall               macro          0.6
# 13 f_meas               macro          0.609
#
# $lca$nmi
# [1] 0.7288964
#
# $lca$ari
# [1] 0.4899206

