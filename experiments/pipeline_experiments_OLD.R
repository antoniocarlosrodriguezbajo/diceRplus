library(diceRplus)
library(mclust)
library(R.matlab)
library(scales)

algorithms = c("km", "gmm", "sc", "cmeans", "pam")


mean_abs_correlation <- function(data) {
  cor_mat <- cor(data, use = "pairwise.complete.obs")
  round(mean(abs(cor_mat[upper.tri(data)])),4)
}

get_best_top_features <- function(data) {
  # Start MATLAB server
  Matlab$startServer()
  matlab <- Matlab()
  print(matlab)

  # Check if MATLAB server is running
  isOpen <- open(matlab)
  if (!isOpen) {
    print("MATLAB server is not running: waited 30 seconds.")
  }
  print(matlab)

  # Set the dataset variable
  setVariable(matlab, X = data)

  # Evaluate MATLAB command to display variables
  evaluate(matlab, "whos")

  evaluate(matlab, "rng('default');")
  evaluate(matlab, "rng(42);")

  num_top_features = 1e10

  for (alpha in seq(0, 0.5, by = 0.1)) {

    evaluate(matlab, "clear RANKED WEIGHT SUBSET WEIGHT_SUM WEIGHT_MEAN;")


    # Run the feature selection method in MATLAB
    cmd <- sprintf('[RANKED, WEIGHT, SUBSET] = InfFS_U(X, %.1f);', alpha)
    evaluate(matlab, cmd)

    ranked  = getVariable(matlab, "RANKED")
    weight  = getVariable(matlab, "WEIGHT")
    subset  = getVariable(matlab, "SUBSET")

    top_features_alpha <- as.vector(subset$SUBSET)
    num_features_alpha = length(top_features_alpha)
    print(alpha)
    print(num_features_alpha)

    if (num_features_alpha < num_top_features) {
      alpha_best = alpha
      top_features = top_features_alpha
      num_top_features = num_features_alpha
    }
  }

  close(matlab)

  return(list(alpha_best=alpha_best,
              top_features=top_features,
              num_features=num_top_features))
}

# Calcute internal_metrics based on consensus_evaluate of diceR
calculate_internal_metrics <- function(data,cluster_labels,nk) {
  num_labels <- length(cluster_labels)

  # Create structure for consensus_evaluate
  cc_data = array(cluster_labels, dim = c(num_labels, 1, 1, 1))
  row_names = rownames(data)
  dimnames(cc_data) = list(
    row_names,  # Primer nivel de nombres: nombres de las filas de Meat$x
    "R1",       # Repetition (just a placeholder)
    "Algo", # Clustering algorithm (just a placeholder)
    as.character(nk)
  )
  # Evaluation
  result_evaluation = consensus_evaluate(data, cc_data)
  metrics_df = result_evaluation$ii[[1]]
  metrics_list = lapply(metrics_df[1, -1], function(x) unname(x))
  return(metrics_list)
}

calculate_external_metrics <- function (consensus_labels, ref_labels) {
  # Evaluate results
  eval_results = lapply(consensus_labels, function(labels) list(
    confmat = ev_confmat(labels, ref_labels),
    nmi     = ev_nmi(labels, ref_labels),
    ari     = adjustedRandIndex(labels, ref_labels)
  ))

  return(eval_results)
}

calculate_external_metrics_from_ids <- function(experiment_ids, ref_labels) {
  # Load all experiment data
  experiments_data_all <- load_experiments()

  # Filter only the requested experiments
  filtered_df <- experiments_data_all[
    experiments_data_all$experiment_id %in% experiment_ids,
  ]

  # Prepare a named list of clustering labels (for each experiment)
  consensus_labels <- setNames(
    lapply(filtered_df$labels_clustering, as.integer),
    filtered_df$experiment_id
  )

  # Compute external metrics
  eval_results <- calculate_external_metrics(consensus_labels, ref_labels)

  metrics_df <- do.call(rbind, lapply(names(eval_results), function(id) {
    row <- eval_results[[id]]
    data.frame(
      experiment_id = id,
      nmi = row$nmi,
      ari = row$ari,
      acc = row$confmat$.estimate[[1]],
      stringsAsFactors = FALSE
    )
  }))

  return(metrics_df)
}

evaluate_experiment_clusterings <- function(experiment_ids, data, nk, seed = 101) {
  set.seed(seed)

  # Load all experiment metadata (assumed to include labels and experiment IDs)
  experiments_data_all <- load_experiments()

  # Filter to include only the specified experiment IDs
  filtered_df <- experiments_data_all[
    experiments_data_all$experiment_id %in% experiment_ids,
  ]

  clustering_list <- list()

  # Prepare each experiment's labels in the format expected by consensus_evaluate
  for (i in seq_len(nrow(filtered_df))) {
    labels <- as.integer(filtered_df$labels_clustering[[i]])
    name <- as.character(filtered_df$experiment_id[[i]])
    num_elements <- length(labels)

    # Create a 4D array with appropriate dimension names
    clustering_array <- array(
      labels,
      dim = c(num_elements, 1, 1, 1),
      dimnames = list(NULL, "R1", name, as.character(nk))
    )

    clustering_list[[name]] <- clustering_array
  }

  # Run internal validation of all selected clusterings
  result <- do.call(consensus_evaluate, c(
    list(data), clustering_list,
    list(n = 1, k.method = nk, trim = TRUE, reweigh = TRUE)
  ))

  return(result)
}



get_best_ensemble_internal_metrics <- function(experiment_ids, data, nk, seed = 101) {
  # Load all experiment metadata
  experiments_data_all <- load_experiments()

  # Filter dataset to only include the experiments of interest
  filtered_df <- experiments_data_all[
    experiments_data_all$experiment_id %in% experiment_ids,
  ]

  # Create a unique identifier for each clustering + consensus combination
  filtered_df$combo_method <- paste0(
    filtered_df$clustering_method, "+", filtered_df$consensus_method
  )

  # Get all distinct clustering methods
  clustering_methods <- unique(filtered_df$clustering_method)

  # Initialize objects to store best combinations and their corresponding experiment IDs
  best_combos <- list()
  combo_to_experiment_id <- list()

  set.seed(seed)

  # Iterate over each clustering method
  for (cm in clustering_methods) {
    # Subset experiments for the current clustering method
    subset_df <- filtered_df[filtered_df$clustering_method == cm, ]
    cc_data_list <- list()

    # Convert each experiment's labels into the format expected by consensus_evaluate
    for (i in seq_len(nrow(subset_df))) {
      combo <- paste0(subset_df$clustering_method[i], "+", subset_df$consensus_method[i])
      labels <- as.integer(subset_df$labels_clustering[[i]])
      num_elements <- length(labels)

      cc_data <- array(
        labels,
        dim = c(num_elements, 1, 1, 1),
        dimnames = list(NULL, "R1", combo, as.character(nk))
      )

      cc_data_list[[combo]] <- cc_data
    }

    # Evaluate internal metrics of this clustering method
    result <- do.call(consensus_evaluate, c(
      list(data), cc_data_list,
      list(n = 1, k.method = nk, trim = TRUE, reweigh = TRUE)
    ))

    # Get the best-performing combination
    best_combo <- result$trim.obj$alg.keep
    best_combos[[best_combo]] <- cc_data_list[[best_combo]]

    # Record the corresponding experiment_id for this best combo
    best_id <- subset_df$experiment_id[subset_df$combo_method == best_combo][1]
    combo_to_experiment_id[[best_combo]] <- best_id
  }

  # Run a final evaluation on the best combinations
  final_result <- do.call(consensus_evaluate, c(
    list(data), best_combos,
    list(n = 1, k.method = nk, trim = TRUE, reweigh = TRUE)
  ))

  # Identify the final winning combination and its experiment ID
  final_combo <- final_result$trim.obj$alg.keep
  final_experiment_id <- combo_to_experiment_id[[final_combo]]

  # Attach that experiment ID to the final result
  final_result$best_experiment_id <- final_experiment_id

  return(final_result)
}


run_clustering_experiments <- function(data_all = data,
                                       dataset = dataset_name,
                                       top_features = top_features,
                                       nk = nk,
                                       algorithms = algorithms,
                                       cons.funs = c("kmodes", "majority", "CSPA", "LCE", "LCA"),
                                       UFS_method = "Inf-FS2020",
                                       alpha = 0.1,
                                       seed = 101) {
  if (is.null(top_features)) {
    data = data_all
  } else {
    data = data_all[,top_features]
  }
  set.seed(seed)

  experiment_ids = c()

  for (algo in algorithms) {
    print(algo)
    execution_time = system.time({
      results_dice = dice(data = data,
                          k.method = nk,
                          nk = nk,
                          algorithms = algo,
                          cons.funs = cons.funs,
                          evaluate = TRUE,
                          seed = seed,
                          seed.data = seed,
                          progress = FALSE,
                          verbose = FALSE)
    })["elapsed"]

    for (consensus_method in colnames(results_dice$clusters)) {
      consensus_labels = as.vector(results_dice$clusters[ ,consensus_method])

      # Calculate internal metrics
      internal_metrics = calculate_internal_metrics(data,consensus_labels,nk)

      # Store the experiment
      description = paste0("Ensemble Clustering: ", algo, " + ", UFS_method , " + ",
                           consensus_method, " + " ,
                           " seed ", as.character(seed))
      exp_data = experiment_logger(
        description = description,
        dataset = dataset,
        clustering_method = algo,
        clustering_method_params = list(nk = nk),
        UFS_method = UFS_method,
        UFS_method_params = list(alpha = alpha),
        num_features = ncol(data),
        features = top_features,
        consensus_method = consensus_method,
        execution_time = as.numeric(execution_time),
        labels_clustering = consensus_labels,
        internal_metrics = internal_metrics
      )
      save_experiment(exp_data)
      experiment_ids = c(experiment_ids, exp_data$experiment_id)
    }

  }
  return(experiment_ids)
}


run_RPClu_experiments <- function(data_all = data,
                                  dataset = dataset_name,
                                  top_features = top_features,
                                  nk = nk,
                                  B = 100,
                                  B.star = 10,
                                  UFS_method = "Inf-FS2020",
                                  alpha = 0.1,
                                  seed = 101) {
  if (is.null(top_features)) {
    data = data_all
  } else {
    data = data_all[,top_features]
  }
  execution_time <- system.time({
    RPClu_results <- RPGMMClu_parallel(data,
                                       g = nk,
                                       B = B,
                                       B.star = B.star,
                                       seed = seed)
  })["elapsed"]

  labels_clustering = RPClu_results$ensemble$label.vec

  # Store the experiment
  description = paste0("Ensemble Clustering: ", "RPClu", " + ", UFS_method , " + ",
                       "DWH", " + " ,
                       " seed ", as.character(seed))

  exp_data = experiment_logger(
    description = description,
    dataset = dataset,
    clustering_method = "RPClu",
    clustering_method_params = list(nk = nk, B = B),
    UFS_method = UFS_method,
    UFS_method_params = list(alpha = alpha),
    num_features = ncol(data),
    features = top_features,
    consensus_method = "DWH",
    consensus_method_params = list(nk = nk, B = B, B.star = B.star),
    execution_time = as.numeric(execution_time),
    labels_clustering = labels_clustering,
    internal_metrics = calculate_internal_metrics(data,labels_clustering,nk),
  )
  save_experiment(exp_data)

  return(exp_data$experiment_id)
}



run_pipeline <- function(data,
                             dataset_name,
                             top_features,
                             nk,
                             algorithms,
                             UFS_method = "Inf-FS2020",
                             alpha = 0.1,
                             n_reps=1,
                             seed=100) {
  # Best clustering with UFS
  print("Best clustering with UFS")
  experiment_UFS_ids <- run_clustering_experiments(data_all = data,
                                                   dataset = dataset_name,
                                                   top_features = top_features,
                                                   nk = nk,
                                                   algorithms = algorithms,
                                                   UFS_method = "Inf-FS2020",
                                                   alpha = alpha,
                                                   seed = seed)

  # Get best ensemble (clustering + consensus) according to internal metrics
  result <- get_best_ensemble_internal_metrics(experiment_ids = experiment_UFS_ids,
                                               data = data[,top_features],
                                               nk = nk)

  best_experiment_UFS_id <- result$best_experiment_id
  best_overall <- result$trim.obj$alg.keep

  parts <- strsplit(best_overall, "\\+")[[1]]
  best_clustering <- parts[1]
  best_consensus <- parts[2]

  print(sprintf("Best Clustering UFS: %s, Best Consensus: %s, Experiment id: %d",
                best_clustering, best_consensus,best_experiment_UFS_id))

  # Best clustering without UFS
  print("Best clustering without UFS")
  experiment_NO_UFS_ids <- run_clustering_experiments(data_all = data,
                                                   dataset = dataset_name,
                                                   top_features = NULL,
                                                   nk = nk,
                                                   algorithms = algorithms,
                                                   UFS_method = "NO UFS",
                                                   seed = seed)

  # Get best ensemble (clustering + consensus) according to internal metrics
  result <- get_best_ensemble_internal_metrics(experiment_ids = experiment_NO_UFS_ids,
                                               data = data[,top_features],
                                               nk = nk)

  best_experiment_NO_UFS_id <- result$best_experiment_id
  best_overall <- result$trim.obj$alg.keep

  parts <- strsplit(best_overall, "\\+")[[1]]
  best_clustering <- parts[1]
  best_consensus <- parts[2]

  print(sprintf("Best Clustering NO UFS: %s, Best Consensus: %s, Experiment id: %d",
                best_clustering, best_consensus,best_experiment_NO_UFS_id))


  # RPClu with UFS
  print("RPClu with UFS")
  experiment_id <- run_RPClu_experiments(data_all = data,
                                         dataset = dataset_name,
                                         top_features = top_features,
                                         nk = nk,
                                         B = 100,
                                         B.star = 10,
                                         UFS_method = "Inf-FS2020",
                                         alpha = 0.1,
                                         seed = 101)
  print(experiment_id)

  return(list(best_experiment_UFS_id = best_experiment_UFS_id,
              best_experiment_NO_UFS_id = best_experiment_NO_UFS_id,
              best_experiment_RPClu = experiment_id
              ))
}




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
             result$best_experiment_RPClu)

em <- calculate_external_metrics_from_ids(exp_ids, as.integer(data$y))
em$ari


result <- evaluate_experiment_clusterings(
  experiment_ids = c(265,287,298),
  data = data$x,
  nk = 5,
  seed = 101
)

# "Best Clustering UFS: cmeans, Best Consensus: CSPA, Experiment id: 265"
# "Best Clustering NO UFS: sc, Best Consensus: LCA, Experiment id: 287"
# [1] "RPClu with UFS"
# [1] 298
# > em$ari
# [1] 0.3298860 0.2425831 0.6897873

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
             result$best_experiment_RPClu)

em <- calculate_external_metrics_from_ids(exp_ids, as.integer(data$y))
em$ari

# "Best Clustering UFS: sc, Best Consensus: LCE, Experiment id: 312"
# "Best Clustering NO UFS: cmeans, Best Consensus: CSPA, Experiment id: 341"
# [1] "RPClu with UFS"
# [1] 349
# > em$ari
# [1] 0.7884032 0.8518607 1.0000000


########################################################
########################################################
# lung_cancer
########################################################
########################################################

data(lung_cancer)

data <- lung_cancer

mean(data$x)
sd(data$x)

mean_abs_correlation(data$x)
#0.1467

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
             result$best_experiment_RPClu)
em <- calculate_external_metrics_from_ids(exp_ids, as.integer(data$y))
em$ari

# Best Clustering UFS: cmeans, Best Consensus: LCE, Experiment id: 368"
# "Best Clustering NO UFS: pam, Best Consensus: CSPA, Experiment id: 397"
# [1] "RPClu with UFS"
# [1] 400
# > em$ari
# [1] 0.2092839 0.4814490 0.4013991

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
             algorithms = setdiff(algorithms, "ap"),
             UFS_method = "Inf-FS2020",
             alpha = results_UFS$alpha_best,
             n_reps=1,
             seed=100)

exp_ids <- c(result$best_experiment_UFS_id,
             result$best_experiment_NO_UFS_id,
             result$best_experiment_RPClu)
em <- calculate_external_metrics_from_ids(exp_ids, as.integer(data$y))
em$ari

# "Best Clustering UFS: cmeans, Best Consensus: kmodes, Experiment id: 416"
# "Best Clustering NO UFS: sc, Best Consensus: majority, Experiment id: 437"
# [1] "RPClu with UFS"
# [1] 451
# > em$ari
# [1] 0.1597852 0.1071034 0.5379900



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
# 0.1643

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
             result$best_experiment_RPClu)
em <- calculate_external_metrics_from_ids(exp_ids, as.integer(data$y))
em$ari

"Best Clustering UFS: gmm, Best Consensus: kmodes, Experiment id: 457"
"Best Clustering NO UFS: km, Best Consensus: LCE, Experiment id: 480"
[1] "RPClu with UFS"
[1] 502

########################################################
########################################################
# leukemia
########################################################
########################################################

data(leukemia)

data <- leukemia

mean(data$x)
sd(data$x)
min(data$x)
max(data$x)

mean_abs_correlation(data$x)
# 0.1508

data$y <- ifelse(data$y == -1, 1, 2)


# UFS
results_UFS <- get_best_top_features(data$x)

result <- run_pipeline(data = data$x,
             dataset_name = "Leukemia",
             top_features = results_UFS$top_features,
             nk = 2,
             algorithms = algorithms,
             UFS_method = "Inf-FS2020",
             alpha = results_UFS$alpha_best,
             n_reps=1,
             seed=100)

exp_ids <- c(result$best_experiment_UFS_id,
             result$best_experiment_NO_UFS_id,
             result$best_experiment_RPClu)
em <- calculate_external_metrics_from_ids(exp_ids, as.integer(data$y))
em$acc
