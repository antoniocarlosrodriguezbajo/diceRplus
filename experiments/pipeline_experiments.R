library(diceRplus)
library(mclust)
library(R.matlab)

algorithms = c("nmf", "hc", "diana", "km", "pam", "ap", "sc",
               "gmm", "block", "som", "cmeans")

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
calculate_internal_metrics <- function(data,cluster_labels) {
  num_labels <- length(cluster_labels)

  # Create structure for consensus_evaluate
  cc_data = array(cluster_labels, dim = c(num_labels, 1, 1, 1))
  row_names = rownames(data)
  dimnames(cc_data) = list(
    row_names,  # Primer nivel de nombres: nombres de las filas de Meat$x
    "R1",       # Repetition (just a placeholder)
    "Algo", # Clustering algorithm (just a placeholder)
    "99"         # Number of clusters (just a placeholder)
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

  # Optionally, return as a data frame for easy inspection
  metrics_df <- do.call(rbind, lapply(names(eval_results), function(id) {
    row <- eval_results[[id]]
    data.frame(
      experiment_id = id,
      nmi = row$nmi,
      ari = row$ari,
      stringsAsFactors = FALSE
    )
  }))

  return(metrics_df)
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

    # Evaluate internal metrics for all combinations of this clustering method
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
                          evaluate = TRUE,
                          seed = seed,
                          seed.data = seed,
                          progress = TRUE,
                          verbose = FALSE)
    })["elapsed"]

    for (consensus_method in colnames(results_dice$clusters)) {
      consensus_labels = as.vector(results_dice$clusters[ ,consensus_method])

      # Calculate internal metrics
      internal_metrics = calculate_internal_metrics(data,consensus_labels)

      # Store the experiment
      description = paste0("Ensemble Clustering: ", algo, " + ", UFS_method , " + ",
                           consensus_method, " + " ,
                           " seed ", as.character(seed))
      exp_data = experiment_logger(
        description = description,
        dataset = dataset,
        clustering_method = algo,
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


run_pipeline <- function(data,
                         dataset_name,
                         top_features,
                         nk,
                         algorithms,
                         UFS_method = "Inf-FS2020",
                         alpha = 0.1,
                         n_reps=20,
                         seed=100) {
  # Best clustering with UFS
  print("Best clustering with UFS")
  experiment_ids <- run_clustering_experiments(data = data,
                                               dataset = dataset_name,
                                               top_features = top_features,
                                               nk = nk,
                                               algorithms = algorithms,
                                               UFS_method = "Inf-FS2020",
                                               alpha = 0.1,
                                               seed = seed)

  # Get best ensemble (clustering + consensus) according to internal metrics
  result <- get_best_ensemble_internal_metrics(experiment_ids = experiment_ids,
                                               data = data[,top_features],
                                               nk = nk)

  best_experiment_id <- result$best_experiment_id
  best_overall <- result$trim.obj$alg.keep

  parts <- strsplit(best_overall, "\\+")[[1]]
  best_clustering <- parts[1]
  best_consensus <- parts[2]

  print(sprintf("Best Clustering: %s, Best Consensus: %s, id %d",
                best_clustering, best_consensus,best_experiment_id))
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

run_pipeline(data = data$x,
             dataset_name = "Meat",
             top_features = results_UFS$top_features,
             nk = 5,
             algorithms = algorithms,
             UFS_method = "Inf-FS2020",
             alpha = results_UFS$alpha_best,
             n_reps=20,
             seed=100)

calculate_external_metrics_from_ids(2174, as.integer(data$y))
