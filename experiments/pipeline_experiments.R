library(diceRplus)
library(R.matlab)

algorithms = c("nmf", "hc", "diana", "km", "pam", "ap", "sc",
               "gmm", "block", "som", "cmeans", "hdbscan")

algorithms = c("gmm","cmeans")


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

run_clustering_experiments <- function(data_all = data,
                                       dataset = dataset_name,
                                       top_features = top_features,
                                       nk = nk,
                                       algorithms = algorithms,
                                       UFS_method = "Inf-FS2020",
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
      description = paste0("Ensemble Clustering + ", algo, " + ",
                           consensus_method, " + " ,
                           " seed ", as.character(seed))
      exp_data = experiment_logger(
        description = description,
        dataset = dataset,
        clustering_method = algo,
        clustering_method_params = list(nk = nk),
        UFS_method = UFS_method,
        num_features = ncol(data),
        features = top_features,
        consensus_method = consensus_method,
        execution_time = as.numeric(execution_time),
        labels_clustering = consensus_labels,
        internal_metrics = internal_metrics
      )
      save_experiment(exp_data)
    }



  }



  return(NULL)
}


run_pipeline <- function(data,
                         dataset_name,
                         top_features,
                         nk,
                         algorithms,
                         n_reps=20,
                         seed=100) {
  # Best clustering with UFS
  experiments_id <- run_clustering_experiments(data = data,
                                               dataset = dataset_name,
                                               top_features = top_features,
                                               nk = nk,
                                               algorithms = algorithms,
                                               UFS_method = "Inf-FS2020",
                                               seed = seed)


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
             n_reps=20,
             seed=100)


