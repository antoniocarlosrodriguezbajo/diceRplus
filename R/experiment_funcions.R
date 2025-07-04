#' Mean Absolute Pairwise Correlation
#'
#' Computes the mean absolute value of the upper-triangular elements of the pairwise correlation matrix
#' for a given numeric dataset. Useful for quantifying overall redundancy or similarity among features.
#'
#' @param data A numeric data frame or matrix. Each column is treated as a variable; each row is an observation.
#'
#' @return A single numeric value: the mean absolute correlation between variable pairs, rounded to four decimals.
#'
#' @details The function uses Pearson correlation via \code{cor()}, ignoring missing values via
#' \code{use = "pairwise.complete.obs"}. Only the upper triangle (excluding the diagonal) of the correlation matrix is considered.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(100 * 5), ncol = 5)
#' mean_abs_correlation(X)
#'
#' @export
mean_abs_correlation <- function(data) {
  cor_mat <- cor(data, use = "pairwise.complete.obs")
  round(mean(abs(cor_mat[upper.tri(data)])),4)
}


#' Select Optimal Feature Subset Using Inf-FS via MATLAB
#'
#' Connects to a MATLAB server to run the \code{InfFS\_U} algorithm over a range of \code{alpha} values and
#' identifies the subset of top features that minimizes (or maximizes) the number of selected variables.
#'
#' @param data A numeric matrix or data frame where columns are features and rows are samples.
#' @param minimum Logical. If \code{TRUE} (default), selects the smallest feature subset across the \code{alpha} range;
#' if \code{FALSE}, selects the largest subset.
#' @param range A numeric vector of length 2 specifying the lower and upper bounds of \code{alpha} to explore.
#' Values are incremented by 0.1 in the search.
#'
#' @return A list with three elements:
#' \describe{
#'   \item{\code{alpha_best}}{The alpha value yielding the optimal subset (minimum or maximum number of features).}
#'   \item{\code{top_features}}{The indices of selected features for that alpha.}
#'   \item{\code{num_features}}{The number of selected features.}
#' }
#'
#' @details This function requires a running MATLAB server and a valid installation of the Inf-FS
#' method accessible from MATLAB. Feature selection is performed in MATLAB by evaluating
#' \code{[RANKED, WEIGHT, SUBSET] = InfFS\_U(X, alpha)} at multiple \code{alpha} values.
#'
#' @examples
#' \dontrun{
#' X <- matrix(rnorm(1000), ncol = 20)
#' get_best_top_features(X, minimum = TRUE, range = c(0.0, 0.3))
#' }
#'
#' @export
get_best_top_features <- function(data, minimum=TRUE, range=c(0,0.5)) {
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

  if (minimum) {
    num_top_features = 1e10
  } else {
    num_top_features = -1e10
  }

  for (alpha in seq(range[1], range[2], by = 0.1)) {

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

    update_needed <- (minimum && num_features_alpha < num_top_features) ||
      (!minimum && num_features_alpha > num_top_features)

    if (update_needed) {
      alpha_best <- alpha
      top_features <- top_features_alpha
      num_top_features <- num_features_alpha
    }
  }

  close(matlab)

  return(list(alpha_best=alpha_best,
              top_features=top_features,
              num_features=num_top_features))
}

#' Compute Internal Clustering Metrics via diceR Consensus Evaluation
#'
#' Computes a set of internal validation metrics (e.g., connectivity, Dunn index, silhouette)
#' for a single clustering solution using \code{consensus_evaluate()} from the diceR package.
#'
#' @param data A numeric data frame or matrix of features (samples in rows, variables in columns).
#' @param cluster_labels A vector of cluster assignments for each observation.
#' @param nk An integer specifying the number of clusters used. Used as a dimension label.
#'
#' @return A named list containing internal clustering metric values as computed by diceR:
#' typically including \code{connectivity}, \code{dunn}, and \code{silhouette}.
#'
#' @details This function wraps a single clustering solution into the expected four-dimensional array
#' format required by \code{consensus_evaluate()}, assigning dummy labels for repetition and algorithm.
#'
#' @examples
#' \dontrun{
#' library(diceRplus)
#' data(iris)
#' cl_labels <- kmeans(iris[, -5], centers = 3)$cluster
#' calculate_internal_metrics(iris[, -5], cl_labels, nk = 3)
#' }
#'
#' @export
# Calculate internal_metrics based on consensus_evaluate of diceR
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

#' Compute External Clustering Metrics
#'
#' Calculates external validation scores between predicted clustering labels and reference (ground truth) labels
#' across multiple solutions. Each solution is evaluated using adjusted Rand index (ARI), normalized mutual information (NMI),
#' and the confusion matrix.
#'
#' @param consensus_labels A list of clustering label vectors (e.g., from consensus clustering runs).
#' Each element must be a vector of cluster assignments.
#' @param ref_labels A vector of true reference labels to evaluate against.
#'
#' @return A list of lists. Each sublist corresponds to one clustering solution and contains:
#' \describe{
#'   \item{\code{confmat}}{A confusion matrix comparing predicted and true labels.}
#'   \item{\code{nmi}}{Normalized mutual information between predicted and true labels.}
#'   \item{\code{ari}}{Adjusted Rand Index between predicted and true labels.}
#' }
#'
#' @details Requires \code{ev_confmat()} and \code{ev_nmi()} to be defined or imported in your package.
#' Uses \code{mclust::adjustedRandIndex()} for ARI.
#'
#' @examples
#' \dontrun{
#' library(mclust)
#' pred_list <- list(
#'   method1 = sample(1:3, 100, replace = TRUE),
#'   method2 = kmeans(iris[, -5], centers = 3)$cluster
#' )
#' calculate_external_metrics(pred_list, iris$Species)
#' }
#'
#' @export
calculate_external_metrics <- function (consensus_labels, ref_labels) {
  # Evaluate results
  eval_results = lapply(consensus_labels, function(labels) list(
    confmat = ev_confmat(labels, ref_labels),
    nmi     = ev_nmi(labels, ref_labels),
    ari     = adjustedRandIndex(labels, ref_labels)
  ))

  return(eval_results)
}

#' Compute External Clustering Metrics from Experiment IDs
#'
#' Retrieves clustering labels from stored experiment results and computes external evaluation metrics
#' (NMI, ARI, and accuracy) by comparing them with reference ground truth labels.
#'
#' @param experiment_ids A character vector of experiment IDs to evaluate.
#' @param ref_labels A vector of ground truth class labels for the observations.
#'
#' @return A data frame where each row corresponds to an experiment and includes:
#' \describe{
#'   \item{\code{experiment_id}}{The unique identifier for the experiment.}
#'   \item{\code{nmi}}{Normalized Mutual Information between clustering and reference labels.}
#'   \item{\code{ari}}{Adjusted Rand Index between clustering and reference labels.}
#'   \item{\code{acc}}{Accuracy computed from the confusion matrix (typically maximum matching).}
#' }
#'
#' @details This function loads all available experiment metadata using \code{load_experiments()},
#' extracts clustering labels for the selected experiments, and evaluates performance using
#' \code{calculate_external_metrics()}. Accuracy is derived from the confusion matrix result.
#'
#' @examples
#' \dontrun{
#' ids <- c("exp001", "exp002", "exp003")
#' reference <- iris$Species
#' calculate_external_metrics_from_ids(ids, reference)
#' }
#'
#' @export
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

#' Internal Validation of Selected Clusterings via diceR Consensus Evaluation
#'
#' Evaluates a set of clustering results (identified by experiment IDs) using internal validation metrics
#' provided by the \code{consensus_evaluate()} function from the diceR package.
#'
#' @param experiment_ids A character vector of experiment identifiers corresponding to stored clustering results.
#' @param data A numeric matrix or data frame of the input dataset on which clustering was performed.
#' @param nk An integer indicating the number of clusters used (\code{k}) in each experiment.
#' @param seed An optional integer for reproducibility of the evaluation process. Default is \code{101}.
#'
#' @return An object returned by \code{consensus_evaluate()}, typically a list containing internal metrics
#' and processed clustering ensemble results.
#'
#' @details This function retrieves clustering labels for the specified experiments using \code{load_experiments()},
#' formats them into a 4D array structure expected by \code{consensus_evaluate()}, and runs internal validation
#' (e.g., connectivity, silhouette, Dunn index). Each clustering solution is treated as a separate algorithm run
#' under a fixed value of \code{k}.
#'
#' @examples
#' \dontrun{
#' ids <- c("exp01", "exp02")
#' result <- evaluate_experiment_clusterings(ids, iris[, -5], nk = 3)
#' result$ii  # View internal index scores
#' }
#'
#' @export
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


#' Select Best Consensus Combination Based on Internal Validation
#'
#' Evaluates all consensus clustering combinations associated with a set of experiment IDs
#' and selects the best-performing ensemble per clustering method, using internal metrics.
#' Then compares the best from each method and returns the overall best performer.
#'
#' @param experiment_ids A character vector of experiment IDs (as returned by \code{load_experiments()})
#' to be evaluated.
#' @param data A numeric data matrix or data frame used as input to clustering, with observations as rows.
#' @param nk Integer. The number of clusters expected in the data.
#' @param seed Integer seed for reproducibility. Default is \code{101}.
#'
#' @return A list returned by \code{consensus_evaluate()} for the final evaluation round, including:
#' \describe{
#'   \item{\code{trim.obj}}{Internal metric evaluation with best combination identified.}
#'   \item{\code{best_experiment_id}}{The experiment ID corresponding to the best consensus combination across all methods.}
#' }
#'
#' @details This function first identifies the best consensus method (based on internal validation)
#' within each clustering algorithm. It then re-evaluates the top-scoring configurations to determine
#' the overall most stable and internally valid combination.
#'
#' Requires that \code{load_experiments()} returns a data frame with columns:
#' \code{experiment_id}, \code{clustering_method}, \code{consensus_method}, and \code{labels_clustering}.
#'
#' @examples
#' \dontrun{
#' ids <- c("exp001", "exp002", "exp003", "exp004")
#' data_matrix <- iris[, -5]
#' result <- get_best_ensemble_internal_metrics(ids, data_matrix, nk = 3)
#' result$best_experiment_id
#' result$trim.obj$ii
#' }
#'
#' @export
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

#' Select Best Consensus Combination Based on Internal Validation
#'
#' Evaluates all consensus clustering combinations associated with a set of experiment IDs
#' and selects the best-performing ensemble per clustering method, using internal metrics.
#' Then compares the best from each method and returns the overall best performer.
#'
#' @param experiment_ids A character vector of experiment IDs (as returned by \code{load_experiments()})
#' to be evaluated.
#' @param data A numeric data matrix or data frame used as input to clustering, with observations as rows.
#' @param nk Integer. The number of clusters expected in the data.
#' @param seed Integer seed for reproducibility. Default is \code{101}.
#'
#' @return A list returned by \code{consensus_evaluate()} for the final evaluation round, including:
#' \describe{
#'   \item{\code{trim.obj}}{Internal metric evaluation with best combination identified.}
#'   \item{\code{best_experiment_id}}{The experiment ID corresponding to the best consensus combination across all methods.}
#' }
#'
#' @details This function first identifies the best consensus method (based on internal validation)
#' within each clustering algorithm. It then re-evaluates the top-scoring configurations to determine
#' the overall most stable and internally valid combination.
#'
#' Requires that \code{load_experiments()} returns a data frame with columns:
#' \code{experiment_id}, \code{clustering_method}, \code{consensus_method}, and \code{labels_clustering}.
#'
#' @examples
#' \dontrun{
#' ids <- c("exp001", "exp002", "exp003", "exp004")
#' data_matrix <- iris[, -5]
#' result <- get_best_ensemble_internal_metrics(ids, data_matrix, nk = 3)
#' result$best_experiment_id
#' result$trim.obj$ii
#' }
#'
#' @export
run_clustering_ensemble_experiments <- function(data_all = data,
                                                dataset = dataset_name,
                                                top_features = top_features,
                                                nk = nk,
                                                algorithms = algorithms,
                                                cons.funs = c("kmodes", "majority", "CSPA", "LCE", "LCA"),
                                                UFS_method = "Inf-FS2020",
                                                alpha = 0.1,
                                                seed = 101) {

  # Feature selection
  if (is.null(top_features)) {
    data <- data_all
  } else {
    data <- data_all[, top_features]
  }

  set.seed(seed)
  experiment_ids <- c()

  # Run DICE ensemble clustering with all algorithms
  execution_time <- system.time({
    results_dice <- dice(data = data,
                         k.method = nk,
                         nk = nk,
                         algorithms = algorithms,
                         cons.funs = cons.funs,
                         evaluate = TRUE,
                         seed = seed,
                         seed.data = seed,
                         progress = FALSE,
                         verbose = FALSE)
  })["elapsed"]

  # Iterate over consensus methods
  for (consensus_method in colnames(results_dice$clusters)) {
    consensus_labels <- as.vector(results_dice$clusters[, consensus_method])

    # Calculate internal clustering metrics
    internal_metrics <- calculate_internal_metrics(data, consensus_labels, nk)

    # Create description and store experiment
    description <- paste0("Ensemble Clustering: ", paste(algorithms, collapse = ", "),
                          " + ", UFS_method, " + ",
                          consensus_method, " + seed ", seed)

    exp_data <- experiment_logger(
      description = description,
      dataset = dataset,
      clustering_method = paste(algorithms, collapse = ", "),
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
    experiment_ids <- c(experiment_ids, exp_data$experiment_id)
  }

  return(experiment_ids)
}

#' Run Ensemble Clustering Experiments and Log Results
#'
#' Executes ensemble clustering on a dataset using multiple base algorithms and consensus functions
#' (via the \code{dice()} function), evaluates internal metrics, and stores the results as logged experiments.
#'
#' @param data_all A numeric matrix or data frame with the full feature set. Defaults to \code{data}.
#' @param dataset A character string identifying the name of the dataset (used in metadata).
#' @param top_features A vector of feature indices or names used for clustering. If \code{NULL}, all features are used.
#' @param nk Integer. The number of clusters for the clustering algorithms.
#' @param algorithms A character vector of base clustering algorithms supported by \code{dice()} (e.g., \code{"kmeans"}, \code{"spectral"}).
#' @param cons.funs A character vector of consensus functions. Defaults to \code{c("kmodes", "majority", "CSPA", "LCE", "LCA")}.
#' @param UFS_method A string indicating the unsupervised feature selection method used. Default is \code{"Inf-FS2020"}.
#' @param alpha A numeric value passed as a hyperparameter to the UFS method. Default is \code{0.1}.
#' @param seed An integer used to set the random seed for reproducibility. Default is \code{101}.
#'
#' @return A character vector of experiment IDs corresponding to each clustering + consensus combination.
#'
#' @details The function prepares data based on selected features, runs ensemble clustering using \code{dice()},
#' computes internal clustering metrics using \code{calculate_internal_metrics()}, and logs the results using
#' \code{experiment_logger()} and \code{save_experiment()}. Each consensus method result is stored as a separate experiment.
#'
#' @examples
#' \dontrun{
#' selected_feats <- c(1, 5, 7, 12)
#' algorithms <- c("kmeans", "spectral")
#' ids <- run_clustering_ensemble_experiments(
#'   data_all = iris[, -5],
#'   dataset = "iris",
#'   top_features = selected_feats,
#'   nk = 3,
#'   algorithms = algorithms
#' )
#' }
#'
#' @export
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

#' Run Base Clustering Experiments with Multiple Algorithms and Consensus Functions
#'
#' Executes ensemble clustering for each specified base algorithm individually, using a range of consensus functions.
#' Each combination is internally evaluated and logged as a separate experiment for further analysis.
#'
#' @param data_all A numeric matrix or data frame containing the original dataset. Defaults to \code{data}.
#' @param dataset A character string representing the dataset's name (used for metadata logging).
#' @param top_features A numeric or character vector of selected feature indices or names. If \code{NULL}, all features are used.
#' @param nk Integer specifying the number of clusters to use across all runs.
#' @param algorithms A character vector of clustering algorithms to evaluate (e.g., \code{"kmeans"}, \code{"spectral"}).
#' @param cons.funs A character vector of consensus functions to apply. Defaults to \code{c("kmodes", "majority", "CSPA", "LCE", "LCA")}.
#' @param UFS_method A character string indicating the unsupervised feature selection method used. Defaults to \code{"Inf-FS2020"}.
#' @param alpha A numeric value passed to the UFS method as a tunable hyperparameter. Default is \code{0.1}.
#' @param seed Integer seed for random number generation. Ensures reproducibility. Default is \code{101}.
#'
#' @return A character vector of experiment IDs representing the stored results from each clustering + consensus combination.
#'
#' @details This function loops through each base clustering algorithm in \code{algorithms}, applies it via \code{dice()},
#' and evaluates each resulting consensus solution using internal metrics. Each result is logged through
#' \code{experiment_logger()} and saved via \code{save_experiment()}. This is useful for benchmarking individual
#' base learners across different consensus strategies.
#'
#' @examples
#' \dontrun{
#' ids <- run_clustering_experiments(
#'   data_all = iris[, -5],
#'   dataset = "iris",
#'   top_features = NULL,
#'   nk = 3,
#'   algorithms = c("kmeans", "pam")
#' )
#' }
#'
#' @export
run_RPClu_experiments <- function(data_all = data,
                                  dataset = dataset_name,
                                  top_features = top_features,
                                  nk = nk,
                                  B = 500,
                                  B.star = 50,
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


#' Execute Full Clustering Pipeline with and without Feature Selection
#'
#' Runs a comprehensive clustering pipeline that compares the performance of multiple ensemble clustering strategies—
#' with and without unsupervised feature selection (UFS)—and includes RPClu as a competing method. Results are logged
#' as experiments with associated metadata and internal validation metrics.
#'
#' @param data A numeric matrix or data frame with observations in rows and features in columns.
#' @param dataset_name A character string used to label the dataset across experiments.
#' @param top_features A vector of selected feature indices or names to use for clustering with UFS.
#' @param nk Integer. Number of clusters to use in all clustering methods.
#' @param algorithms A character vector of clustering algorithms to evaluate (e.g., \code{"kmeans"}, \code{"spectral"}).
#' @param UFS_method A string indicating the unsupervised feature selection technique used (default: \code{"Inf-FS2020"}).
#' @param alpha Numeric hyperparameter passed to the UFS method (default: \code{0.1}).
#' @param n_reps Not currently used. Reserved for future extension with multiple seeds or replicates. Default is \code{1}.
#' @param seed Integer seed for reproducibility (default: \code{100}).
#'
#' @return A named list containing three experiment IDs:
#' \describe{
#'   \item{\code{best_experiment_UFS_id}}{ID of the best clustering experiment using UFS.}
#'   \item{\code{best_experiment_NO_UFS_id}}{ID of the best experiment without UFS.}
#'   \item{\code{best_experiment_UFS_RPClu_id}}{ID of the RPClu experiment using UFS.}
#' }
#'
#' @details This pipeline performs three evaluations:
#' \enumerate{
#'   \item Clustering with unsupervised feature selection and selection of the best ensemble.
#'   \item Clustering without feature selection, using all available features.
#'   \item Execution of RPClu ensemble clustering using selected features.
#' }
#' Internal clustering metrics are computed for all approaches using \code{calculate_internal_metrics()},
#' and the results are saved via \code{save_experiment()}.
#'
#' @examples
#' \dontrun{
#' run_pipeline(
#'   data = iris[, -5],
#'   dataset_name = "iris",
#'   top_features = c(1, 2, 3),
#'   nk = 3,
#'   algorithms = c("kmeans", "spectral")
#' )
#' }
#'
#' @export
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
                                         B = 500,
                                         B.star = 50,
                                         UFS_method = "Inf-FS2020",
                                         alpha = alpha,
                                         seed = 101)
  print(experiment_id)


  return(list(best_experiment_UFS_id = best_experiment_UFS_id,
              best_experiment_NO_UFS_id = best_experiment_NO_UFS_id,
              best_experiment_UFS_RPClu_id = experiment_id
  ))
}


#' Plot Evaluation Metrics vs Number of Features
#'
#' Generates line plots for clustering evaluation metrics (internal or external)
#' across varying numbers of features, facilitating visual comparison of performance.
#' Metrics are labeled to indicate whether higher values are preferred (+) or lower (-).
#'
#' @param data A data frame or matrix (currently not used inside the function).
#' @param internal_metrics Logical. If `TRUE`, uses internal metrics; otherwise uses external metrics.
#' @param title Character string. Title of the plot.
#' @param xAxis_text Character string. Label for the x-axis (usually describing number of features).
#' @param file_save Character string. Filename for saving the plot (EPS format), saved under `experiments/` directory.
#'
#' @return A ggplot object showing metric trends across feature sets.
#' @export
#'
#' @examples
#' plot_metrics_vs_num_features(
#'   data = NULL,
#'   internal_metrics = TRUE,
#'   title = "Internal Metrics vs Number of Features",
#'   xAxis_text = "Number of Selected Features",
#'   file_save = "internal_metrics_plot.eps"
#' )
plot_metrics_vs_num_features <- function (data, internal_metrics, title, xAxis_text, file_save) {
  # Define which metrics are better when higher (+) and which are better when lower (-)
  positive_metrics <- c("calinski_harabasz", "dunn", "pbm", "tau", "gamma", "silhouette", "ensemble_ari")
  negative_metrics <- c("c_index", "davies_bouldin", "mcclain_rao", "g_plus", "sd_dis", "ray_turi", "Compactness", "Connectivity")

  # Convert the list of metrics into columns in a dataframe
  if (internal_metrics) {
    metrics_df <- do.call(rbind, lapply(experiments_infFS$internal_metrics, as.data.frame))
  } else {
    metrics_df <- do.call(rbind, lapply(experiments_infFS$external_metrics, as.data.frame))
  }

  # Add the 'num_features' column
  metrics_df$num_features <- experiments_infFS$num_features

  # Convert dataframe to long format for easier plotting
  long_metrics <- pivot_longer(metrics_df, cols = -num_features, names_to = "metric", values_to = "value")

  # Remove the 's_dbw' metric from the dataframe
  long_metrics <- long_metrics[long_metrics$metric != "s_dbw", ]

  # Assign labels based on the metric type
  long_metrics$metric_label <- ifelse(long_metrics$metric %in% positive_metrics,
                                      paste0(long_metrics$metric, " (+)"),
                                      ifelse(long_metrics$metric %in% negative_metrics,
                                             paste0(long_metrics$metric, " (-)"),
                                             long_metrics$metric))

  # Generate plots with updated labels
  plot <- ggplot(long_metrics, aes(x = num_features, y = value)) +
    geom_line() +
    facet_wrap(~ metric_label, scales = "free_y") +
    theme_minimal() +
    labs(title = title,
         x = xAxis_text,
         y = "Metric value") +
    theme(plot.title = element_text(hjust = 0.5))

  ggsave(paste0("experiments/", file_save), plot = plot,
         width = 8, height = 6, device = "eps")
  return(plot)
}

#' Plot Comparison of Evaluation Metrics vs Number of Features
#'
#' Compares evaluation metrics (internal or external) across two different datasets
#' based on the number of selected features. Generates line plots to show trends,
#' with metric labels indicating whether higher (+) or lower (-) values are preferred.
#'
#' @importFrom dplyr bind_rows
#' @importFrom tidyr pivot_longer
#'
#' @param data1 A list containing metric results and `num_features` for the first dataset.
#' @param legend_data1 Character string. Legend label for the first dataset.
#' @param data2 A list containing metric results and `num_features` for the second dataset.
#' @param legend_data2 Character string. Legend label for the second dataset.
#' @param legend_title Character string. Title for the legend that differentiates the datasets.
#' @param internal_metrics Logical. If `TRUE`, internal metrics are used; otherwise external metrics.
#' @param title Character string. Title of the plot.
#' @param xAxis_text Character string. Label for the x-axis.
#' @param file_save Character string. Filename to save the plot (EPS format) in the `experiments/` folder.
#'
#' @return A ggplot object visualizing selected clustering metrics across features and datasets.
#' @export
#'
#' @examples
#' plot_metrics_vs_num_features2(
#'   data1 = list(internal_metrics = list(...), num_features = c(...)),
#'   legend_data1 = "Method A",
#'   data2 = list(internal_metrics = list(...), num_features = c(...)),
#'   legend_data2 = "Method B",
#'   legend_title = "Method",
#'   internal_metrics = TRUE,
#'   title = "Internal Metrics Comparison",
#'   xAxis_text = "Number of Selected Features",
#'   file_save = "metrics_comparison.eps"
#' )
plot_metrics_vs_num_features2 <- function(data1, legend_data1,
                                          data2, legend_data2,
                                          legend_title,
                                          internal_metrics,
                                          title, xAxis_text, file_save,
                                          baseline = NULL) {
  # Define which metrics are better when higher (+) and which are better when lower (-)
  positive_metrics <- c("calinski_harabasz", "dunn", "pbm", "tau", "gamma", "silhouette", "ensemble_ari")
  negative_metrics <- c("c_index", "davies_bouldin", "mcclain_rao", "g_plus", "sd_dis", "ray_turi", "Compactness", "Connectivity")

  # Function to extract metrics from data and add source label
  process_data <- function(data, source_name) {
    metrics_df <- if (internal_metrics) {
      do.call(rbind, lapply(data$internal_metrics, as.data.frame))
    } else {
      do.call(rbind, lapply(data$external_metrics, as.data.frame))
    }

    metrics_df$num_features <- data$num_features
    metrics_df$source <- source_name  # Label for differentiation
    return(metrics_df)
  }

  # Process both datasets with their corresponding legend labels
  metrics_df1 <- process_data(data1, legend_data1)
  metrics_df2 <- process_data(data2, legend_data2)

  # Combine datasets
  combined_metrics <- bind_rows(metrics_df1, metrics_df2)

  # Convert to long format
  long_metrics <- pivot_longer(combined_metrics, cols = -c(num_features, source), names_to = "metric", values_to = "value")

  # Remove the 's_dbw' metric
  # long_metrics <- long_metrics[long_metrics$metric != "s_dbw", ]
  long_metrics <- long_metrics[long_metrics$metric %in% c("silhouette", "davies_bouldin", "calinski_harabasz", "c_index", "ensemble_ari"), ]


  # Assign metric labels
  long_metrics$metric_label <- ifelse(long_metrics$metric %in% positive_metrics,
                                      paste0(long_metrics$metric, " (+)"),
                                      ifelse(long_metrics$metric %in% negative_metrics,
                                             paste0(long_metrics$metric, " (-)"),
                                             long_metrics$metric))

  # Generate plot with colors for each dataset and custom legend title
  plot <- ggplot(long_metrics, aes(x = num_features, y = value, color = source)) +
    geom_line() +
    facet_wrap(~ metric_label, scales = "free_y") +
    theme_minimal() +
    labs(title = title,
         x = xAxis_text,
         y = "Metric value",
         color = legend_title) +  # Custom legend title
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "bottom" )

  if (!is.null(baseline)) {
    plot <- plot +
      geom_hline(yintercept = baseline, linetype = "dashed", color = "#1E3A5F") +
      annotate("text",
               x = min(combined_metrics$num_features),
               y = baseline,
               label = paste0("Baseline = ", baseline),
               hjust = 0, vjust = -0.5,
               color = "#1E3A5F", size = 4)
  }

  # Save plot
  ggsave(paste0("experiments/", file_save), plot = plot,
         width = 8, height = 6, device = "eps")

  return(plot)
}
