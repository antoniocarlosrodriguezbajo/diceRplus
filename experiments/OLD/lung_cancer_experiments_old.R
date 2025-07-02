library(diceRplus)
library(dplyr)
data(lung_cancer)

# Verify normalizacion
mean(lung_cancer$x)
sd(lung_cancer$x)


# Calcute internal_metrics based on consensus_evaluate of diceR
calculate_internal_metrics <- function(data,cluster_labels) {
  num_labels <- length(cluster_labels)

  # Create structure for consensus_evaluate
  cc_data <- array(cluster_labels, dim = c(num_labels, 1, 1, 1))
  # row_names <- rownames(data)
  row_names <- if (!is.null(rownames(data))) rownames(data) else seq_len(nrow(data))
    dimnames(cc_data) <- list(
    row_names,  # Primer nivel de nombres: nombres de las filas de Meat$x
    "R1",       # Repetition
    "RPGMMClu", # Clustering algorithm
    "5"         # Number of clusters
  )
  # Evaluation
  result_evaluation <- consensus_evaluate(data, cc_data)
  metrics_df <- result_evaluation$ii[[1]]
  metrics_list <- lapply(metrics_df[1, -1], function(x) unname(x))
  return(metrics_list)
}

calculate_redundancy_jaccard <- function(X, ranking, porcentages = c(0.05, 0.10, 0.15, 0.20, 0.30),
                                         k_Jaccard = 5, alpha = 0.5) {
  # Ensure ranking is integer and within valid range
  ranking <- as.integer(ranking)
  n_features <- ncol(X)
  n_samples <- nrow(X)

  # Initialize list to store results
  results <- list()

  # Loop over each specified percentage
  for (pct in porcentages) {
    # Calculate the number of features to select
    n_feat <- max(1, round(n_features * pct))

    # Select the top features according to the ranking
    selected <- ranking[1:n_feat]
    X_sel <- X[, selected, drop = FALSE]

    ## --- Compute Redundancy ---
    # Compute the correlation matrix between features
    cor_mat <- cor(X_sel)
    # Take upper triangle (excluding diagonal) to avoid redundancy and self-correlation
    cor_vals <- cor_mat[upper.tri(cor_mat)]
    # Redundancy is the average absolute correlation
    redundancy <- mean(abs(cor_vals))

    ## --- Compute Jaccard Index ---
    # Compute similarity matrix on selected features (Euclidean by default)
    dist_sel <- as.matrix(dist(X_sel))
    sel_neighbors <- apply(dist_sel, 1, function(row) order(row)[2:(k_Jaccard+1)])

    dist_all <- as.matrix(dist(X))
    all_neighbors <- apply(dist_all, 1, function(row) order(row)[2:(k_Jaccard+1)])

    jaccard_vec <- numeric(n_samples)
    for (i in 1:n_samples) {
      inter <- length(intersect(sel_neighbors[,i], all_neighbors[,i]))
      union <- length(union(sel_neighbors[,i], all_neighbors[,i]))
      jaccard_vec[i] <- ifelse(union == 0, 0, inter / union)
    }
    jaccard_score <- mean(jaccard_vec)

    ## --- Compute Explained Variance ---
    # Perform PCA on selected features
    pca_result <- prcomp(X_sel, center = TRUE, scale. = TRUE)
    # Explained variance ratio
    explained_variance <- sum(pca_result$sdev^2) / sum(prcomp(X, center = TRUE, scale. = TRUE)$sdev^2)

    ## --- Compute Entropy ---
    # Normalize the data for entropy calculation
    X_norm <- apply(X_sel, 2, function(x) (x - min(x)) / (max(x) - min(x) + 1e-10)) # Avoid division by zero
    entropy_values <- apply(X_norm, 2, function(x) -sum(x * log(x + 1e-10))) # Shannon entropy
    entropy <- mean(entropy_values) # Average entropy across features

    ## --- Composite Index ---
    composite_index <- alpha * (1 - redundancy) + (1 - alpha) * jaccard_score

    # Store results
    results[[length(results) + 1]] <- data.frame(
      porcentage = pct,
      n_features = n_feat,
      redundancy = redundancy,
      jaccard = jaccard_score,
      explained_variance = explained_variance,
      entropy = entropy,
      composite_index = composite_index
    )
  }

  # Combine all results into a single data.frame
  result_df <- do.call(rbind, results)
  rownames(result_df) <- NULL
  return(result_df)
}


# results_all: data.frame con columnas method, porcentage, redundancy, jaccard

get_best_ufs_percent <- function(results_all) {
  # 1. Para cada método, encuentra el porcentaje con el mayor Jaccard
  library(dplyr)
  best_per_method <- results_all %>%
    group_by(method) %>%
    # En caso de empate en jaccard, selecciona el de menor redundancia
    arrange(desc(jaccard), redundancy) %>%
    slice(1) %>%
    ungroup()

  # 2. Entre los mejores de cada método, elige el de mayor Jaccard (si empate, menor redundancia)
  best_overall <- best_per_method %>%
    arrange(desc(jaccard), redundancy) %>%
    slice(1)

  # 3. Devuelve el resultado
  return(list(
    best_method = best_overall$method,
    best_percentage = best_overall$porcentage,
    best_jaccard = best_overall$jaccard,
    best_redundancy = best_overall$redundancy
  ))
}

# Ejemplo de uso:
# results_all <- data.frame(method = rep(...), porcentage = ..., redundancy = ..., jaccard = ...)
# resultado <- get_best_ufs_percent(results_all)
# print(resultado)


# Calcute scores based on internal metrics
calculate_internal_scores <- function(metrics_list, weights=NULL) {
  # Convert metrics_list to a dataframe
  metrics_df <- as.data.frame(do.call(rbind, metrics_list))
  # Convert column to numeric
  metrics_df <- as.data.frame(lapply(metrics_df, function(x) as.numeric(unlist(x))))
  # Restore methods' name
  rownames(metrics_df) <- names(metrics_list)
  # Normalize metrics' values
  normalize <- function(x) (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  metrics_norm <- as.data.frame(lapply(metrics_df, normalize))
  # Restore names
  rownames(metrics_norm) <- rownames(metrics_df)
  # If NA set then set metric to 0
  metrics_norm[is.na(metrics_norm)] <- 0

  # Invert metrics where lower values indicate better performance
  metrics_norm$davies_bouldin <- 1 - metrics_norm$davies_bouldin
  metrics_norm$Compactness <- 1 - metrics_norm$Compactness
  metrics_norm$Connectivity <- 1 - metrics_norm$Connectivity
  metrics_norm$c_index <- 1 - metrics_norm$c_index
  metrics_norm$mcclain_rao <- 1 - metrics_norm$mcclain_rao
  metrics_norm$sd_dis <- 1 - metrics_norm$sd_dis
  metrics_norm$ray_turi <- 1 - metrics_norm$ray_turi

  # Use default weights if none are provided
  if (is.null(weights)) {
    weights <- c(
      calinski_harabasz = 0.15,
      dunn = 0.2,
      pbm = 0.05,
      tau = 0,
      gamma = 0,
      c_index = 0,
      mcclain_rao = 0,
      sd_dis = 0,
      ray_turi = 0,
      g_plus = 0,
      silhouette = 0.35,
      s_dbw = 0,
      davies_bouldin = 0.1,
      Compactness = 0.1,
      Connectivity = 0.05
    )
  }
  # Calculate scores with weights
  # Match weights to columns
  valid_weights <- weights[names(weights) %in% colnames(metrics_norm)]
  metrics_use <- metrics_norm[, names(valid_weights), drop = FALSE]
  scores <- rowSums(sweep(metrics_use, 2, valid_weights, "*"))
  return(scores)
}

select_UFS_methods <- function(data,
                               threshold_variance  = 0.1,           # Default threshold for low variance
                               threshold_variance_PCA = 0.9,        # Default threshold for explained variance in PCA
                               threshold_high_dimensionality = 1000,# Default threshold for high-dimensional data
                               threshold_correlation = 0.5,         # Default threshold for high correlation
                               threshold_low_variance = 0.3,        # Default threshold for proportion of low-variance features
                               threshold_pca_components = 0.2       # Default threshold for number of PCA components
) {
  # Correlation
  cor_mat <- cor(data, use = "pairwise.complete.obs")
  mean_cor <- mean(abs(cor_mat[upper.tri(cor_mat)]))

  # Variance
  var_vec <- apply(data, 2, var, na.rm = TRUE)
  low_var_prop <- mean(var_vec < quantile(var_vec, threshold_variance))

  # Identify the minimum number of components (PCA) needed
  # to capture at least a threshold of the variance
  pca <- prcomp(data, scale. = TRUE)
  var_exp <- summary(pca)$importance[2, ]
  cum_var_threshold <- as.numeric(which(cumsum(var_exp) > threshold_variance_PCA)[1])

  # Select UFS candidates
  ufs_candidates <- c()

  if(mean_cor > threshold_correlation) ufs_candidates <- c(ufs_candidates, "MCFS", "InfFS", "Laplacian")
  if(low_var_prop > threshold_low_variance) ufs_candidates <- c(ufs_candidates, "UDFS", "NDFS", "InfFS")
  if(cum_var_threshold < ncol(data) * threshold_pca_components) ufs_candidates <- c(ufs_candidates, "UDFS", "NDFS", "InfFS")
  if(ncol(data) > threshold_high_dimensionality) ufs_candidates <- c(ufs_candidates, "Laplacian", "SPEC", "InfFS")

  # Cover at least 3 types of UFS
  ufs_candidates <- unique(c(ufs_candidates, "UDFS", "InfFS", "Laplacian", "NDFS"))

  return(list(ufs_candidates = ufs_candidates, cum_var_threshold = cum_var_threshold))
}

run_UFS_methods_experiments <- function(data, UFS_results, N, g, B, B.star, weights=NULL) {
  scores <- list()
  metrics_list <- list()

  for (method in names(UFS_results$Results)) {
    print(method)
    top_features <- UFS_results$Results[[method]]$Result[[1]][1:N]
    execution_time <- system.time(out.clu_result <- RPGMMClu_parallel(data$x[, top_features],
                                                                        data$y,
                                                                        g=g,
                                                                        B=B,
                                                                        B.star=B.star,
                                                                        verb=TRUE))["elapsed"]
    # Calculate internal metrics
    internal_metrics <- calculate_internal_metrics(data$x[, top_features],
                                                   out.clu_result$ensemble$label.vec)

    cat("ARI:", out.clu_result$ensemble$ari, "\n")
    #cat("BIC:", out.clu_result$ensemble$bic.final, "\n")
    #print(internal_metrics)
    metrics_list[[method]] <- internal_metrics
  }

  scores <- calculate_internal_scores(metrics_list, weights)
  return(scores)

}

run_optimized_experiments <- function(data, UFS_features, N_start, N_max, step, g, B, B.star, metrics_to_track, patience=30) {
  best_metrics <- list()
  best_N <- N_start
  no_improvement_count <- 0
  metrics_list <- list()

  for (N in seq(N_start, N_max, by = step)) {
    print(paste("Evaluating N =", N))

    top_features <- UFS_features[1:N]
    execution_time <- system.time(out.clu_result <- RPGMMClu_parallel(data$x[, top_features],
                                                                      data$y,
                                                                      g = g,
                                                                      B = B,
                                                                      B.star = B.star,
                                                                      verb = TRUE))["elapsed"]

    # Print ARI
    ari_value <- out.clu_result$ensemble$ari
    print(paste("ARI for N =", N, ":", ari_value))

    internal_metrics <- calculate_internal_metrics(data$x[, top_features], out.clu_result$ensemble$label.vec)
    metrics_list[[paste0("N_", N)]] <- internal_metrics

    improvement <- FALSE
    for (metric in metrics_to_track) {
      if (!is.null(best_metrics[[metric]]) && internal_metrics[[metric]] > best_metrics[[metric]]) {
        best_metrics[[metric]] <- internal_metrics[[metric]]
        improvement <- TRUE
      } else if (is.null(best_metrics[[metric]])) {
        best_metrics[[metric]] <- internal_metrics[[metric]]
      }
    }

    if (improvement) {
      best_N <- N
      no_improvement_count <- 0
    } else {
      no_improvement_count <- no_improvement_count + 1
    }

    if (no_improvement_count >= patience) {
      print(paste("Stopping at N =", best_N, "due to no improvement"))
      break
    }
  }

  return(list(N = N, metrics_list=metrics_list))
}


# Parameter Num clusters
# UDFS: Yes
# NDFS: No
# InfFS: No
# Laplacian: No
# SPEC: No
# MCFS: No

# Run just once
#UFS_results <- runUFS_manual_params(lung_cancer$x, ufs_candidates)
#save(UFS_results, file = "experiments/UFS_lung_cancer.RData")

# Load UFS_results
load("experiments/UFS_lung_cancer.RData")

result <- select_UFS_methods(lung_cancer$x)
ufs_candidates <- result$ufs_candidates
cum_var_threshold <- result$cum_var_threshold

for (method in ufs_candidates) {
  print(method)
  UFS_features <- UFS_results$Results[[method]]$Result[[1]][1:ncol(lung_cancer$x)]
  results <- calculate_redundancy_jaccard(lung_cancer$x, UFS_features, k_Jaccard=20)
  print(results)
}

dice_results <- list() # To store results per method

for (method in ufs_candidates) {
  print(method)
  candidates <- ufs_candidates[[i]]

  # Select top features using names or indices, depending on your UFS_results structure
  # If Result[[1]] is a vector of names, this will subset by name; if indices, it still works in R
  UFS_features <- UFS_results$Results[[method]]$Result[[1]][1:cum_var_threshold]

  # Subset the candidate data for selected features
  selected_data <- lung_cancer$x[, UFS_features]

  # Run dice clustering with GMM and kmodes consensus
  dice_obj <- dice(
    data = selected_data,
    nk = 5,            # Range of clusters to try (can be changed)
    reps = 100,           # Number of bootstrap repetitions (adjust as needed)
    algorithms = "gmm",  # Use Gaussian Mixture Model for clustering [[5]], [[12]]
    cons.funs = "kmodes",# Consensus method; can also use "majority", "CSPA", etc. [[12]], [[13]]
    progress = TRUE,      # Show progress bar
    trim = TRUE,
    reweigh = TRUE,
    ref.cl =  as.integer(factor(lung_cancer$y))
  )

  # Store the result under the method's name
  dice_results[[method]] <- dice_obj
}

lapply(dice_results, function(x) x$indices$ii)
adjustedRandIndex(lung_cancer$y, dice_results$InfFS$clusters[,"kmodes"])

N <- cum_var_threshold
g <- 5
B <- 250
B.star <- round(B/10)


cat("N:", N, "\ng:", g, "\nB:", B, "\nB.star:", B.star, "\n")

method = "Laplacian"
top_features <- UFS_results$Results[[method]]$Result[[1]][1:N]
execution_time <- system.time(out.clu_baseline <- RPGMMClu_noens_parallel(lung_cancer$x[, top_features],
                                                                    lung_cancer$y,
                                                                    g=5,
                                                                    B=B,
                                                                    B.star=B.star,
                                                                    verb=TRUE))["elapsed"]

E <- as.matrix(out.clu_baseline$label.vec)
k = 5
# Etiquetas de consenso por varios métodos
consensus_labels_kmodes <- k_modes(E, is.relabelled = TRUE, seed = 1)
consensus_labels_majority <- majority_voting(E, is.relabelled = TRUE)
consensus_labels_cspa <- CSPA(E, k)
consensus_labels_lce <- LCE(E, k, sim.mat = "cts")
consensus_labels_lca <- LCA(E, is.relabelled = TRUE, seed = 1)


adjustedRandIndex(consensus_labels_kmodes, lung_cancer$y)
adjustedRandIndex(consensus_labels_majority, lung_cancer$y)
adjustedRandIndex(consensus_labels_lce, lung_cancer$y)
adjustedRandIndex(consensus_labels_lca, lung_cancer$y)

calculate_internal_metrics(lung_cancer$x, consensus_labels_kmodes)
calculate_internal_metrics(lung_cancer$x, consensus_labels_majority)
calculate_internal_metrics(lung_cancer$x, consensus_labels_lce)

y <- consensus_matrix(E)
PAC(y, lower = 0.05, upper = 0.95)

# Consensus clustering for multiple algorithms
set.seed(911)
x <- matrix(rnorm(500), ncol = 10)
CC <- consensus_cluster(x, nk = 3:4, reps = 10, algorithms = c("ap", "km"),
                        progress = FALSE)

# Obtener la matriz de clases de consenso usando, por ejemplo, majority_voting
# Para un valor específico de k (por ejemplo, k=3):
cons_cl <- majority_voting(CC[, , 1, 1, drop = FALSE], is.relabelled = FALSE)
cons_cl_mat <- matrix(cons_cl, ncol = 1)
# Crear un vector de referencia
set.seed(1)
ref.cl <- sample(1:4, 50, replace = TRUE)

z <- consensus_evaluate(
  data = x,
  cons.cl = matrix(cons_cl, ncol = 1),  # matriz de clases de consenso
  ref.cl = ref.cl,        # referencia, si la tienes
  n = 1,
  trim = TRUE
)
str(z, max.level = 2)


# Simulación de datos
set.seed(1)
x <- matrix(rnorm(100*10), nrow = 100)
cons_cl <- sample(1:4, 100, replace = TRUE)
ref.cl <- sample(1:4, 100, replace = TRUE)

# Evaluación SIN trimming (lo recomendado para una sola columna)
z <- consensus_evaluate(
  data = x,
  cons.cl = matrix(cons_cl, ncol = 1),
  ref.cl = ref.cl
)

cons.cl <- cbind(
  kmodes = consensus_labels_kmodes,
  majority = consensus_labels_majority
)

cons.cl <- as.matrix(cons.cl)

storage.mode(cons.cl) <- "integer"

result_eval <- consensus_evaluate(
  data = lung_cancer$x,
  cons.cl = cons.cl,
  ref.cl = as.integer(lung_cancer$y),
  k.method = 5
)

# Duplicar artificialmente el vector de consenso
cons_cl_mat <- cbind(cons_cl, cons_cl)
z <- consensus_evaluate(
  data = x,
  cons.cl = cons_cl_mat,
  ref.cl = ref.cl
)

# Ver resumen de resultados:
str(result_eval, max.level = 2)

print(dim(lung_cancer$x))
print(dim(cons.cl))
print(length(lung_cancer$y))
print(str(cons.cl))
print(head(cons.cl))
print(str(lung_cancer$y))

scores <- run_UFS_methods_experiments(lung_cancer, UFS_results, N, g, B, B.star)
print(scores)
# Method with best score
best_method <- names(which.max(scores))
cat("The best method based on internal metrics score is:", best_method, "\n")

UFS_features <- UFS_results$Results[[best_method]]$Result[[1]][1:ncol(lung_cancer$x)]

N_start <- round(cum_var_threshold/2, -1)
N_max   <- round(cum_var_threshold*2, -1)
step <- 10
g <- 5
B <- 1000
B.star <- round(B/10)
metrics_to_track = list("calinski_harabasz", "dunn", "silhouette")

cat("N_start:", N_start,
    "\nN_max:", N_max,
    "\nstep:", step,
    "\ng:", g,
    "\nB:", B,
    "\nB.star:", B.star,
    "\nmetrics_to_track:", paste(metrics_to_track, collapse=", "), "\n")


results <- run_optimized_experiments(lung_cancer, UFS_features, N_start, N_max, step, g, B, B.star, metrics_to_track)

weights <- c(
  calinski_harabasz = 0.15,
  dunn = 0.2,
  pbm = 0.05,
  tau = 0,
  gamma = 0,
  c_index = 0,
  mcclain_rao = 0,
  sd_dis = 0,
  ray_turi = 0,
  g_plus = 0,
  silhouette = 0.35,
  s_dbw = 0,
  davies_bouldin = 0.1,
  Compactness = 0.1,
  Connectivity = 0.05
)


internal_scores <- calculate_internal_scores(results$metrics_list)
internal_scores

################################
################################
################################

scores <- list()
metrics_list <- list()

for (method in names(UFS_results$Results)) {
  print(method)
  top_features <- UFS_results$Results[[method]]$Result[[1]][1:N]
  execution_time <- system.time(out.clu_baseline <- RPGMMClu_parallel(lung_cancer$x[, top_features],
                                                                      lung_cancer$y,
                                                                      g=5,
                                                                      B=B,
                                                                      B.star=B.star,
                                                                      verb=TRUE))["elapsed"]
  # Calculate internal metrics
  internal_metrics <- calculate_internal_metrics(lung_cancer$x[, top_features],
                                                 out.clu_baseline_UFS$ensemble$label.vec)

  cat("ARI:", out.clu_baseline$ensemble$ari, "\n")
  cat("BIC:", out.clu_baseline$ensemble$bic.final, "\n")
  print(internal_metrics)
  metrics_list[[method]] <- internal_metrics
}

# Convertir cada columna a numérica
metrics_df <- as.data.frame(lapply(metrics_df, function(x) as.numeric(unlist(x))))

# Restaurar nombres de los métodos
rownames(metrics_df) <- names(metrics_list)

# Ahora sí, aplicar la normalización
normalize <- function(x) (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
metrics_norm <- as.data.frame(lapply(metrics_df, normalize))

# Restaurar nombres de los métodos en la tabla normalizada
rownames(metrics_norm) <- rownames(metrics_df)

metrics_norm[is.na(metrics_norm)] <- 0

# Invert metrics where lower values indicate better performance
metrics_norm$davies_bouldin <- 1 - metrics_norm$davies_bouldin
metrics_norm$Compactness <- 1 - metrics_norm$Compactness
metrics_norm$Connectivity <- 1 - metrics_norm$Connectivity
metrics_norm$c_index <- 1 - metrics_norm$c_index
metrics_norm$mcclain_rao <- 1 - metrics_norm$mcclain_rao
metrics_norm$sd_dis <- 1 - metrics_norm$sd_dis
metrics_norm$ray_turi <- 1 - metrics_norm$ray_turi



# Default weights
weights <- c(
  calinski_harabasz = 0.2,
  dunn = 0.15,
  pbm = 0.2,
  tau = 0,
  gamma = 0,
  c_index = 0,
  mcclain_rao = 0,
  sd_dis = 0,
  ray_turi = 0,
  g_plus = 0,
  silhouette = 0.15,
  s_dbw = 0,
  davies_bouldin = 0.1,
  Compactness = 0.1,
  Connectivity = 0.1
)

scores <- rowSums(sweep(metrics_norm, 2, weights, "*"))  # Multiplica por pesos

# Method with best score
best_method <- names(which.max(scores))
cat("The best method based on internal metrics is:", best_method, "\n")


# Fixed UFS: Lapacian
method <- names(UFS_results$Results)[4]
print(method)

# B <- 500
# B.star <- 50


# Iterate over different values of N (from 20 to 100 in steps of 10)
for (N in seq(20, cum_var_90*4, by = 10)) {
  print(N)
  top_features <- UFS_results$Results[[method]]$Result[[1]][1:N]

  # Run RPGMMClu_parallel
  execution_time <- system.time(out.clu_baseline <- RPGMMClu_parallel(lung_cancer$x[, top_features],
                                                                          lung_cancer$y,
                                                                          g = 5,
                                                                          B = B,
                                                                          B.star = B.star))["elapsed"]
  # Calculate internal metrics
  internal_metrics <- calculate_internal_metrics(lung_cancer$x_filtered,
                                                     out.clu_baseline_UFS$ensemble$label.vec)

  # Log the experiment
  exp_data <- experiment_logger(
    description = paste("Clustering with RPGMMClu - UFS:", method, "- Features:", N),
    dataset = "Meat",
    ensemble_method = "RPGMMClu_parallel",
    ensemble_method_params = list(g = 5, B = B, B.star = B.star),
    UFS_method = method,
    UFS_method_params = list(default = 'YES'),
    num_features = N,
    features = top_features,
    execution_time = as.numeric(execution_time),
    labels_clustering = out.clu_baseline_UFS$ensemble$label.vec,
    #internal_metrics = internal_metrics_UFS,
    external_metrics = list(ensemble_ari = out.clu_baseline_UFS$ensemble$ari[[1]])
  )
  cat("N:", N, "\nARI:", out.clu_baseline$ensemble$ari, "\n")
  print(internal_metrics)

  # save_experiment(exp_data)
}



# Fixed UFS: InfFS
method <- names(UFS_results$Results)[1]
print(method)
N=1000
top_features <- UFS_results$Results[[method]]$Result[[1]][1:N]

B <- 500
B.star=50

# Run RPGMMClu
execution_time <- system.time(out.clu_baseline <- RPGMMClu_parallel(lung_cancer$x[, top_features],
                                                                    lung_cancer$y,
                                                                    g=5,
                                                                    B=B,
                                                                    B.star=B.star,
                                                                    verb=TRUE))["elapsed"]


# 1. Armar los datos
Ns <- c(60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220)
internal_metrics <- data.frame(
  N = Ns,
  calinski_harabasz = c(57.87697, 51.59859, 46.37122, 42.01917, 34.62241, 37.68687, 42.05708, 41.16215, 41.55947, 39.58954, 32.70785, 38.30188, 38.75257, 37.43411, 37.34429, 37.45155, 36.25126),
  dunn = c(0.2334573, 0.2441958, 0.232152, 0.2319141, 0.2276068, 0.2533581, 0.3098038, 0.3144582, 0.3287489, 0.3098608, 0.2653664, 0.3245465, 0.3379175, 0.3328103, 0.3300723, 0.3292157, 0.3331432),
  pbm = c(100.0493, 91.91092, 93.25978, 90.91565, 94.79159, 95.76134, 90.43483, 93.59372, 94.7651, 86.96051, 91.75608, 88.97748, 98.00736, 91.2206, 90.4391, 89.49346, 87.29137),
  davies_bouldin = c(1.611013,1.742295,1.822133,2.059175,1.88417,1.749407,1.721634,1.733397,1.660232,1.701205,1.798366,1.764388,1.735302,1.776442,1.796006,1.672096,1.923498),
  silhouette = c(0.2087521,0.1935313,0.196098,0.1958769,0.2084995,0.2343272,0.2134696,0.2179294,0.2224218,0.2028203,0.2191014,0.2013259,0.2135716,0.2101525,0.197515,0.2023932,0.178178),
  Compactness = c(15.27308,16.42041,17.30482,18.10112,19.27864,19.40635,19.39378,19.72146,19.98496,20.49254,21.48226,21.15103,21.29808,21.6439,21.85902,22.05195,22.3309),
  Connectivity = c(115.9583,116.5611,119.6143,103.0925,101.419,76.41627,80.40079,83.98492,65.13611,91.43968,80.50119,86.67619,78.69206,77.41587,88.41984,75.84246,96.20794)
)
ARI <- c(
  0.282262972975235, 0.176268096030117, 0.297607963558047, 0.325196772173436, 0.533146358353904, 0.467886784908866,
  0.385844113947494, 0.381068527414226, 0.400375652692011, 0.380177834183966, 0.595739317125691, 0.386131132419878,
  0.393785510381912, 0.387607522713476, 0.387307502016007, 0.442162127047557, 0.356817526400704
)
internal_metrics\$ARI <- ARI

# 2. Normalizar los índices internos (opcional, pero recomendable)
internal_metrics_norm <- as.data.frame(scale(internal_metrics[,2:8]))
internal_metrics_norm\$ARI <- internal_metrics\$ARI
internal_metrics_norm\$N <- internal_metrics\$N

# 3. Definir nombres de índices internos
nombres_indices <- c("calinski_harabasz", "dunn", "pbm",
                     "davies_bouldin", "silhouette", "Compactness", "Connectivity")

# 4. Función de optimización
optimizar_pesos <- function(metrics, nombres_indices, ARI) {
  fn_objetivo <- function(pesos) {
    # Normaliza pesos
    pesos <- abs(pesos) / sum(abs(pesos))
    scores <- as.numeric(as.matrix(metrics[, nombres_indices]) %*% pesos)
    -cor(scores, ARI, use = "complete.obs")
  }
  pesos_iniciales <- rep(1/length(nombres_indices), length(nombres_indices))
  res <- optim(pesos_iniciales, fn_objetivo)
  pesos_opt <- abs(res\$par) / sum(abs(res\$par))
  names(pesos_opt) <- nombres_indices
  return(pesos_opt)
}

# 5. Ejecutar la optimización
pesos_ajustados <- optimizar_pesos(internal_metrics_norm, nombres_indices, internal_metrics_norm\$ARI)

print(pesos_ajustados)

[1] "Evaluating N = 60"
[1] "ARI for N = 60 : 0.290050154366465"
[1] "Evaluating N = 70"
[1] "ARI for N = 70 : 0.187788574383325"
[1] "Evaluating N = 80"
[1] "ARI for N = 80 : 0.304401387120705"
[1] "Evaluating N = 90"
[1] "ARI for N = 90 : 0.387829065298249"
[1] "Evaluating N = 100"
[1] "ARI for N = 100 : 0.37757872028926"
[1] "Evaluating N = 110"
[1] "ARI for N = 110 : 0.425436171792138"
[1] "Evaluating N = 120"
[1] "ARI for N = 120 : 0.383183749730404"
[1] "Evaluating N = 130"
[1] "ARI for N = 130 : 0.415012805006948"
[1] "Evaluating N = 140"
[1] "ARI for N = 140 : 0.385836396928524"
[1] "Evaluating N = 150"
[1] "ARI for N = 150 : 0.381006308540034"
[1] "Evaluating N = 160"
[1] "ARI for N = 160 : 0.388813277239087"
[1] "Evaluating N = 170"
[1] "ARI for N = 170 : 0.394334009870363"
[1] "Evaluating N = 180"
[1] "ARI for N = 180 : 0.395252549678362"
[1] "Evaluating N = 190"
[1] "ARI for N = 190 : 0.378635676718746"
[1] "Evaluating N = 200"
[1] "ARI for N = 200 : 0.378123862916848"
[1] "Evaluating N = 210"
[1] "ARI for N = 210 : 0.394133466312472"
[1] "Evaluating N = 220"
[1] "ARI for N = 220 : 0.403822782791122"

> B
[1] 1000

> internal_scores
N_60      N_70      N_80      N_90     N_100     N_110     N_120     N_130     N_140     N_150     N_160     N_170     N_180     N_190
0.6806323 0.3398183 0.1564344 0.5631280 0.5252703 0.6416572 0.3180427 0.7094260 0.5597850 0.4893351 0.3955159 0.3512323 0.5776050 0.4966191
N_200     N_210     N_220
0.4856763 0.4479814 0.5761606

