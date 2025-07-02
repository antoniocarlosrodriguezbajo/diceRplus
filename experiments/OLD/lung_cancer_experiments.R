library(diceRplus)
data(lung_cancer)
data(ALLAML)
data(Meat)
data(prostate_ge)
library(mclust)

data <- prostate_ge
k=2
UFS_results_file <- "experiments/UFS_prostate_ge.RData"


# Verify normalizacion

mean(data$x)
sd(data$x)
data$x<- scale(data$x)


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

  if(mean_cor > threshold_correlation) ufs_candidates <- c(ufs_candidates, "MCFS", "Inf-FS2020", "Laplacian")
  if(low_var_prop > threshold_low_variance) ufs_candidates <- c(ufs_candidates, "UDFS", "NDFS", "Inf-FS2020")
  if(cum_var_threshold < ncol(data) * threshold_pca_components) ufs_candidates <- c(ufs_candidates, "UDFS", "NDFS", "Inf-FS2020")
  if(ncol(data) > threshold_high_dimensionality) ufs_candidates <- c(ufs_candidates, "Laplacian", "SPEC", "Inf-FS2020")

  # Cover at least 3 types of UFS
  ufs_candidates <- unique(c(ufs_candidates, "UDFS", "InfFS", "Laplacian", "NDFS"))

  return(list(ufs_candidates = ufs_candidates, cum_var_threshold = cum_var_threshold))
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




result <- select_UFS_methods(data$x)
ufs_candidates <- result$ufs_candidates
cum_var_threshold <- result$cum_var_threshold

# Run just once
UFS_results <- runUFS_manual_params(data$x, ufs_candidates)
save(UFS_results, file = UFS_results_file)


# Load UFS_results
load(UFS_results_file)
names(UFS_results$Results)

N <- cum_var_threshold
g <- k
B <- 300
B.star <- round(B/10)

if (N < 100) N <- 100
N<- 1000

cat("N:", N, "\ng:", g, "\nB:", B, "\nB.star:", B.star, "\n")


# List of consensus functions
consensus_methods <- list(
  kmodes   = function(E, k)      k_modes(E, is.relabelled = FALSE, seed = 1),
  majority = function(E, k)      majority_voting(E, is.relabelled = FALSE),
  cspa     = function(E, k)      CSPA(E, k),
  lce      = function(E, k)      LCE(E, k, sim.mat = "cts"),
  lca      = function(E, k)      LCA(E, is.relabelled = FALSE, seed = 1)
)




cc_list <- list()
methods <- names(UFS_results$Results)
for (method in methods) {
  tryCatch({

    print(method)

    top_features <- UFS_results$Results[[method]]$Result[[1]][1:N]

    execution_time <- system.time(out.clusters <- RPGMMClu_noens_parallel(data$x[, top_features],
                                                                          data$y,
                                                                          g=g,
                                                                          B=B,
                                                                          B.star=B.star,
                                                                          verb=TRUE))["elapsed"]

    E <- as.matrix(out.clusters$label.vec)
    E_4d <- array(E, dim = c(dim(E)[1], dim(E)[2], 1, 1))
    rownames <- paste0("Sample", seq_len(dim(E)[1]))
    repsnames <- paste0("R", seq_len(dim(E)[2]))
    algonames <- "Algorithm"
    knames <- as.character(k)
    dimnames(E_4d) <- list(rownames, repsnames, algonames, knames)

    E_list <- list(kmodes=E, majority=E, cspa=E_4d, lce=E, lca=E)

    consensus_labels <- mapply(function(f, EE) f(EE, k), consensus_methods, E_list, SIMPLIFY = FALSE)

    cat("Method:", method, " B:", B, "\n")
    print("ARI")
    ari_scores <- sapply(consensus_labels, adjustedRandIndex, data$y)
    print(ari_scores)

    cc_this <- mapply(
      function(lbls, name) prepare_consensus_evaluation(data$x, lbls, paste(method, name, B, sep = "+")),
      consensus_labels, names(consensus_labels), SIMPLIFY = FALSE
    )

    cc_list <- c(cc_list, cc_this)

  }, error = function(e) {
    cat("Error con método:", method, "\n")
    print(e$message)
  })
}

set.seed(123)
result_evaluation <- do.call(consensus_evaluate, c(
  list(data$x),
  cc_list,
  list(ref.cl = as.integer(data$y), trim = TRUE, n = 1, reweigh=TRUE)
))

print(result_evaluation$trim.obj$top.list)
print(result_evaluation$trim.obj$rank.matrix)
print(result_evaluation$ii)

