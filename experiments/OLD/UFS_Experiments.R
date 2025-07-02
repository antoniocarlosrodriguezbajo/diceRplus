library(diceRplus)

UFS_methods_fast <- list(
  "InfFS" = FALSE,
  "Laplacian" = TRUE,
  "LLCFS" = TRUE,
  "CFS" = FALSE,
  "FSASL" = TRUE,
  "DGUFS" = TRUE,
  "UFSOL" = TRUE,
  "SPEC" = TRUE,
  "UDFS" = TRUE,
  "SRCFS" = TRUE,
  "RNE" = TRUE,
  "RUFS" = TRUE,
  "NDFS" = TRUE,
  "EGCFS" = TRUE,
  "CNAFS" = TRUE,
  "Inf-FS2020" = TRUE
)


# 2. Auxiliary Functions
redundancy <- function(X) {
  cor_mat <- abs(cor(X))
  cor_mat[lower.tri(cor_mat, diag = TRUE)] <- NA
  mean(cor_mat, na.rm = TRUE)
}

jaccard_similarity <- function(Xorig, Xsel, k = 5) {
  # Compute k-NN neighborhoods for Jaccard between original and selected features
  orig_dist <- as.matrix(dist(Xorig))
  red_dist <- as.matrix(dist(Xsel))
  n <- nrow(orig_dist)
  orig_nb <- lapply(1:n, function(i) order(orig_dist[i, ])[2:(k+1)])
  red_nb  <- lapply(1:n, function(i) order(red_dist[i, ])[2:(k+1)])
  mean(sapply(1:n, function(i) {
    length(intersect(orig_nb[[i]], red_nb[[i]])) / length(union(orig_nb[[i]], red_nb[[i]]))
  }))
}

# 3. Exploratory Analysis
analyze_dataset <- function(X) {
  # a. Redundancy
  mean_cor <- redundancy(X)
  # b. Variance
  vars <- apply(X, 2, var)
  low_var_thresh <- quantile(vars, 0.05)
  low_var_pct <- mean(vars < low_var_thresh)
  # c. PCA
  pca <- prcomp(X, center = TRUE, scale. = TRUE)
  cum_var <- cumsum(pca$sdev^2) / sum(pca$sdev^2)
  num_comp90 <- which(cum_var > 0.9)[1]
  dim_eff <- num_comp90 / ncol(X)
  # d. Number of variables
  nvar <- ncol(X)
  list(mean_cor = mean_cor, low_var_pct = low_var_pct, dim_eff = dim_eff, nvar = nvar)
}

# 4. UFS Method Selection Logic
select_ufs_methods <- function(analysis) {
  methods <- character()
  # High redundancy
  if (analysis$mean_cor > 0.5) methods <- c(methods, "RUFS", "UDFS", "Inf-FS2020", "NDFS")
  # High percentage of low-variance features
  if (analysis$low_var_pct > 0.3) methods <- c(methods, "RUFS", "UDFS", "SOCFS")
  # Low effective dimensionality: redundancy
  if (analysis$dim_eff < 0.2) methods <- c(methods, "SPEC", "Inf-FS2020")
  # High effective dimensionality: need many features
  if (analysis$dim_eff > 0.5) methods <- c(methods, "NDFS")
  # Number of variables
  if (analysis$nvar > 1000) methods <- c(methods, "Laplacian", "SPEC", "Inf-FS2020")
  else if (analysis$nvar > 100) methods <- c(methods, "UDFS", "NDFS", "FSASL")
  else methods <- c(methods, "FSASL", "UFSOL", "FMIUFS", "FRUAR")
  # Ensure diversity
  methods_diversity <- c("UDFS", "Inf-FS2020", "Laplacian", "NDFS")
  methods <- unique(c(methods,methods_diversity))
  methods <- intersect(methods, names(UFS_methods_fast))
  methods
}

# 5. Main Pipeline: Selection and Comparison
# Automatically choose subset sizes based on the number of features, as per the article [[4]]
calculate_subset_sizes <- function(X, base = 2, min_subset_size = 50) {
  n_features = ncol(X)

  if (n_features < 100) {
    subset_sizes <- c(5, 10, 20, 30, 40)
  } else if (n_features < 1000) {
    subset_sizes <- c(50, 100, 150, 200, 300)
  } else {
    max_power <- floor(log(n_features, base))
    subset_sizes <- unique(base^(0:max_power))
    # Añade el total de características si no está incluido
    if (!(n_features %in% subset_sizes)) {
      subset_sizes <- c(subset_sizes, n_features-1)
    }
    # Filtra los tamaños que no exceden n_features
    subset_sizes <- subset_sizes[subset_sizes <= n_features]
    subset_sizes <- subset_sizes[subset_sizes >= min_subset_size]
  }
  return(subset_sizes)
}


ufs_pipeline <- function(X, UFS_Results, subset_sizes) {
  analysis <- analyze_dataset(X)
  cat("Exploratory analysis:\n"); print(analysis)
  methods <- select_ufs_methods(analysis)
  cat("Selected UFS methods:", paste(methods, collapse = ", "), "\n")
  results <- list()
  for (method in methods) {
    for (n_sel in subset_sizes[subset_sizes < ncol(X)]) {
      # Get feature indices from your UFS results object
      idxs <- UFS_Results$Results[[method]]$Result[[1]][1:n_sel]
      Xsel <- X[, idxs, drop = FALSE]
      red <- redundancy(Xsel)
      jac <- jaccard_similarity(X, Xsel)
      results[[paste(method, n_sel, sep = "_")]] <- data.frame(
        method = method, n_selected = n_sel, redundancy = red, jaccard = jac
      )
    }
  }
  do.call(rbind, results)
}

# 6. Example Usage
# Suppose X is your input data matrix (rows=samples, cols=features)
# Suppose UFS_Results is your results object as described

# X <- ... # your data
# UFS_Results <- ... # your UFS results object

# res <- ufs_pipeline(X, UFS_Results, subset_sizes = c(5, 10, 20, 30, 40))
# print(res)

# 7. Final Selection: Methods with lowest redundancy and highest Jaccard per subset size
# library(dplyr)
# best <- res %>%
#   group_by(n_selected) %>%
#   arrange(redundancy, desc(jaccard)) %>%
#   slice(1)
# print(best)


data(ALLAML)
# Check normalization
min(colMeans(ALLAML$x))
max(colMeans(ALLAML$x))
min(apply(ALLAML$x, 2, sd))
max(apply(ALLAML$x, 2, sd))

analysis <- analyze_dataset(ALLAML$x)
cat("Exploratory analysis:\n"); print(analysis)
UFS_methods <- select_ufs_methods(analysis)
cat("Sugested UFS methods:\n"); print(UFS_methods)

# Run just once
#UFS_Results <- runUFS(ALLAML$x, UFS_Methods_fast[UFS_methods])
#save(UFS_Results, file = "experiments/UFS_ALLAML.RData")

#UFS_Results <- runUFS(ALLAML$x, UFS_methods_fast)
#save(UFS_Results, file = "experiments/UFS_ALLAML.RData")


load("experiments/UFS_ALLAML.RData")

subset_sizes <- calculate_subset_sizes(ALLAML$x)

res <- ufs_pipeline(ALLAML$x,UFS_Results,subset_sizes)
print(res)

res_sorted <- res[order(res$n_selected, res$redundancy, -res$jaccard), ]
best <- res_sorted[!duplicated(res_sorted$n_selected), ]

# Fixed UFS: InfFS
method <- "NDFS"
print(method)
N <- 200
B <- 200
B.star=20

for (method in names(UFS_Results$Results)){
# for (method in c("DGUFS")){
  print(method)
  if (method == "CFS") { # RNE??
    top_features <- UFS_Results$Results[["CFS"]]$Result[1:N]
  } else {
    top_features <- UFS_Results$Results[[method]]$Result[[1]][1:N]
  }
  ALLAML$x_filtered <- ALLAML$x[, top_features]
  # Run RPGMMClu
  execution_time <- system.time(out.clu_baseline <- RPGMMClu_parallel(ALLAML$x[, top_features],
                                                                    ALLAML$y,
                                                                    g=2,
                                                                    B=B,
                                                                    B.star=B.star,
                                                                    verb=TRUE))["elapsed"]

  print(out.clu_baseline$ensemble)
}

set.seed(123)
top_features <- UFS_Results$Results[["Laplacian"]]$Result[[1]][1:N]
kmeans_result <- kmeans(ALLAML$x[, top_features], centers = 2)
kmeans_result$cluster
library(mclust)
adjustedRandIndex(ALLAML$y, kmeans_result$cluster)

accuracy <- mean(ALLAML$y == kmeans_result$cluster)
accuracy
library(clue)

# Crear la matriz de confusión entre clusters y etiquetas reales
conf_mat <- table(kmeans_result$cluster, ALLAML$y)

# Resolver la asignación óptima de etiquetas
optimal_assignment <- solve_LSAP(conf_mat, maximum = TRUE)

# Relabeling: Asignar etiquetas correctas al clustering
pred_labels_reassigned <- optimal_assignment[kmeans_result$cluster]
accuracy <- mean(ALLAML$y == pred_labels_reassigned)
accuracy

#####

# Crear la matriz de confusión entre clusters y etiquetas reales
conf_mat <- table(out.clu_baseline$ensemble$label.vec, ALLAML$y)

# Resolver la asignación óptima de etiquetas
optimal_assignment <- solve_LSAP(conf_mat, maximum = TRUE)

# Relabeling: Asignar etiquetas correctas al clustering
pred_labels_reassigned <- optimal_assignment[out.clu_baseline$ensemble$label.vec]
accuracy <- mean(ALLAML$y == pred_labels_reassigned)
accuracy
