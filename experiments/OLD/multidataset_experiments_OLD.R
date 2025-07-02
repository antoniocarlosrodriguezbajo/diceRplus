library(diceRplus)
library(mclust)
library(cluster)

UFS_Methods <- list(
  "Inf-FS2020" = TRUE
)

algorithms = c("nmf", "hc", "diana", "km", "pam", "ap", "sc",
               "gmm", "block", "som", "cmeans", "hdbscan")



library(MASS)  # Required for generalized inverse

library(MASS)  # Required for generalized inverse

InfFS_U <- function(X_train, alpha=0.1, verbose = FALSE) {
  if (verbose) {
    cat("\n+ Feature selection method: Inf-FS Unsupervised (2020)\n")
  }

  # 1) Compute Spearman correlation
  if (verbose) cat("1) Priors/weights estimation \n")

  corr_ij <- cor(X_train, method = "spearman", use = "pairwise.complete.obs")
  corr_ij[is.na(corr_ij)] <- 0
  corr_ij[is.infinite(corr_ij)] <- 0
  corr_ij <- 1 - abs(corr_ij)

  # Standard deviation estimation
  STD <- apply(X_train, 2, sd)
  STDMatrix <- pmax(outer(STD, STD, FUN = "pmax"), STD)
  STDMatrix <- STDMatrix - min(STDMatrix)
  sigma_ij <- STDMatrix / max(STDMatrix, na.rm = TRUE)
  sigma_ij[is.na(sigma_ij)] <- 0
  sigma_ij[is.infinite(sigma_ij)] <- 0

  # 2) Build the graph G = <V,E>
  if (verbose) cat("2) Building the graph G = <V,E> \n")

  N <- ncol(X_train)
  eps <- 5e-06 * N
  factor <- 1 - eps

  A <- alpha * sigma_ij + (1 - alpha) * corr_ij
  rho <- max(rowSums(A))

  # Substochastic Rescaling
  A <- A / (max(rowSums(A)) + eps)
  stopifnot(max(rowSums(A)) < 1)

  # Inf-FS Core: Infinite Paths
  I <- diag(N)
  r <- factor / rho
  y <- I - (r * A)
  S <- ginv(y)  # Using generalized inverse to match MATLAB's inv()

  # 5) Compute energy scores
  WEIGHT <- rowSums(S)

  # 6) Rank features according to importance
  RANKED <- order(WEIGHT, decreasing = TRUE)

  e <- rep(1, N)
  t <- S %*% e

  # Histogram binning and threshold setting
  nbins <- round(0.5 * N)
  counts <- hist(t, breaks = nbins, plot = FALSE)$counts
  thr <- mean(counts)  # Matching MATLAB threshold logic
  size_sub <- sum(counts > thr)

  cat(sprintf("Inf-FS (U) Nb. Features Selected = %.4f (%.2f%%)\n", size_sub, 100 * (size_sub / N)))

  SUBSET <- RANKED[1:size_sub]

  return(list(RANKED = RANKED, WEIGHT = WEIGHT, SUBSET = SUBSET))
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

generate_random_projections <- function(x, d = NULL, c = 10, B = 1000, seed = 101) {
  p <- ncol(x)
  if (is.null(d)) {
    d <- ceiling(c * log(p))  # Ajuste para usar p
  }
  print(d)
  set.seed(seed)
  RPbase <- generateRP(p, d, B)  # Genera B matrices de proyección aleatorias

  # Aplicación de las proyecciones en un solo paso
  projections <- array(x %*% RPbase, dim = c(nrow(x), d, B))

  return(projections)
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




############################
# Meat
############################
data(Meat)

data <- Meat

mean(data$x)
sd(data$x)

results_UFS <- InfFS_U(data$x, alpha=0.2)

# Run just once
UFS_results_file <- "experiments/UFS_Meat_Inf-FS2020.RData"
# UFS_results <- runUFS(data$x, UFS_Methods)
# save(UFS_results, file = UFS_results_file)

load(UFS_results_file)

method = "Inf-FS2020"
top_features <- UFS_results$Results[[method]]$Result[[3]]

top_features <- results_UFS$SUBSET

sprintf("The %s method selected %d features.", method, length(top_features))


g <- 5
B <- 100
B.star <- 100
execution_time <- system.time(out.clu_baseline <- RPGMMClu_noens_parallel(data$x[, top_features],
                                                                          NULL,
                                                                          g=5,
                                                                          B=B,
                                                                          B.star=B.star,
                                                                          verb=TRUE))["elapsed"]

# Internal Metrics
E <- as.matrix(out.clu_baseline$label.vec)

E_4d <- array(E, dim = c(dim(E)[1], dim(E)[2], 1, 1))
rownames <- paste0("Sample", seq_len(dim(E)[1]))
repsnames <- paste0("RP1", seq_len(dim(E)[2]))
algonames <- "Algorithm"
knames <- as.character(k)
dimnames(E_4d) <- list(rownames, repsnames, algonames, knames)

E_list <- list(kmodes=E, majority=E, cspa=E_4d, lce=E, lca=E)

k = 5
# List of consensus functions
consensus_methods <- list(
  kmodes   = function(E, k)      k_modes(E, is.relabelled = FALSE, seed = 1),
  majority = function(E, k)      majority_voting(E, is.relabelled = FALSE),
  cspa     = function(E, k)      CSPA(E, k),
  lce      = function(E, k)      LCE(E, k, sim.mat = "cts"),
  lca      = function(E, k)      LCA(E, is.relabelled = FALSE, seed = 1)
)

consensus_labels <- mapply(function(f, EE) f(EE, k), consensus_methods, E_list, SIMPLIFY = FALSE)

cc_list <- mapply(
  function(lbls, name) prepare_consensus_evaluation(data$x, lbls, paste(method, name, B, sep = "+")),
  consensus_labels, names(consensus_labels), SIMPLIFY = FALSE
)

set.seed(123)
result_evaluation <- do.call(consensus_evaluate, c(
  list(data$x),
  cc_list,
  list(ref.cl = as.integer(data$y), trim = TRUE, n = 1, reweigh=TRUE)
))

print(result_evaluation$trim.obj$top.list)
print(result_evaluation$trim.obj$rank.matrix)
print(result_evaluation$ii)

print(result_evaluation$ei)

adjustedRandIndex(cc_list$lca, as.integer(data$y))

# Direct with consensus labels

# Etiquetas de consenso por varios métodos
consensus_labels_kmodes <- k_modes(E, is.relabelled = FALSE, seed = 1)
consensus_labels_majority <- majority_voting(E, is.relabelled = FALSE)
consensus_labels_lce <- LCE(E, k, sim.mat = "cts")
consensus_labels_lca <- LCA(E, is.relabelled = FALSE, seed = 1)

adjustedRandIndex(consensus_labels_kmodes, as.integer(data$y))
adjustedRandIndex(consensus_labels_majority, as.integer(data$y))
adjustedRandIndex(consensus_labels_lce, as.integer(data$y))
adjustedRandIndex(consensus_labels_lca, as.integer(data$y))


# NMI
ev_nmi(consensus_labels_kmodes, as.integer(data$y))
ev_nmi(consensus_labels_majority, as.integer(data$y))
ev_nmi(consensus_labels_lce, as.integer(data$y))
ev_nmi(consensus_labels_lca, as.integer(data$y))

# Otras métricas externas (accuracy, kappa, etc.)
ev_confmat(consensus_labels_kmodes, as.integer(data$y))
ev_confmat(consensus_labels_majority, as.integer(data$y))
ev_confmat(consensus_labels_lce,as.integer(data$y))
ev_confmat(consensus_labels_lca, as.integer(data$y))



g <- 5
B <- 100
execution_time <- system.time(out.clu_baseline <- RPClu_parallel(data$x[, top_features],
                                                                 clust_fun = "gmm",
                                                                 g=5,
                                                                 B=B,
                                                                 verb=TRUE))["elapsed"]

E <- as.matrix(out.clu_baseline$clusterings)

k=5
# Etiquetas de consenso por varios métodos
consensus_labels_kmodes <- k_modes(E, is.relabelled = FALSE, seed = 1)
consensus_labels_majority <- majority_voting(E, is.relabelled = FALSE)
consensus_labels_lce <- LCE(E, k, sim.mat = "cts")
consensus_labels_lca <- LCA(E, is.relabelled = FALSE, seed = 1)

adjustedRandIndex(consensus_labels_kmodes, as.integer(data$y))
adjustedRandIndex(consensus_labels_majority, as.integer(data$y))
adjustedRandIndex(consensus_labels_lce, as.integer(data$y))
adjustedRandIndex(consensus_labels_lca, as.integer(data$y))


E <- as.matrix(out.clu_baseline$clusterings)

k <- 5
ref_labels <- as.integer(data$y)

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


eval_results

############################
# ALLAML
############################
data(ALLAML)

data <- ALLAML

mean(data$x)
sd(data$x)


top_features <- cgufs_selection(data$x, k1=2,k2=2,quantile = 0.9)
sprintf("The %s method selected %d features.", "cgufs", length(top_features))


# Run dice
dice.obj <- dice(
  data =  data$x[,top_features],
  reps=100,
  p.item = 1,
  nk = 2,
  algorithms = c("km"),
  ref.cl = as.integer(data$y),
  evaluate = TRUE,
  progress = TRUE,
  verbose = TRUE,

)

dice.obj$indices$ei

# Run just once
UFS_results_file <- "experiments/UFS_ALLAML_Inf-FS2020.RData"
# UFS_results <- runUFS(data$x, UFS_Methods)
# save(UFS_results, file = UFS_results_file)

# Run just once
# DGUFS default parameters
# Number of clusters: 2
# Alpha: 0.5
# Beta: 0.9
# Dimension of selected features: 150
#ufs_candidates = list("DGUFS")
#UFS_results_file <- "experiments/UFS_ALLAML_ge_DGUFS.RData"
#UFS_results <- runUFS_manual_params(data$x, ufs_candidates)
#save(UFS_results, file = UFS_results_file)

load(UFS_results_file)

top_features <- UFS_results$Results[[method]]$Result[[3]]
sprintf("The %s method selected %d features.", method, length(top_features))

top_features <- cgufs_selection(data$x[,top_features], k1=2,k2=2,quantile = 0.9)
sprintf("The %s method selected %d features.", "cgufs", length(top_features))

# Generar proyecciones aleatorias
rp <- generate_random_projections(data$x[, top_features], B=2000)

k <- 2  # o el número de clusters de tu interés
labels_list <- lapply(1:2000, function(i) {
  km <- kmeans(rp[,,i], centers = k)
  km$cluster
})

compactness_scores <- sapply(1:2000, function(i) {
  compactness(rp[,,i], labels_list[[i]])
})

best_idx <- order(compactness_scores)[1:200]

cl_matrix <- do.call(cbind, labels_list[best_idx])

# k-modes
final_clusters <- k_modes(cl_matrix)
# o majority voting
final_clusters2 <- majority_voting(cl_matrix)

adjustedRandIndex(final_clusters2, as.integer(data$y))




# Calcular coeficientes de silueta con k-means
siluette_scores <- sapply(1:dim(rp)[3], function(i) {
  dat <- rp[,,i]
  km_model <- kmeans(dat, centers=2)  # Ajusta k-means con 2 clusters
  cl <- km_model$cluster  # Extrae etiquetas de cluster

  sil <- silhouette(cl, dist(dat))  # Calcula coeficiente de silueta
  mean(sil[, 3])  # Retorna el promedio de los valores de silueta
})

# Seleccionar las 50 mejores proyecciones
best_idx <- order(siluette_scores, decreasing=TRUE)[1:50]  # Mejores por silueta
best_projections <- lapply(best_idx, function(i) rp[,,i])
rp_combined <- do.call(cbind, best_projections)

# Run dice
dice.obj <- dice(
  data =  rp_combined,
  nk = 2,
  p.item = 1,
  reps = 1,
  algorithms = "km",
  ref.cl = as.integer(data$y),
  evaluate = TRUE,
  plot = FALSE,
  trim = TRUE,
  reweigh = TRUE,
  n = 2,
  progress = TRUE,
  verbose = TRUE
)

dice.obj$indices$ei

# NMF_Brunet  9.887782e-02 0.7222222


g <- 2
B <- 100
B.star <- round(B/10)
execution_time <- system.time(out.clu_baseline <- RPGMMClu_noens_parallel(data$x[, top_features],
                                                                          data$y,
                                                                          g=2,
                                                                          B=B,
                                                                          B.star=B.star,
                                                                          verb=TRUE))["elapsed"]

E <- as.matrix(out.clu_baseline$label.vec)
k = 2
# Etiquetas de consenso por varios métodos
consensus_labels_kmodes <- k_modes(E, is.relabelled = FALSE, seed = 1)
consensus_labels_majority <- majority_voting(E, is.relabelled = FALSE)
# consensus_labels_cspa <- CSPA(E, k)
consensus_labels_lce <- LCE(E, k, sim.mat = "cts")
consensus_labels_lca <- LCA(E, is.relabelled = FALSE, seed = 1)


adjustedRandIndex(consensus_labels_kmodes, as.integer(data$y))
adjustedRandIndex(consensus_labels_majority, as.integer(data$y))
adjustedRandIndex(consensus_labels_lce, as.integer(data$y))
adjustedRandIndex(consensus_labels_lca, as.integer(data$y))


# NMI
ev_nmi(consensus_labels_kmodes, as.integer(data$y))
ev_nmi(consensus_labels_majority, as.integer(data$y))
ev_nmi(consensus_labels_lce, as.integer(data$y))
ev_nmi(consensus_labels_lca, as.integer(data$y))

# Otras métricas externas (accuracy, kappa, etc.)
ev_confmat(consensus_labels_kmodes, as.integer(data$y))
ev_confmat(consensus_labels_majority, as.integer(data$y))
ev_confmat(consensus_labels_lce,as.integer(data$y))
ev_confmat(consensus_labels_lca, as.integer(data$y))



g <- 2
B <- 100
B.star <- round(B/10)
B <- 200
B.star <- 100
execution_time <- system.time(out.clu_baseline <- RPClu_parallel(data$x[, top_features],
                                                                 clust_fun = "km",
                                                                 g=2,
                                                                 B=B,
                                                                 verb=TRUE))["elapsed"]

E <- as.matrix(out.clu_baseline$clusterings)

k <- 2
ref_labels <- as.integer(data$y)

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
  nmi     = ev_nmi(labels, ref_labels)
))

eval_results

############################
# leukemia
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

# Run dice
dice.obj <- dice(
  data =  data$x[, top_features],
  nk = 5,
  p.item = 1,
  reps = 1,
  algorithms = algorithms,
  cons.funs = c("majority"),
  ref.cl = as.integer(data$y),
  evaluate = TRUE,
  plot = FALSE,
  progress = TRUE,
  verbose = TRUE
)

dice.obj$indices$ei
# GMM  0.47031289 0.7339901
# BLOCK 0.02107671 0.665024
# SOM 0.36508799 0.6551724

############################
# lymphoma
############################
data(lymphoma)

data <- lymphoma

mean(data$x)
sd(data$x)

top_features <- cgufs_selection(data$x, k1=3,k2=3,quantile = 0.8)
sprintf("The %s method selected %d features.", "cgufs", length(top_features))


# Run dice
dice.obj <- dice(
  data =  data$x[,top_features],
  nk = 3,
  algorithms = c("km","cmeans"),
  ref.cl = as.integer(data$y+1),
  evaluate = TRUE,
  progress = TRUE,
  verbose = TRUE,

)

dice.obj$indices$ei



# Run just once
UFS_results_file <- "experiments/UFS_lymphoma_Inf-FS2020.RData"
# UFS_results <- runUFS(data$x, UFS_Methods)
# save(UFS_results, file = UFS_results_file)

load(UFS_results_file)

method = "Inf-FS2020"
top_features <- UFS_results$Results[[method]]$Result[[3]]
sprintf("The %s method selected %d features.", method, length(top_features))

# Run dice
dice.obj <- dice(
  data =  data$x[, top_features],
  nk = 3,
  p.item = 1,
  reps = 1,
  algorithms = algorithms,
  cons.funs = c("majority"),
  ref.cl = as.integer(data$y+1),
  evaluate = TRUE,
  plot = FALSE,
  progress = TRUE,
  verbose = TRUE
)

dice.obj$indices$ei

# SC 0.47753735 0.6935484
# SOM   0.40328275 0.6935484


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


# Run dice
dice.obj <- dice(
  data =  data$x[, top_features],
  nk = 2,
  p.item = 1,
  reps = 1,
  algorithms = algorithms,
  cons.funs = c("majority"),
  ref.cl = as.integer(data$y),
  evaluate = TRUE,
  plot = FALSE,
  progress = TRUE,
  verbose = TRUE
)

# External Metrics
dice.obj$indices$ei

# PAM_Euclidean 0.34553179 0.7843137
# KM            0.21329961 0.7450980
# AP            0.34553179 0.7843137
# CMEANS        0.18426776 0.7450980


############################
# COIL20
############################
data(COIL20)

data <- COIL20

mean(data$x)
sd(data$x)

# Run just once
UFS_results_file <- "experiments/UFS_COIL20_Inf-FS2020.RData"
# UFS_results <- runUFS(data$x, UFS_Methods)
# save(UFS_results, file = UFS_results_file)

load(UFS_results_file)

method = "Inf-FS2020"
top_features <- UFS_results$Results[[method]]$Result[[3]]
sprintf("The %s method selected %d features.", method, length(top_features))

# Run dice
dice.obj <- dice(
  data =  data$x[, top_features],
  nk = 20,
  reps = 1,
  algorithms = algorithms,
  cons.funs = c("majority"),
  ref.cl = as.integer(data$y),
  evaluate = TRUE,
  plot = FALSE,
  progress = TRUE,
  verbose = TRUE
)

dice.obj$indices$ei
# MEDIOCRES


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


############################
# warpPIE10P
############################

data(warpPIE10P)

data <- warpPIE10P

data$x <- scale(data$x)

mean(data$x)
sd(data$x)

# Run just once
UFS_results_file <- "experiments/UFS_warpPIE10P_Inf-FS2020.RData"
UFS_results <- runUFS(data$x, UFS_Methods)
save(UFS_results, file = UFS_results_file)

load(UFS_results_file)
