library(diceRplus)
data(prostate_ge)

data <- prostate_ge

mean(data$x)
sd(data$x)

UFS_results_file <- "experiments/UFS_prostate_ge.RData"

UFS_Methods$MCFS <- NULL
UFS_Methods$SOGFS <- NULL
UFS_Methods$UDFS <- NULL
UFS_Methods$FMIUFS <- NULL
UFS_Methods$UAR_HKCMI <- NULL
UFS_Methods$U2FS <- NULL
UFS_Methods$NDFS <- NULL

# DGUFS default parameters
# Number of clusters: 2
# Alpha: 0.5
# Beta: 0.9
# Dimension of selected features: 2

# Run just once
ufs_candidates = list("Inf-FS2020")
UFS_results_file <- "experiments/UFS_prostate_ge_Inf-FS2020.RData"
UFS_results <- runUFS_manual_params(data$x, ufs_candidates)
save(UFS_results, file = UFS_results_file)

# Run just once
UFS_results <- runUFS(data$x, UFS_Methods)
save(UFS_results, file = UFS_results_file)

load(UFS_results_file)

cc_list <- list()
k_modes_classes <- list()
N <- 500
k = 2

for (method in names(UFS_results$Results)){
  print(method)
  top_features <- UFS_results$Results[[method]]$Result[[3]]

  cc_list[[method]] <- consensus_cluster(data$x[, top_features], nk = k,
                                         p.item = 1,reps = 3, algorithms = "km",
                                         progress = FALSE)

}

ref_method <- names(UFS_results$Results)[1]
ref_matrix <- cc_list[[ref_method]][,,1,1]  # Matriz: filas = muestras, columnas = repeticiones

for (method in names(cc_list)) {
  # Extrae la matriz de etiquetas para todas las repeticiones
  mat <- cc_list[[method]][,,1,1]  # dimensiones: muestras x repeticiones
  # Reetiqueta cada repetición usando la correspondiente de la referencia
  for (rep in seq_len(ncol(mat))) {
    mat[, rep] <- relabel_class(mat[, rep], ref_matrix[, 1])
  }
  # Sobrescribe en la lista
  cc_list[[method]][,,1,1] <- mat
}




load(UFS_results_file)
UFS_results$Results$CFS <- NULL
UFS_results$Results$RNE <- NULL
UFS_results$Results$CNAFS <- NULL
UFS_results$Results$EGCFS <- NULL
cc_list <- list()
N <- 50
k = 2

for (method in names(UFS_results$Results)){
  print(method)
  top_features <- UFS_results$Results[[method]]$Result[[1]][1:N]

  cc_list[[method]] <- consensus_cluster(data$x[, top_features], nk = k,
                                         p.item = 1,reps = 1, algorithms = "km",
                                         progress = FALSE)
}

# Suponiendo que cc_list ya está creado y contiene los resultados de consensus_cluster
ref_method <- names(UFS_results$Results)[1]
ref_labels <- cc_list[[ref_method]][,1,1,1]  # Extrae las etiquetas del primer clustering del método de referencia

for (method in names(cc_list)) {
  # Extrae las etiquetas del clustering para el método actual
  labels <- cc_list[[method]][,1,1,1]
  # Relabel para que coincidan lo mejor posible con el método de referencia
  relabeled <- relabel_class(labels, ref_labels)
  # Guarda el resultado relabelado
  cc_list[[method]][,1,1,1] <- relabeled
}


ufs_names <- names(cc_list)

for (i in seq_along(cc_list)) {
  dn <- dimnames(cc_list[[i]])
  # Ajust name alg + UFS method
  dn[[3]] <- paste0("KM-", ufs_names[i])
  dimnames(cc_list[[i]]) <- dn
}

results_combined <- do.call(
  consensus_evaluate,
  c(
    list(
      data = data$x,
      ref.cl = as.integer(data$y),
      n = 5,
      trim = TRUE
    ),
    cc_list
  )
)

# External metrics
results_combined$ei

top_features <- UFS_results$Results[[method]]$Result[[3]]
# Run dice
dice.obj <- dice(
  data =  data$x[, top_features],
  nk = 2,
  p.item = 1,
  reps = 100,
  algorithms = c("km", "pam"),
  ref.cl = as.integer(data$y),
  evaluate = TRUE,
  plot = FALSE,
  progress = TRUE,
  verbose = TRUE
)

top_features <- UFS_results$Results[[method]]$Result[[3]]
# Run dice
dice.obj <- dice(
  data = data$x[, top_features],
  nk = 2,
  algorithms = c("pam"),
  ref.cl = as.integer(data$y),
  evaluate = TRUE,
  plot = FALSE,
  progress = TRUE,
  verbose = TRUE
)


dice.obj$indices$ei
#############################
#############################


# Combine the results of UFS methods
res_combined <- do.call(consensus_evaluate, c(
  list(data = data$x, ref.cl = as.integer(data$y)), cc_list)
)
str(res_combined)

res_combined$ei

# Evaluate and compare the results
eval <- do.call(consensus_evaluate, c(list(data$x), cc_list))
str(eval, max.level = 2)















# Run just once
UFS_Results <- runUFS(data$x, UFS_Methods)
save(UFS_Results, file = "experiments/UFS_Meat.RData")


x_prep <- prepare_data(data$x, scale = TRUE, type = "robust", min.var = 0.01)
k <- 2

# List of algorithms to use
algos <- c("nmf", "hc", "diana", "km", "pam", "ap", "sc", "gmm", "block", "som", "cmeans", "hdbscan")
algos <- c("nmf", "hc", "diana", "km", "pam", "ap", "sc", "gmm", "block", "hdbscan")

# Run dice
dice.obj <- dice(
  data =x_prep,
  nk = 2,
  reps = 10,
  algorithms = algos,
  ref.cl = as.integer(data$y),
  evaluate = TRUE,
  trim = TRUE,
  reweigh = TRUE,
  n=5,
  plot = FALSE,
  progress = TRUE,
  verbose = TRUE
)

dice.obj$indices$trim$top.list

dice.obj$indices$ei[["2"]]


E.new <- dice.obj$indices$trim$trim.obj$E.new

dice.obj$clusters[,"kmodes"]

ii2 <- dice.obj$indices$ii[["2"]]
ii2["kmodes", ]
ei2 <- dice.obj$indices$ei[["2"]]
ei2[, ]


CC <- consensus_cluster(data$x, nk = k, p.item = 0.8, reps = 5,
                        algorithms = c("hc", "pam", "diana"))
co <- capture.output(str(CC))
strwrap(co, width = 80)

CC[, , "HC_Euclidean",]
dice.obj$indices$ii

dice.obj$indices$ei


CC_imputed <- apply(CC, 2:4, impute_knn, data = data$x, seed = 1)

CC_imputed[, , "HC_Euclidean",]


pam <- CC[, , "PAM_Euclidean", "2", drop = FALSE]
cm <- consensus_matrix(pam)

ccomb_matrix <- consensus_combine(CC, element = "matrix")
ccomb_class <- consensus_combine(CC, element = "class")

ccomb_matrix[["2"]][["HC_Euclidean"]]
ccomb_class[["2"]]


super_cc <- consensus_matrix(ccomb_class)



ccomp <- consensus_evaluate(data$x, CC, plot = FALSE)

ccomp
