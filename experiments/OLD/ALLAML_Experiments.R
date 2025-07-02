library(diceRplus)
data("ALLAML")

UFS_methods <- list(
  "LLCFS"
)


#UFS_Results <- runUFS_manual_params(ALLAML$x, UFS_methods)
#save(UFS_Results, file = "experiments/UFS_ALLAM_LLCFS.RData")

# load("experiments/UFS_ALLAML.RData")
load("experiments/UFS_ALLAM_LLCFS.RData")

method = "LLCFS"
ALLAML$x <- ALLAML$x + 9

# Bucle para iterar sobre valores de N desde 50 hasta 500 en incrementos de 50
results_list <- list()

for (N in seq(25, 300, by = 25)) {

  # Seleccionar las características superiores según el valor de N
  print(N)
  top_features <- UFS_Results$Results[[method]]$Result[[1]][1:N]

  # Ejecutar DICE con los top features seleccionados
  dice.obj <- dice(
    data = ALLAML$x[, top_features],
    nk = 2,
    algorithms = c("nmf"),
    nmf.method = c("brunet"),
    cons.funs = c("CSPA"),
    progress = TRUE,
    verbose = FALSE,
    distance = "manhattan",
    ref.cl = as.integer(as.factor(ALLAML$y))
  )

  # Almacenar los resultados en una lista para referencia posterior
  results_list[[as.character(N)]] <- list(
    ei = dice.obj$indices$ei,
    ii = dice.obj$indices$ii
  )
}


library(ggplot2)
library(tidyr)
library(dplyr)

plot_metrics_vs_num_features <- function(results_list, title, xAxis_text, file_save) {
  # Define which metrics are better when higher (+) and which are better when lower (-)
  positive_metrics <- c("calinski_harabasz", "dunn", "pbm", "tau", "gamma", "silhouette", "ensemble_ari")
  negative_metrics <- c("c_index", "davies_bouldin", "mcclain_rao", "g_plus", "sd_dis", "ray_turi", "Compactness", "Connectivity")

  # Extract internal metrics from results_list with proper character-based indexing
  metrics_df <- do.call(rbind, lapply(names(results_list), function(N) {
    df <- as.data.frame(results_list[[as.character(N)]][["ii"]][["2"]][1,])  # Ensure N is treated as a character
    df$num_features <- as.numeric(N)  # Convert N to numeric for plotting
    return(df)
  }))

  metrics_df <- metrics_df %>%
    select(-Algorithms)

  # Convert dataframe to long format for easier plotting
  long_metrics <- pivot_longer(metrics_df, cols = -num_features, names_to = "metric", values_to = "value")

  # Remove the 's_dbw' metric if present
  long_metrics <- long_metrics[long_metrics$metric != "s_dbw", ]

  # Assign labels based on the metric type
  long_metrics$metric_label <- ifelse(long_metrics$metric %in% positive_metrics,
                                      paste0(long_metrics$metric, " (+)"),
                                      ifelse(long_metrics$metric %in% negative_metrics,
                                             paste0(long_metrics$metric, " (-)"),
                                             long_metrics$metric))

  # Generate the plot with metric labels
  plot <- ggplot(long_metrics, aes(x = num_features, y = value)) +
    geom_line() +
    facet_wrap(~ metric_label, scales = "free_y") +
    theme_minimal() +
    labs(title = title,
         x = xAxis_text,
         y = "Metric value") +
    theme(plot.title = element_text(hjust = 0.5))

  # Save the plot as an EPS file
  ggsave(paste0("experiments/", file_save), plot = plot,
         width = 8, height = 6, device = "eps")

  return(plot)
}

gamma_values <- sapply(names(results_list), function(N) {
  results_list[[N]][["ii"]][["2"]][,"Connectivity"] # Extract gamma column
})

# Example usage of the function
plot_metrics_vs_num_features(
  results_list = results_list,  # The list containing internal metrics (ii) for different N values
  title = "External Metrics vs Number of Features",  # Plot title
  xAxis_text = "Number of Features (N)",  # X-axis label
  file_save = "metrics_vs_features_ALLAML.eps"  # Filename for saving the plot
)

# Define which metrics are better when higher (+) and which are better when lower (-)
positive_metrics <- c("calinski_harabasz", "dunn", "pbm", "tau", "gamma", "silhouette", "ensemble_ari")
negative_metrics <- c("c_index", "davies_bouldin", "mcclain_rao", "g_plus", "sd_dis", "ray_turi", "Compactness", "Connectivity")

# Extract internal metrics from results_list with proper character-based indexing
metrics_df_ii <- do.call(rbind, lapply(names(results_list), function(N) {
  df <- as.data.frame(results_list[[as.character(N)]][["ii"]][["2"]][1,])  # Ensure N is treated as a character
  df$num_features <- as.numeric(N)  # Convert N to numeric for plotting
  return(df)
}))

metrics_df_ei <- do.call(rbind, lapply(names(results_list), function(N) {
  df <- as.data.frame(results_list[[as.character(N)]][["ei"]][["2"]][1,])  # Ensure N is treated as a character
  df$num_features <- as.numeric(N)  # Convert N to numeric for plotting
  return(df)
}))

# Select relevant columns
df_silhouette <- metrics_df_ii[, c("num_features", "silhouette")]
df_accuracy <- metrics_df_ei[, c("num_features", "accuracy")]

# Merge both dataframes by `num_features`
metrics_df <- left_join(df_gamma, df_accuracy, by = "num_features")

metrics_df <- metrics_df %>%
  select(-Algorithms)

# Convert to long format for ggplot
long_metrics <- pivot_longer(metrics_df, cols = -num_features, names_to = "metric", values_to = "value")

# Plot both metrics
ggplot(long_metrics, aes(x = num_features, y = value, color = metric)) +
  geom_line(size = 1) +
  theme_minimal() +
  labs(title = "Gamma & Accuracy vs Number of Features",
       x = "Number of Features (N)",
       y = "Value",
       color = "Metric") +
  theme(plot.title = element_text(hjust = 0.5))


###
data(Meat)

UFS_methods <- list(
  "LLCFS"
)

# UFS_Results <- runUFS_manual_params(Meat$x, UFS_methods)
# save(UFS_Results, file = "experiments/UFS_Meat_LLCFS.RData")
load("experiments/UFS_Meat_LLCFS.RData")

method = "LLCFS"
print(method)
N <- 200
B <- 100
B.star=10
top_features <- UFS_Results$Results[[method]]$Result[[1]][1:N]

dice.obj <- dice(
  data = Meat$x[, top_features],
  #data = ALLAML$x,
  p.item = 0.85,
  reps = 10,
  nk = 5,
  algorithms = c("nmf"),
  nmf.method = c("lee"),
  progress = TRUE,
  verbose = FALSE,
  #ref.cl = as.integer(as.factor(Meat$y))
)

dice.obj$indices$ei

dice.obj$indices$ii

library(mclust)
adjustedRandIndex(dice.obj$clusters[,"CSPA"], as.integer(as.factor(Meat$y)))

cor_matrix <- cor(metrics_df_ii[, c("gamma", "calinski_harabasz", "silhouette", "davies_bouldin", "Compactness", "Connectivity")],
                  metrics_df_ei[, "accuracy"],
                  use = "complete.obs")
print(cor_matrix)

best_N <- metrics_df_ei$num_features[which.max(metrics_df_ei$accuracy)]
best_accuracy <- max(metrics_df_ei$accuracy)

cat("Best accuracy is at N =", best_N, "with a value of", best_accuracy, "\n")
library(ggplot2)

# Merge internal metrics with accuracy data
merged_df <- left_join(metrics_df_ii, metrics_df_ei, by = "num_features")

merged_df <- merged_df %>%
  select(-Algorithms.x,-Algorithms.y)

# Convert to long format
long_df <- pivot_longer(merged_df, cols = -num_features, names_to = "metric", values_to = "value")

# Plot internal metrics vs accuracy
ggplot(long_df, aes(x = num_features, y = value, color = metric)) +
  geom_line() +
  theme_minimal() +
  labs(title = "Accuracy & Internal Metrics vs Number of Features",
       x = "Number of Features (N)",
       y = "Metric Value",
       color = "Metric") +
  theme(plot.title = element_text(hjust = 0.5))


library(NMF)

# Aplicar NMF sobre la matriz de datos
nmf_result <- nmf(ALLAML$x, rank = 2, method = "brunet")

# Extraer la matriz W (representación de muestras en bases latentes)
W_matrix <- basis(nmf_result)

# Extraer la matriz H (peso de cada característica en cada base latente)
H_matrix <- coef(nmf_result)

# Ver coeficientes más altos
top_features <- apply(H_matrix, 2, function(x) order(x, decreasing = TRUE)[1:10])
print(top_features)



library(NMF)

# Aplicar NMF
res <- nmf(ALLAML$x, rank = 2, method = "brunet")
H <- coef(res)
colnames(H) <- paste0("Var", 1:ncol(H))  # Nombres automáticos: Var1, Var2, ...

colnames(ALLAML$x) <- paste0("Var", 1:ncol(ALLAML$x))  # Nombres automáticos: Var1, Var2, ...

N <- 150  # Número de características más importantes
top_features <- names(sort(apply(H, 2, max), decreasing = TRUE)[1:N])
top_features



for (N in seq(25, 300, by = 25)) {

  # Seleccionar las características superiores según el valor de N
  print(N)
  top_features <- names(sort(apply(H, 2, max), decreasing = TRUE)[1:N])
  # top_features <- UFS_Results$Results[[method]]$Result[[1]][1:N]

  # Ejecutar DICE con los top features seleccionados
  dice.obj <- dice(
    data = ALLAML$x[, top_features],
    nk = 2,
    p.item = 0.9,
    reps = 100,
    trim = TRUE,
    reweigh = TRUE,
    algorithms = c("nmf"),
    nmf.method = c("brunet"),
    cons.funs = c("CSPA"),
    progress = TRUE,
    verbose = FALSE,
    ref.cl = as.integer(as.factor(ALLAML$y))
  )

  # Almacenar los resultados en una lista para referencia posterior
  results_list[[as.character(N)]] <- list(
    ei = dice.obj$indices$ei,
    ii = dice.obj$indices$ii
  )
}

library(irlba)  # Para SVD eficiente en matrices grandes



calc_svd_entropy <- function(matrix, window_size = 5) {
  svd_result <- irlba(matrix, nv = window_size)
  singular_values <- svd_result$d
  total_energy <- sum(singular_values^2)
  v_j <- (singular_values^2) / total_energy
  entropy <- -sum(v_j * log(v_j)) / log(length(singular_values))
  return(entropy)
}

svd_entropy_selection <- function(data_matrix, n_features = 300, threshold = 0.90) {

  # 2. Clustering jerárquico de características
  dist_matrix <- dist(t(data_matrix), method = "euclidean")
  hclust_features <- hclust(dist_matrix, method = "ward.D2")
  feature_clusters <- cutree(hclust_features, k = n_features %/% 10)

  # 3. Selección por SVD-entropía
  selected_features <- c()
  for(cluster_id in unique(feature_clusters)) {
    cluster_data <- data_matrix[, feature_clusters == cluster_id]
    if(ncol(cluster_data) > 1) {
      entropy <- calc_svd_entropy(cluster_data)
      if(entropy < threshold) {
        selected_features <- c(selected_features,
                               colnames(cluster_data)[which.max(apply(cluster_data, 2, sd))])
      }
    }
  }

  return(head(unique(selected_features), n_features))
}

colnames(ALLAML$x) <- paste0("Gene_", 1:7129)

# Selección de características
genes_seleccionados <- svd_entropy_selection(ALLAML$x, n = 300)

# Verificar resultados
length(genes_seleccionados)  # Debería retornar 300 genes

# Clustering con PAM + Manhattan
pam_result <- pam(matriz_reducida,
                  metric = "manhattan",
                  stand = TRUE)  # Estandarización crucial

# Validación con silhouette
sil_score <- silhouette(pam_result$clustering, dist(matriz_reducida, "manhattan"))
mean(sil_score[,3])  # Valores >0.5 indican estructura sólida
