# Anderlucci's meat with UFS integrated with diceRplus
# Anderlucci's: https://rdrr.io/cran/RPEClust/f/
# UFS: https://github.com/farhadabedinzadeh/AutoUFSTool

library(tidyr)
library(ggplot2)
library(dplyr)
library(diceRplus)


# Calcute internal_metrics based on consensus_evaluate of diceR
calculate_internal_metrics <- function(data,cluster_labels) {
  num_labels <- length(cluster_labels)

  # Create structure for consensus_evaluate
  cc_data <- array(cluster_labels, dim = c(num_labels, 1, 1, 1))
  row_names <- rownames(data)
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

plot_metrics_vs_num_features2 <- function(data1, legend_data1,
                                          data2, legend_data2,
                                          legend_title,
                                          internal_metrics,
                                          title, xAxis_text, file_save) {
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
  #long_metrics <- long_metrics[long_metrics$metric != "s_dbw", ]
  long_metrics <- long_metrics[long_metrics$metric %in% c("silhouette", "davies_bouldin", "calinski_harabasz", "c_index"), ]


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

  # Save plot
  ggsave(paste0("experiments/", file_save), plot = plot,
         width = 8, height = 6, device = "eps")

  return(plot)
}

# Meat dataset
data(Meat)

############################################################
# Experiment in Anderlucci's paper (original implementation)
############################################################
B <- 1000
B.star=100

# Run RPGMMClu
execution_time <- system.time(out.clu_baseline <- RPGMMClu(Meat$x,
                                                           Meat$y,
                                                           g=5,
                                                           B=B,
                                                           B.star=B.star,
                                                           verb=TRUE))["elapsed"]

# Calculate internal metrics
internal_metrics <- calculate_internal_metrics(Meat$x,
                                               out.clu_baseline$ensemble$label.vec)
# Log the experiment
exp_data <- experiment_logger(
  description = "Clustering with RPGMMClu - Experiment in Anderlucci's paper",
  dataset = "Meat",
  ensemble_method = "RPGMMClu",
  ensemble_method_params = list(g = 5, B = B, B.star = B.star),
  execution_time = as.numeric(execution_time),
  labels_clustering = out.clu_baseline$ensemble$label.vec,
  internal_metrics = internal_metrics,
  external_metrics = list(ensemble_ari = out.clu_baseline$ensemble$ari[[1]])
)
save_experiment(exp_data)

############################################################
# Experiment in Anderlucci's paper (parallel implementation)
############################################################
B <- 1000
B.star=100

# Run RPGMMClu_parallel
execution_time <- system.time(out.clu_baseline_p <- RPGMMClu_parallel(Meat$x,
                                                                      Meat$y,
                                                                      g=5,
                                                                      B=B,
                                                                      B.star=B.star))["elapsed"]

# Calculate internal metrics
internal_metrics_p <- calculate_internal_metrics(Meat$x,
                                               out.clu_baseline_p$ensemble$label.vec)
# Log the experiment
exp_data <- experiment_logger(
  description = "Clustering with RPGMMClu - Experiment in Anderlucci's paper (parallel implementation)",
  dataset = "Meat",
  ensemble_method = "RPGMMClu_parallel",
  ensemble_method_params = list(g = 5, B = B, B.star = B.star),
  execution_time = as.numeric(execution_time),
  labels_clustering = out.clu_baseline_p$ensemble$label.vec,
  internal_metrics = internal_metrics_p,
  external_metrics = list(ensemble_ari = out.clu_baseline_p$ensemble$ari[[1]])
)
save_experiment(exp_data)


############################################################
# Experiment with UFS
############################################################
# UFS Methods list
# If a method is set to `TRUE`,
# the AutoHotkey script will send `1` to MATLAB as input for default parameters.
UFS_Methods <- list(
  "InfFS" = FALSE,
  "Laplacian" = TRUE,
  "MCFS" = TRUE,
  "LLCFS" = TRUE,
  "CFS" = FALSE,
  "FSASL" = TRUE,
  "DGUFS" = TRUE,
  "UFSOL" = TRUE,
  "SPEC" = TRUE,
  # "SOCFS" = TRUE, # Won't work
  "SOGFS" = TRUE,
  "UDFS" = TRUE,
  "SRCFS" = TRUE,
  "FMIUFS" = TRUE,
  "UAR_HKCMI" = TRUE,
  "RNE" = TRUE,
  # "FRUAFR" = TRUE, # Won't work
  "U2FS" = TRUE,
  "RUFS" = TRUE,
  "NDFS" = TRUE,
  "EGCFS" = TRUE,
  "CNAFS" = TRUE,
  "Inf-FS2020" = TRUE
)
# Run just once
# UFS_Results <- runUFS(Meat$x, UFS_Methods)
# save(UFS_Results, file = "experiments/UFS_Meat.RData")

load("experiments/UFS_Meat.RData")

############################################################
# Experiments with InfFS and N most relevant features
# (B,B*):
#        (1000,100)
#        (500,50)
#        (100,10)
############################################################
B <- 100
B.star <- 10

# Fixed UFS: InfFS
method <- names(UFS_Results$Results)[1]
print(method)

# Iterate over different values of N (from 20 to 100 in steps of 10)
for (N in seq(25, ncol(Meat$x), by = 25)) {
  print(N)
  top_features <- UFS_Results$Results[[method]]$Result[[1]][1:N]
  Meat$x_filtered <- Meat$x[, top_features]

  # Run RPGMMClu_parallel
  execution_time <- system.time(out.clu_baseline_UFS <- RPGMMClu_parallel(Meat$x_filtered,
                                                                          Meat$y,
                                                                          g = 5,
                                                                          B = B,
                                                                          B.star = B.star))["elapsed"]
  # Calculate internal metrics
  internal_metrics_UFS <- calculate_internal_metrics(Meat$x_filtered,
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
    internal_metrics = internal_metrics_UFS,
    external_metrics = list(ensemble_ari = out.clu_baseline_UFS$ensemble$ari[[1]])
  )

  save_experiment(exp_data)
}


## Analysis (infFS)

experiments_data_all <- load_experiments()

# B=100, B*=10
filter_100 <- experiments_data_all$experiment_id >= 206 &
  experiments_data_all$experiment_id <= 247

# B=500, B*=50
filter_500 <- experiments_data_all$experiment_id >= 64 &
              experiments_data_all$experiment_id <= 105

# B=1000, B*=100
filter_1000 <- experiments_data_all$experiment_id >= 22 &
  experiments_data_all$experiment_id <= 63

experiments_data <- experiments_data_all[filter_100 | filter_500 | filter_1000,]

# Define which metrics should be minimized
metrics_to_minimize <- c("c_index", "davies_bouldin", "mcclain_rao", "g_plus",
                         "s_dbw", "ray_turi", "Compactness", "Connectivity")

# Extract the best internal metrics based on whether they should be maximized or minimized
best_metrics <- lapply(names(experiments_data$internal_metrics[[1]]), function(metric) {
  values <- sapply(experiments_data$internal_metrics, function(x) x[[metric]])

  # Determine whether to find the maximum or minimum value
  if (metric %in% metrics_to_minimize) {
    best_row <- which.min(values)  # If it's a metric to minimize, find the minimum
  } else {
    best_row <- which.max(values)  # Otherwise, find the maximum
  }

  experiment_id <- experiments_data$experiment_id[best_row]

  list(metric = metric, experiment_id = experiment_id, best_value = values[best_row])
})

# Convert to a data frame for better visualization
df_results <- do.call(rbind, lapply(best_metrics, as.data.frame))

# Extract additional columns from experiments_data
df_results$ensemble_method_params <- sapply(df_results$experiment_id, function(exp_id) experiments_data$ensemble_method_params[experiments_data$experiment_id == exp_id])
df_results$num_features <- sapply(df_results$experiment_id, function(exp_id) experiments_data$num_features[experiments_data$experiment_id == exp_id])

df_results <- df_results[order(df_results$metric), ]
print(df_results)


# Count occurrences of the best experiments across all metric
best_experiment_counts <- table(sapply(best_metrics, function(x) x$experiment_id))

# Convert to a data frame and sort by descending frequency
ranking_df <- as.data.frame(best_experiment_counts)
colnames(ranking_df) <- c("experiment_id", "best_metric_count")
ranking_df <- ranking_df[order(ranking_df$best_metric_count, decreasing = TRUE), ]

# Extract additional columns based on experiment_id
ranking_df$ensemble_method_params <- sapply(ranking_df$experiment_id, function(exp_id) experiments_data$ensemble_method_params[experiments_data$experiment_id == exp_id])
ranking_df$num_features <- sapply(ranking_df$experiment_id, function(exp_id) experiments_data$num_features[experiments_data$experiment_id == exp_id])
ranking_df$ensemble_ari <- sapply(ranking_df$experiment_id, function(exp_id) experiments_data$external_metrics[[which(experiments_data$experiment_id == exp_id)]]$ensemble_ari)

# Print the updated ranking
print(ranking_df)


## Experiment 66 holds the top position

# Graphics

plotInfFS_internal <- plot_metrics_vs_num_features2(experiments_data_all[filter_500,], "B=500,B*=50",
                                                    experiments_data_all[filter_100,], "B=100,B*=10",
                                                    "RPEClu params",
                                                   internal = TRUE,
                                                   "Number of selected features (Inf-FS) vs internal metrics",
                                                   "Number of features (meat dataset)",
                                                   "infFS_meat_internal.eps")
print(plotInfFS_internal)


plotInfFS_external <- plot_metrics_vs_num_features2(experiments_data_all[filter_500,], "B=500,B*=50",
                                                    experiments_data_all[filter_100,], "B=100,B*=10",
                                                   internal = FALSE,
                                                   "RPEClu params",
                                                   "Number of selected features (Inf-FS) vs ARI",
                                                   "Number of features (meat dataset)",
                                                   "infFS_meat_external.eps")
print(plotInfFS_external)

# Add baseline
plotInfFS_external <- plotInfFS_external +
  geom_hline(yintercept = 0.32, linetype = "dashed", color = "#1E3A5F") +  # Línea de referencia
  annotate("text", x = min(long_metrics$num_features), y = 0.32,
           label = "Baseline = 0.32", hjust = 0, vjust = -0.5, color = "#1E3A5F", size = 2) +  # Texto más pequeño
  facet_wrap(~ metric_label, scales = "free_y")

print(plotInfFS_external)
ggsave("experiments/infFS_meat_external.eps",
       plot = plotInfFS_external, device = "eps")


plotInfFS_external <- plotInfFS_external +
  theme(text = element_text(size = 6),  # Reduce el tamaño del texto
        axis.text = element_text(size = 5),  # Ajusta las etiquetas de los ejes
        plot.title = element_text(size = 7))  # Modifica el tamaño del título
plotInfFS_external <- plotInfFS_external +
  theme(legend.position = "bottom",  # Mantiene la leyenda abajo
        legend.margin = margin(t = -10, r = 0, b = 0, l = 0))  # Reduce el espacio superior
ggsave("experiments/infFS_meat_external.eps",
       plot = plotInfFS_external, device = "eps",
       width = 5, height = 3, scale = 0.6)  # Reduce la escala del gráfico

# Execution times
UFS_excution_time <- UFS_Results$ExecutionTimes[[method]]
infFS_500_execution_time <- sum(experiments_data_all[filter_500,"execution_time"])
infFS_1000_execution_time <- sum(experiments_data_all[filter_1000,"execution_time"])
print(round(UFS_excution_time))
print(round(infFS_500_execution_time))
print(round(infFS_1000_execution_time))


