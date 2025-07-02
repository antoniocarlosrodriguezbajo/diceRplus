library(diceRplus)
library(mclust)
library(scales)
library(R.matlab)

algorithms = c("diana", "km", "pam", "ap", "sc",
               "gmm",  "cmeans", "hdbscan")


algorithms = c("gmm", "pam")



UFS_Methods <- list(
  "Inf-FS2020" = TRUE
)


mean_abs_correlation <- function(data) {
  cor_mat <- cor(data, use = "pairwise.complete.obs")
  round(mean(abs(cor_mat[upper.tri(data)])),4)
}

get_best_clustering <- function(data, nk, algorithms,
                                cons.funs = c("kmodes", "majority", "CSPA", "LCE", "LCA"),
                                n = 1,
                                p.item = 0.8,
                                reps = 10) {
  set.seed(123)
  result <- NULL
  best_clustering_alg <- NULL
  best_consensus_alg <- NULL
  best_clustering <- NULL

  tryCatch({
    dice.obj <- dice(
      data = data,
      nk = nk,
      p.item = p.item,
      rep = reps,
      algorithms = algorithms,
      cons.funs = cons.funs,
      trim = TRUE,
      n = n
    )

    best_clustering_alg <- dice.obj$indices$trim$alg.keep

    cc_data_list <- list()

    for (j in colnames(dice.obj$clusters)) {
      column <- dice.obj$clusters[, j]
      cat("Method:", j, "\n")
      num_labels <- length(column)
      cc_data <- array(column, dim = c(num_labels, 1, 1, 1))
      row_names <- rownames(column)
      dimnames(cc_data) <- list(
        row_names,
        "R1",
        j,
        as.character(nk)
      )
      cc_data_list[[j]] <- cc_data
    }

    result <- do.call(consensus_evaluate, c(list(data), cc_data_list, list(n = 1, trim = TRUE)))
    best_consensus_alg <- result$trim.obj$alg.keep
    best_clustering <- as.vector(result$trim.obj$E.new[[1]])
  },
  error = function(e) {
    message("Error: ", conditionMessage(e))
  })

  return(list(
    best_clustering_alg = best_clustering_alg,
    best_consensus_alg = best_consensus_alg,
    best_clustering = best_clustering
  ))
}


get_best_clustering_simple <- function(data, nk, algorithms) {

  set.seed(123)

  cc = consensus_cluster(data, nk = nk, reps = 1, p.item = 1, algorithms = algorithms, progress = TRUE)
  eval_result = consensus_evaluate(data, cc, trim = TRUE, n = 1)
  best_clustering_alg = eval_result$trim.obj$alg.keep
  best_clustering = as.vector(eval_result$trim.obj$E.new[[1]])

  return(list(
    best_clustering_alg = best_clustering_alg,
    best_clustering = best_clustering
  ))
}

get_best_UFS_reduction <- function(data, nk, algorithms) {
  setVariable(matlab, X = data)

  evaluate(matlab, "whos")
  evaluate(matlab, "rng('default');")
  evaluate(matlab, "rng(42);")

  result_list <- list()  # Store top_features and best_clustering by alpha

  for (alpha in seq(0, 1, by = 0.1)) {
    evaluate(matlab, "clear RANKED WEIGHT SUBSET;")

    # Run InfFS_U for current alpha
    cmd <- sprintf('[RANKED, WEIGHT, SUBSET] = InfFS_U(X, %.1f);', alpha)
    evaluate(matlab, cmd)

    ranked  <- getVariable(matlab, "RANKED")
    weight  <- getVariable(matlab, "WEIGHT")
    subset  <- getVariable(matlab, "SUBSET")

    top_features <- as.vector(subset$SUBSET)
    num_features <- length(top_features)

    cat("Alpha =", alpha, "| Selected features:", num_features, "\n")

    # Clustering
    cc <- consensus_cluster(data[, top_features], nk = nk, reps = 1, p.item = 1,
                            algorithms = algorithms, progress = TRUE)
    eval_result <- consensus_evaluate(data[, top_features], cc, trim = TRUE, n = 1)
    best_clustering_alg <- eval_result$trim.obj$alg.keep
    best_clustering <- as.vector(eval_result$trim.obj$E.new[[1]])

    # Store in list
    result_list[[as.character(alpha)]] <- list(
      top_features = top_features,
      best_clustering = best_clustering,
      best_clustering_alg = best_clustering_alg
    )
  }

  return(result_list)
}




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


########################################################
########################################################
# Meat
########################################################
########################################################

data(Meat)

data <- Meat

mean(data$x)
sd(data$x)
table(data$y)
mean_abs_correlation(data$x)




result <- get_best_clustering(data = data$x,
                              nk = 5,
                              algorithms = algorithms)

paste(result$best_clustering_alg, result$best_consensus_alg,  sep = "-")
adjustedRandIndex(result$best_clustering, as.numeric(data$y))



# Set the dataset variable
setVariable(matlab, X = data$x)

# Evaluate MATLAB command to display variables
evaluate(matlab, "whos")

evaluate(matlab, "rng('default');")
evaluate(matlab, "rng(42);")

for (alpha in seq(0, 1, by = 0.1)) {

  evaluate(matlab, "clear RANKED WEIGHT SUBSET;")


  # Run the feature selection method in MATLAB
  cmd <- sprintf('[RANKED, WEIGHT, SUBSET] = InfFS_U(X, %.1f);', alpha)
  evaluate(matlab, cmd)

  ranked  <- getVariable(matlab, "RANKED")
  weight  <- getVariable(matlab, "WEIGHT")
  subset  <- getVariable(matlab, "SUBSET")

  top_features <- as.vector(subset$SUBSET)
  num_features = length(top_features)

  result <- get_best_clustering(data = data$x[,top_features],
                                nk = 5,
                                algorithms = "diana",
                                cons.funs = "LCE")
  print(paste(
    alpha,
    num_features,
    result$best_clustering_alg, result$best_consensus_alg,
    sep = "-"))

  print(adjustedRandIndex(result$best_clustering, as.numeric(data$y)))
  conf_mat <- ev_confmat(result$best_clustering, as.numeric(data$y))
  print(conf_mat$.estimate[1])
}


########################################################
########################################################
# Lymphoma
########################################################
########################################################

data(lymphoma)

data <- lymphoma

mean(data$x)
sd(data$x)
table(data$y)
mean_abs_correlation(data$x)

data$y <- data$y +1


result <- get_best_clustering(data = data$x,
                              nk = 3,
                              algorithms = algorithms)

paste(result$best_clustering_alg, result$best_consensus_alg,  sep = "-")
adjustedRandIndex(result$best_clustering, as.numeric(data$y))


# Set the dataset variable
setVariable(matlab, X = data$x)

# Evaluate MATLAB command to display variables
evaluate(matlab, "whos")

evaluate(matlab, "rng('default');")
evaluate(matlab, "rng(42);")

for (alpha in seq(0, 1, by = 0.1)) {

  evaluate(matlab, "clear RANKED WEIGHT SUBSET;")


  # Run the feature selection method in MATLAB
  cmd <- sprintf('[RANKED, WEIGHT, SUBSET] = InfFS_U(X, %.1f);', alpha)
  evaluate(matlab, cmd)

  ranked  <- getVariable(matlab, "RANKED")
  weight  <- getVariable(matlab, "WEIGHT")
  subset  <- getVariable(matlab, "SUBSET")

  top_features <- as.vector(subset$SUBSET)
  num_features = length(top_features)

  result <- get_best_clustering(data = data$x[,top_features],
                                nk = 3,
                                algorithms = "hc",
                                cons.funs = "CSPA")
  print(paste(
    alpha,
    num_features,
    result$best_clustering_alg, result$best_consensus_alg,
    sep = "-"))

  print(adjustedRandIndex(result$best_clustering, as.numeric(data$y)))

}



########################################################
########################################################
# Prostate GE
########################################################
########################################################

data(prostate_ge)

data <- prostate_ge

mean(data$x)
sd(data$x)

table(data$y)
mean_abs_correlation(data$x)


# setdiff(algorithms, "som")

result <- get_best_clustering(data = data$x,
                              nk = 2,
                              algorithms = setdiff(algorithms, "som"))

paste(result$best_clustering_alg, result$best_consensus_alg,  sep = "-")
adjustedRandIndex(result$best_clustering, as.numeric(data$y))

conf_mat <- ev_confmat(result$best_clustering, as.numeric(data$y))
conf_mat$.estimate[1]

# > paste(result$best_clustering_alg, result$best_consensus_alg,  sep = "-")
# [1] "KM-majority"
# > adjustedRandIndex(result$best_clustering, as.numeric(data$y))
# [1] 0.02304851
# > conf_mat <- ev_confmat(result$best_clustering, as.numeric(data$y))
# > conf_mat$.estimate[1]
# [1] 0.5882353


# Set the dataset variable
setVariable(matlab, X = data$x)

# Evaluate MATLAB command to display variables
evaluate(matlab, "whos")

evaluate(matlab, "rng('default');")
evaluate(matlab, "rng(42);")

for (alpha in seq(0, 1, by = 0.1)) {

  evaluate(matlab, "clear RANKED WEIGHT SUBSET WEIGHT_SUM WEIGHT_MEAN;")


  # Run the feature selection method in MATLAB
  cmd <- sprintf('[RANKED, WEIGHT, SUBSET] = InfFS_U(X, %.1f);', alpha)
  evaluate(matlab, cmd)

  ranked  <- getVariable(matlab, "RANKED")
  weight  <- getVariable(matlab, "WEIGHT")
  subset  <- getVariable(matlab, "SUBSET")

  top_features <- as.vector(subset$SUBSET)
  num_features = length(top_features)

  subset_weights <- weight$WEIGHT[top_features]

  result <- get_best_clustering(data = data$x[,top_features],
                                nk = 2,
                                algorithms = "pam",
                                cons.funs = "majority")
  print(paste(
    alpha,
    num_features,
    result$best_clustering_alg, result$best_consensus_alg,
    sep = "-"))

  print(adjustedRandIndex(result$best_clustering, as.numeric(data$y)))
  conf_mat <- ev_confmat(result$best_clustering, as.numeric(data$y))
  print(conf_mat$.estimate[1])

}


########################################################
########################################################
# lung_cancer
########################################################
########################################################

data(lung_cancer)

data <- lung_cancer

mean(data$x)
sd(data$x)

# data$x<- prepare_data(data$x)
ncol(data$x)

table(data$y)

mean_abs_correlation(data$x)
# 0.1467


result <- get_best_clustering(data = data$x,
                              nk = 5,
                              n = 1,
                              algorithms = algorithms)

paste(result$best_clustering_alg, result$best_consensus_alg,  sep = "-")
adjustedRandIndex(result$best_clustering, as.numeric(data$y))

conf_mat <- ev_confmat(result$best_clustering, as.numeric(data$y))
conf_mat$.estimate[1]


# Set the dataset variable
setVariable(matlab, X = data$x)

# Evaluate MATLAB command to display variables
evaluate(matlab, "whos")

evaluate(matlab, "rng('default');")
evaluate(matlab, "rng(42);")

for (alpha in seq(0, 1, by = 0.1)) {

  evaluate(matlab, "clear RANKED WEIGHT SUBSET WEIGHT_SUM WEIGHT_MEAN;")


  # Run the feature selection method in MATLAB
  cmd <- sprintf('[RANKED, WEIGHT, SUBSET] = InfFS_U(X, %.1f);', alpha)
  evaluate(matlab, cmd)

  ranked  <- getVariable(matlab, "RANKED")
  weight  <- getVariable(matlab, "WEIGHT")
  subset  <- getVariable(matlab, "SUBSET")

  top_features <- as.vector(subset$SUBSET)
  num_features = length(top_features)

  subset_weights <- weight$WEIGHT[top_features]

  result <- get_best_clustering(data = data$x[,top_features],
                                nk = 5,
                                algorithms = "gmm",
                                cons.funs = "LCE")
  print(paste(
    alpha,
    num_features,
    result$best_clustering_alg, result$best_consensus_alg,
    sep = "-"))

  print(adjustedRandIndex(result$best_clustering, as.numeric(data$y)))
  conf_mat <- ev_confmat(result$best_clustering, as.numeric(data$y))
  print(conf_mat$.estimate[1])

}


########################################################
########################################################
# warpPIE10P
########################################################
########################################################

data(warpPIE10P)

data <- warpPIE10P

mean(data$x)
sd(data$x)
min(data$x)
max(data$x)

data$x <- apply(data$x, 2, rescale)

mean(data$x)
sd(data$x)
min(data$x)
max(data$x)

mean_abs_correlation(data$x)


result <- get_best_clustering(data = data$x,
                              nk = 10,
                              n = 1,
                              algorithms = algorithms)

paste(result$best_clustering_alg, result$best_consensus_alg,  sep = "-")
adjustedRandIndex(result$best_clustering, as.numeric(data$y))

conf_mat <- ev_confmat(result$best_clustering, as.numeric(data$y))
conf_mat$.estimate[1]


# Set the dataset variable
setVariable(matlab, X = data$x)

# Evaluate MATLAB command to display variables
evaluate(matlab, "whos")

evaluate(matlab, "rng('default');")
evaluate(matlab, "rng(42);")

for (alpha in seq(0, 1, by = 0.1)) {

  evaluate(matlab, "clear RANKED WEIGHT SUBSET WEIGHT_SUM WEIGHT_MEAN;")


  # Run the feature selection method in MATLAB
  cmd <- sprintf('[RANKED, WEIGHT, SUBSET] = InfFS_U(X, %.1f);', alpha)
  evaluate(matlab, cmd)

  ranked  <- getVariable(matlab, "RANKED")
  weight  <- getVariable(matlab, "WEIGHT")
  subset  <- getVariable(matlab, "SUBSET")

  top_features <- as.vector(subset$SUBSET)
  num_features = length(top_features)

  subset_weights <- weight$WEIGHT[top_features]

  result <- get_best_clustering(data = data$x[,top_features],
                                nk = 10,
                                algorithms = "km",
                                cons.funs = "majority")
  print(paste(
    alpha,
    num_features,
    result$best_clustering_alg, result$best_consensus_alg,
    sep = "-"))

  print(adjustedRandIndex(result$best_clustering, as.numeric(data$y)))
  conf_mat <- ev_confmat(result$best_clustering, as.numeric(data$y))
  print(conf_mat$.estimate[1])

}


cc <- consensus_cluster(data$x, nk = 10, reps = 1, p.item = 1, algorithms = algorithms, progress = FALSE)

eval_result <- consensus_evaluate(data$x, cc, trim = TRUE, n = 1)

eval_result$trim.obj$alg.keep

########################################################
########################################################
# SIMPLE SIMPLE
########################################################
########################################################
########################################################
########################################################
# SIMPLE SIMPLE
########################################################
########################################################


########################################################
########################################################
# Meat
########################################################
########################################################

data(Meat)

data <- Meat

mean(data$x)
sd(data$x)
table(data$y)
mean_abs_correlation(data$x)

result <- get_best_clustering_simple(data$x, nk = 5, algorithms = algorithms)

paste(result$best_clustering_alg,  sep = "-")
adjustedRandIndex(result$best_clustering, as.numeric(data$y))
conf_mat <- ev_confmat(result$best_clustering, as.numeric(data$y))
conf_mat$.estimate[1]

result_UFS <- get_best_UFS_reduction(data$x, nk = 5, algorithms = 'km')

for (alpha in names(result_UFS)) {
  info <- result_UFS[[alpha]]

  cat("Alpha:", alpha,
      "| Features selected:", length(info$top_features),
      "| Best algorithm:", info$best_clustering_alg, "\n")

  print(adjustedRandIndex(info$best_clustering, as.numeric(data$y)))
  conf_mat <- ev_confmat(info$best_clustering, as.numeric(data$y))
  print(conf_mat$.estimate[1])

}

########################################################
########################################################
# Lymphoma
########################################################
########################################################

data(lymphoma)

data <- lymphoma

mean(data$x)
sd(data$x)
table(data$y)
mean_abs_correlation(data$x)

data$y <- data$y +1

result <- get_best_clustering_simple(data$x, nk = 3, algorithms = algorithms)

paste(result$best_clustering_alg,  sep = "-")
adjustedRandIndex(result$best_clustering, as.numeric(data$y))
conf_mat <- ev_confmat(result$best_clustering, as.numeric(data$y))
conf_mat$.estimate[1]

result_UFS <- get_best_UFS_reduction(data$x, nk = 3, algorithms = 'km')

for (alpha in names(result_UFS)) {
  info <- result_UFS[[alpha]]

  cat("Alpha:", alpha,
      "| Features selected:", length(info$top_features),
      "| Best algorithm:", info$best_clustering_alg, "\n")

  print(adjustedRandIndex(info$best_clustering, as.numeric(data$y)))
  conf_mat <- ev_confmat(info$best_clustering, as.numeric(data$y))
  print(conf_mat$.estimate[1])

}


########################################################
########################################################
# Prostate GE
########################################################
########################################################

data(prostate_ge)

data <- prostate_ge

mean(data$x)
sd(data$x)

table(data$y)
mean_abs_correlation(data$x)

result <- get_best_clustering_simple(data$x, nk = 2, algorithms = algorithms)

paste(result$best_clustering_alg,  sep = "-")
adjustedRandIndex(result$best_clustering, as.numeric(data$y))
conf_mat <- ev_confmat(result$best_clustering, as.numeric(data$y))
conf_mat$.estimate[1]

result_UFS <- get_best_UFS_reduction(data$x, nk = 2, algorithms = 'km')

for (alpha in names(result_UFS)) {
  info <- result_UFS[[alpha]]

  cat("Alpha:", alpha,
      "| Features selected:", length(info$top_features),
      "| Best algorithm:", info$best_clustering_alg, "\n")

  print(adjustedRandIndex(info$best_clustering, as.numeric(data$y)))
  conf_mat <- ev_confmat(info$best_clustering, as.numeric(data$y))
  print(conf_mat$.estimate[1])

}

########################################################
########################################################
# lung_cancer
########################################################
########################################################

data(lung_cancer)

data <- lung_cancer

mean(data$x)
sd(data$x)

# data$x<- prepare_data(data$x)
ncol(data$x)

table(data$y)

mean_abs_correlation(data$x)
# 0.1467

result <- get_best_clustering_simple(data$x, nk = 5, algorithms = "gmm")

paste(result$best_clustering_alg,  sep = "-")
adjustedRandIndex(result$best_clustering, as.numeric(data$y))
conf_mat <- ev_confmat(result$best_clustering, as.numeric(data$y))
conf_mat$.estimate[1]

result_UFS <- get_best_UFS_reduction(data$x, nk = 5, algorithms = 'km')

for (alpha in names(result_UFS)) {
  info <- result_UFS[[alpha]]

  cat("Alpha:", alpha,
      "| Features selected:", length(info$top_features),
      "| Best algorithm:", info$best_clustering_alg, "\n")

  print(adjustedRandIndex(info$best_clustering, as.numeric(data$y)))
  conf_mat <- ev_confmat(info$best_clustering, as.numeric(data$y))
  print(conf_mat$.estimate[1])

}


########################################################
########################################################
# warpPIE10P
########################################################
########################################################
data(warpPIE10P)

data <- warpPIE10P

mean(data$x)
sd(data$x)
min(data$x)
max(data$x)

data$x <- apply(data$x, 2, rescale)

data$x <- scale(data$x)

mean(data$x)
sd(data$x)
min(data$x)
max(data$x)


result <- get_best_clustering_simple(data$x, nk = 10, algorithms = algorithms)

paste(result$best_clustering_alg,  sep = "-")
adjustedRandIndex(result$best_clustering, as.numeric(data$y))
conf_mat <- ev_confmat(result$best_clustering, as.numeric(data$y))
conf_mat$.estimate[1]

result_UFS <- get_best_UFS_reduction(data$x, 10, 'km')

for (alpha in names(result_UFS)) {
  info <- result_UFS[[alpha]]

  cat("Alpha:", alpha,
      "| Features selected:", length(info$top_features),
      "| Best algorithm:", info$best_clustering_alg, "\n")

  print(adjustedRandIndex(info$best_clustering, as.numeric(data$y)))
  conf_mat <- ev_confmat(info$best_clustering, as.numeric(data$y))
  print(conf_mat$.estimate[1])

}



########################################################
########################################################
# END
########################################################
########################################################



# Close MATLAB connection
close(matlab)

