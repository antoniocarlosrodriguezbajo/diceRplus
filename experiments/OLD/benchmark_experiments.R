# Check GPU implementation
# https://pypi.org/project/gmm-gpu/
# https://mzdravkov.com/docs/gmm_gpu/index.html


library(diceRplus)
data(lymphoma)


cor_mat <- cor(lymphoma$x, use = "pairwise.complete.obs")
mean_cor <- mean(abs(cor_mat[upper.tri(cor_mat)]))

var_vec <- apply(lymphoma$x, 2, var, na.rm = TRUE)
low_var_prop <- mean(var_vec < quantile(var_vec, 0.1))

pca <- prcomp(lymphoma$x, scale. = TRUE)
var_exp <- summary(pca)$importance[2, ]
cum_var_90 <- which(cumsum(var_exp) > 0.9)[1]

ufs_candidates <- c()

if(mean_cor > 0.5) ufs_candidates <- c(ufs_candidates, "MCFS", "Inf-FS", "LaplacianScore")
if(low_var_prop > 0.3) ufs_candidates <- c(ufs_candidates, "UDFS", "NDFS", "Inf-FS")
if(cum_var_90 < ncol(lymphoma$x) * 0.2) ufs_candidates <- c(ufs_candidates, "UDFS", "NDFS", "Inf-FS")
if(ncol(lymphoma$x) > 3000) ufs_candidates <- c(ufs_candidates, "LaplacianScore", "SPEC", "Inf-FS")


# Asegúrate de cubrir al menos 3 familias
ufs_candidates <- unique(c(ufs_candidates, "UDFS", "Inf-FS", "LaplacianScore", "NDFS"))
ufs_candidates <- unique(ufs_candidates)
print(ufs_candidates)

UFS_Methods <- list(
  "InfFS" = FALSE,
  "Laplacian" = TRUE,
  # "MCFS" = TRUE,
  "LLCFS" = TRUE,
  "CFS" = FALSE,
  "FSASL" = TRUE,
  "DGUFS" = TRUE,
  "UFSOL" = TRUE,
  "SPEC" = TRUE,
  # "SOCFS" = TRUE, # Won't work
  # "SOGFS" = TRUE,
  "UDFS" = TRUE,
  "SRCFS" = TRUE,
  # "FMIUFS" = TRUE,
  # "UAR_HKCMI" = TRUE,
  "RNE" = TRUE,
  # "FRUAFR" = TRUE, # Won't work
  # "U2FS" = TRUE,
  "RUFS" = TRUE,
  "NDFS" = TRUE,
  "EGCFS" = TRUE,
  "CNAFS" = TRUE,
  "Inf-FS2020" = TRUE
)

# Run just once
# UFS_lymphoma <- runUFS(lymphoma$x, UFS_Methods)
# save(UFS_lymphoma, file = "experiments/UFS_lymphoma.RData")

load("experiments/UFS_Results_lymphoma.RData")

# Fixed UFS: InfFS
method <- names(UFS_lymphoma$Results)[5]
print(method)
N=2000
N=1000
top_features <- UFS_lymphoma$Results[[method]]$Result[[1]][1:N]
lymphoma$x_filtered <- lymphoma$x[, top_features]

B <- 100
B.star=10

# Run RPGMMClu
execution_time <- system.time(out.clu_baseline <- RPGMMClu_parallel(lymphoma$x[, top_features],
                                                           lymphoma$y,
                                                           g=3,
                                                           B=B,
                                                           B.star=B.star,
                                                           verb=TRUE))["elapsed"]


# Calculate internal metrics
internal_metrics_p <- calculate_internal_metrics(lymphoma$x[, top_features],
                                                 out.clu_baseline$ensemble$label.vec)
silhouette_values <- cluster::silhouette(out.clu_baseline$ensemble$label.vec, dist(lymphoma$x[, top_features]))
silhouette_score <- mean(silhouette_values[, 3])
print(silhouette_score)



data(prostate_ge)
dim(prostate_ge$x)
min(apply(prostate_ge$x, 2, sd))

# Run just once
#UFS_prostate_ge <- runUFS(prostate_ge$x, UFS_Methods)
#save(UFS_prostate_ge, file = "experiments/UFS_prostate_ge.RData")

load("experiments/UFS_prostate_ge.RData")

# Fixed UFS: InfFS
method <- names(UFS_prostate_ge$Results)[1]
print(method)
N=2000
N=1000
top_features <- UFS_prostate_ge$Results[[method]]$Result[[1]][1:N]
prostate_ge$x_filtered <- prostate_ge$x[, top_features]

B <- 100
B.star=10

# Run RPGMMClu
execution_time <- system.time(out.clu_baseline <- RPGMMClu_parallel(prostate_ge$x[, top_features],
                                                                    prostate_ge$y,
                                                                    g=2,
                                                                    B=B,
                                                                    B.star=B.star,
                                                                    verb=TRUE))["elapsed"]


data(ALLAML)
dim(ALLAML$x)
min(apply(ALLAML$x, 2, sd))

# Run just once
#UFS_ALLAML <- runUFS(ALLAML$x, UFS_Methods)
#save(UFS_ALLAML, file = "experiments/UFS_ALLAML.RData")

load("experiments/UFS_ALLAML.RData")

# Fixed UFS: InfFS
method <- names(UFS_ALLAML$Results)[5]
print(method)
N=2000
N=1000
top_features <- UFS_ALLAML$Results[[method]]$Result[[1]][1:N]
ALLAML$x_filtered <- ALLAML$x[, top_features]

B <- 200
B.star=20

# Run RPGMMClu
execution_time <- system.time(out.clu_baseline <- RPGMMClu_parallel(ALLAML$x_filtered,
                                                                    ALLAML$y,
                                                                    g=3,
                                                                    B=B,
                                                                    B.star=B.star,
                                                                    verb=TRUE))["elapsed"]
############################################################
# COIL20
############################################################

data(COIL20)
dim(COIL20$x)
