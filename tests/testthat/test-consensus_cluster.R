data(hgsc)
hgsc <- hgsc[1:40, 1:30]

test_that("No algorithms means all algorithms, output is an array", {
  skip_if_not_installed("apcluster")
  skip_if_not_installed("blockcluster")
  skip_if_not_installed("cluster")
  skip_if_not_installed("dbscan")
  skip_if_not_installed("e1071")
  skip_if_not_installed("kernlab")
  skip_if_not_installed("kohonen")
  x1 <- consensus_cluster(hgsc, nk = 4, reps = 1, progress = FALSE)
  expect_error(x1, NA)
  expect_true(is.array(x1))
})

test_that("Output can be saved with or without time in file name", {
  x1 <- consensus_cluster(hgsc, nk = 2:4, reps = 5, algorithms = "hc",
                          progress = FALSE, file.name = "CCOutput")
  x2 <- consensus_cluster(hgsc, nk = 2:4, reps = 5, algorithms = "hc",
                          progress = FALSE, file.name = "CCOutput",
                          time.saved = TRUE)
  expect_identical(x1, x2)
  file.remove(list.files(pattern = "CCOutput"))
})

test_that("Custom distance function can be passed", {
  skip_if_not_installed("apcluster")
  assign("my_dist", function(x) stats::dist(x, method = "manhattan"), pos = 1)
  x3 <- consensus_cluster(hgsc, nk = 2, reps = 5,
                          algorithms = c("nmf", "hc", "ap"),
                          distance = c("spear", "my_dist"), nmf.method = "lee",
                          progress = FALSE)
  expect_error(x3, NA)
})

test_that("Able to call only spearman distance", {
  x4a <- consensus_cluster(hgsc, nk = 2, reps = 5, algorithms = "hc",
                          distance = "spear", progress = FALSE)
  expect_error(x4a, NA)

  x4b <- consensus_cluster(hgsc, nk = 2, reps = 5, algorithms = "hc",
                           distance = "spear", abs = FALSE, progress = FALSE)
  expect_error(x4b, NA)
})

test_that("Data preparation on bootstrap samples works", {
  skip_if_not_installed("apcluster")
  x5 <- consensus_cluster(hgsc, nk = 3, reps = 3,
                          algorithms = c("nmf", "hc", "ap"), nmf.method = "lee",
                          prep.data = "sampled", progress = FALSE)
  expect_error(x5, NA)
})

test_that("no scaling means only choose complete cases and high signal vars", {
  x6 <- consensus_cluster(hgsc, nk = 2, reps = 2, algorithms = "hc",
                          scale = FALSE, progress = FALSE)
  expect_error(x6, NA)
})

test_that("t-SNE dimension reduction works", {
  x7 <- consensus_cluster(hgsc, nk = 4, reps = 1, algorithms = c("hc", "km"),
                          type = "tsne", progress = FALSE)
  expect_error(x7, NA)
})

test_that("Able to call Pearson distance", {
  x8 <- consensus_cluster(hgsc, nk = 2, reps = 5, algorithms = "hc",
                          distance = "pearson", progress = FALSE)
  expect_error(x8, NA)
})
