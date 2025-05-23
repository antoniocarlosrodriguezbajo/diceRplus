data(hgsc)
dat <- hgsc[1:40, 1:30]
k <- 4
x <- consensus_cluster(dat, nk = k, reps = 4, progress = FALSE,
                       algorithms = c("nmf", "hc", "diana"), nmf.method = "lee")
x_imputed <- impute_missing(x, dat, nk = k)

skip_if_not_installed("poLCA")

test_that("majority voting works", {
  dt <- array(c(2, 3, 2, 2, 2, 3, 1, 1, 2, 3, 2, 2, 2, 3, 1, 1), c(2, 4, 2))
  expect_equal(majority_voting(dt, is.relabelled = TRUE), c(2, 3))
})

test_that("k modes works with or without missing", {
  x <- array(rep(c(rep(1, 10), rep(2, 10), rep(3, 10)), times = 5), c(30, 6, 5))
  xf <- x
  dim(xf) <- c(30, 30)
  set.seed(1)
  kmo.old <- klaR::kmodes(xf, 3)$cluster
  kmo.new <- k_modes(x)
  expect_equal(kmo.old, kmo.new)

  x[3, , 1] <- NA
  kmo.missing <- k_modes(x)
  expect_false(anyNA(kmo.missing))
})

test_that("k modes only clusters if there are multiple assignment vectors", {
  set.seed(1)
  E <- matrix(sample(seq_len(6), 100, replace = TRUE), ncol = 1)
  expect_equal(E, unname(as.matrix(k_modes(E))))
})

test_that("CSPA works", {
  expect_length(CSPA(x, k = 4), nrow(x))
  expect_equal(dplyr::n_distinct(CSPA(x, k = 4)), 4)
})

test_that("Check LCE with hgsc data with 3 consensus_cluster algorithms", {
  y_cts <- LCE(E = x_imputed, k = k, R = 5, sim.mat = "cts")
  y_srs <- LCE(E = x_imputed, k = k, R = 5, sim.mat = "srs")
  y_asrs <- LCE(E = x_imputed, k = k, R = 5, sim.mat = "asrs")
  expect_length(y_cts, 40)
  expect_length(y_srs, 40)
  expect_length(y_asrs, 40)
  expect_type(y_cts, "integer")
  expect_type(y_srs, "integer")
  expect_type(y_asrs, "integer")
})

test_that("Check LCA works", {
  expect_length(LCA(x_imputed), nrow(x))
  expect_equal(dplyr::n_distinct(LCA(x_imputed)), 4)
})
