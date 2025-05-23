#' Graphical Displays
#'
#' Graph cumulative distribution function (CDF) graphs, relative change in area
#' under CDF curves, heatmaps, and cluster assignment tracking plots.
#'
#' `graph_cdf` plots the CDF for consensus matrices from different algorithms.
#' `graph_delta_area` calculates the relative change in area under CDF curve
#' between algorithms. `graph_heatmap` generates consensus matrix heatmaps for
#' each algorithm in `x`. `graph_tracking` tracks how cluster assignments change
#' between algorithms. `graph_all` is a wrapper that runs all graphing
#' functions.
#'
#' @param x an object from [consensus_cluster()]
#' @param mat same as `x`, or a list of consensus matrices computed from `x` for
#'   faster results
#' @param cl same as `x`, or a matrix of consensus classes computed from `x` for
#'   faster results
#' @return Various plots from \code{graph_*{}} functions. All plots are
#'   generated using `ggplot`, except for `graph_heatmap`, which uses
#'   [NMF::aheatmap()]. Colours used in `graph_heatmap` and `graph_tracking`
#'   utilize [RColorBrewer::brewer.pal()] palettes.
#' @name graphs
#' @author Derek Chiu
#' @export
#' @examples
#' # Consensus clustering for 3 algorithms
#' library(ggplot2)
#' set.seed(911)
#' x <- matrix(rnorm(80), ncol = 10)
#' CC1 <- consensus_cluster(x, nk = 2:4, reps = 3,
#' algorithms = c("hc", "pam", "km"), progress = FALSE)
#'
#' # Plot CDF
#' p <- graph_cdf(CC1)
#'
#' # Change y label and add colours
#' p + labs(y = "Probability") + stat_ecdf(aes(colour = k)) +
#' scale_color_brewer(palette = "Set2")
#'
#' # Delta Area
#' p <- graph_delta_area(CC1)
#'
#' # Heatmaps with column side colours corresponding to clusters
#' CC2 <- consensus_cluster(x, nk = 3, reps = 3, algorithms = "hc", progress =
#' FALSE)
#' graph_heatmap(CC2)
#'
#' # Track how cluster assignments change between algorithms
#' p <- graph_tracking(CC1)
graph_cdf <- function(mat) {
  dat <- get_cdf(mat)
  p <- ggplot(dat, aes(x = !!sym("CDF"), colour = !!sym("k"))) +
    stat_ecdf() +
    facet_wrap(~Method) +
    labs(x = "Consensus Index",
         y = "CDF",
         title = "Consensus Cumulative Distribution Functions")
  print(p)
  return(p)
}

#' @rdname graphs
#' @references https://stackoverflow.com/questions/4954507/calculate-the-area-under-a-curve
#' @export
graph_delta_area <- function(mat) {
  dat <- get_cdf(mat) %>%
    dplyr::group_by(.data$Method, .data$k) %>%
    dplyr::summarize(AUC = sum(diff(seq(0, 1, length.out = length(.data$k))) *
                                 (utils::head(.data$CDF, -1) + utils::tail(.data$CDF, -1))) / 2) %>%
    dplyr::mutate(da = c(.data$AUC[1], diff(.data$AUC) / .data$AUC[-length(.data$AUC)]))
  if (length(unique(dat$k)) > 1) {
    p <- ggplot(dat, aes(x = !!sym("k"), y = !!sym("da"))) +
      geom_line(group = 1) +
      geom_point() +
      facet_wrap(~Method) +
      labs(y = "Relative change in Area under CDF curve",
           title = "Delta Area")
    print(p)
    return(p)
  }
}

#' Calculate CDF for each clustering algorithm at each k
#' @noRd
get_cdf <- function(mat) {
  if (inherits(mat, "array")) {
    mat <- consensus_combine(mat, element = "matrix")
  }
  mat %>%
    purrr::modify_depth(2, ~ .x[lower.tri(.x, diag = TRUE)]) %>%
    purrr::imap(~ purrr::set_names(.x, paste(.y, names(.x), sep = "."))) %>%
    purrr::flatten_dfc() %>%
    tidyr::gather("Group", "CDF", names(.)) %>%
    tidyr::separate("Group", c("k", "Method"), sep = "\\.") %>%
    dplyr::mutate(k = factor(.data$k, levels = as.integer(unique(.data$k))))
}

#' @param main heatmap title. If `NULL` (default), the titles will be taken from
#'   names in `mat`
#'
#' @rdname graphs
#' @export
graph_heatmap <- function(mat, main = NULL) {
  sample_names <- rownames(mat)
  if (inherits(mat, "array")) {
    mat <- consensus_combine(mat, element = "matrix")
  }
  dat <- mat %>%
    purrr::list_flatten(name_spec = "{inner} k={outer}") |>
    purrr::map(~ {
      .x |>
        magrittr::set_colnames(sample_names) |>
        magrittr::set_rownames(sample_names)
    })
  main <- paste(main %||% names(dat), "Consensus Matrix")
  assertthat::assert_that(length(main) == length(purrr::flatten(mat)))
  # Number of clusters of each run
  clusters <- mat |>
    purrr::map_int(length) |>
    purrr::imap(~ rep(.y, each = .x)) |>
    purrr::flatten_chr() |>
    as.numeric()
  # Annotate samples with cluster index
  annotation_col <- purrr::map2(dat, clusters, ~ {
    data.frame(Cluster = paste0("C", hc(stats::dist(.x), k = .y))) |>
      magrittr::set_rownames(sample_names)
  })
  # Annotation colour palette: RColorBrewer Set2
  annotation_pal <- c(
    "#66C2A5",
    "#FC8D62",
    "#8DA0CB",
    "#E78AC3",
    "#A6D854",
    "#FFD92F",
    "#E5C494",
    "#B3B3B3"
  )
  # Map annotations to colour palette
  annotation_colors <- purrr::map2(annotation_col, clusters, ~ {
    list(Cluster = stats::setNames(head(annotation_pal, .y), unique(unlist(.x))))
  })
  # Heatmap colour palette: RColorBrewer PuBuGn
  heatmap_pal <- c(
    "#FFF7FB",
    "#ECE2F0",
    "#D0D1E6",
    "#A6BDDB",
    "#67A9CF",
    "#3690C0",
    "#02818A",
    "#016C59",
    "#014636"
  )
  # Create all heatmaps
  purrr::pwalk(list(dat, annotation_col, annotation_colors, main), ~ {
    pheatmap::pheatmap(
      mat = ..1,
      color = grDevices::colorRampPalette(heatmap_pal)(10),
      clustering_callback = callback,
      clustering_method = "average",
      treeheight_row = 0,
      border_color = NA,
      show_rownames = FALSE,
      show_colnames = FALSE,
      annotation_col = ..2,
      annotation_colors = ..3,
      annotation_names_col = FALSE,
      main = ..4
    )
  })
}

#' @rdname graphs
#' @export
graph_tracking <- function(cl) {
  if (inherits(cl, "array")) {
    cl <- consensus_combine(cl, element = "class")
  }
  dat <- cl %>%
    purrr::imap(~ `colnames<-`(.x, paste(.y, colnames(.x), sep = "."))) %>%
    do.call(cbind, .) %>%
    as.data.frame() %>%
    tidyr::gather("Group", "Class", names(.)) %>%
    tidyr::separate("Group", c("k", "Method"), sep = "\\.") %>%
    cbind(Samples = seq_len(unique(purrr::map_int(cl, nrow)))) %>%
    dplyr::mutate_at(dplyr::vars(c("Class", "Method", "Samples")), factor)
  if (length(unique(dat$k)) > 1) {
    p <- ggplot(dat, aes(x = !!sym("Samples"), y = !!sym("k"))) +
      geom_tile(aes(fill = !!sym("Class"))) +
      facet_wrap(~Method) +
      scale_fill_brewer(palette = "Set2") +
      ggtitle("Tracking Cluster Assignments Across k") +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank())
    print(p)
    return(p)
  }
}

#' @rdname graphs
#' @export
graph_all <- function(x) {
  mat <- consensus_combine(x, element = "matrix")
  cl <- consensus_combine(x, element = "class")
  graph_cdf(mat)
  graph_delta_area(mat)
  graph_heatmap(x)
  graph_tracking(cl)
}

#' Comparing ranked Algorithms vs internal indices (ii) in heatmap
#' @inheritParams dice
#' @param E object in `dice`
#' @param clusters object in `dice`
#' @noRd
algii_heatmap <- function(data, nk, E, clusters, ref.cl = NULL) {
  # Cluster list to keep
  cl.list <- E %>%
    consensus_combine(element = "class") %>%
    magrittr::extract(as.character(nk))

  # Final cluster object construction depends on value of nk
  if (length(nk) > 1) {
    fc <- purrr::map2(cl.list, nk,
                      ~ magrittr::set_colnames(.x, paste_k(colnames(.), .y))) %>%
      purrr::map2(split_clusters(clusters), cbind) %>%
      purrr::map(~ apply(., 2, relabel_class, ref.cl = ref.cl %||% .[, 1])) %>%
      do.call(cbind, .) %>%
      as.data.frame()
  } else {
    fc <- cl.list %>%
      do.call(cbind, .) %>%
      cbind.data.frame(clusters) %>%
      purrr::map_df(relabel_class, ref.cl = ref.cl %||% .[, 1])
  }

  # Internal indices
  ii <- ivi_table(fc, data)

  # Heatmap: order algorithms by ranked ii, remove indices with NaN
  hm <- ii %>%
    dplyr::select(-"Algorithms") %>%
    magrittr::extract(match(consensus_rank(ii, 1)$top.list, rownames(.)),
                      purrr::map_lgl(., ~ all(!is.nan(.x))))

  # Plot heatmap with annotated colours, column scaling, no further reordering
  NMF::aheatmap(
    hm,
    annCol = data.frame(Criteria = c(
      rep("Maximized", 5), rep("Minimized", ncol(hm) - 5)
    )),
    annColors = list(Criteria = stats::setNames(
      c("darkgreen", "deeppink4"), c("Maximized", "Minimized")
    )),
    Colv = NA,
    Rowv = NA,
    scale = "column",
    col = "PiYG",
    main = "Ranked Algorithms on Internal Validity Indices"
  )
  # Heatmap palette: RColorBrewer PiYg
  heatmap_pal <- c(
    "#8E0152",
    "#C51B7D",
    "#DE77AE",
    "#F1B6DA",
    "#FDE0EF",
    "#F7F7F7",
    "#E6F5D0",
    "#B8E186",
    "#7FBC41",
    "#4D9221",
    "#276419"
  )
  pheatmap::pheatmap(
    mat = hm,
    color = grDevices::colorRampPalette(heatmap_pal)(12),
    border_color = NA,
    scale = "column",
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    annotation_col = data.frame(Criteria = c(
      rep("Maximized", 5), rep("Minimized", ncol(hm) - 5)
    )) |>
      magrittr::set_rownames(colnames(hm)),
    annotation_colors = list(Criteria = stats::setNames(
      c("darkgreen", "deeppink4"), c("Maximized", "Minimized")
    )),
    annotation_names_col = FALSE,
    main = "Ranked Algorithms on Internal Validity Indices"
  )
}

#' Split clusters matrix into list based on value of k
#' @noRd
split_clusters <- function(clusters) {
  tc <- t(clusters)
  split.data.frame(x = tc,
                   f = stringr::str_split_fixed(
                     string = rownames(tc),
                     pattern = " ",
                     n = 2
                   )[, 2]) %>%
    purrr::map(t)
}
