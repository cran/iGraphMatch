


# #' @title Matching performance summary
# #'
# #' @description Get a summary of the matching result and measures of the matching performance
# #' based on several evaluation metrics associated with nodes and edges of two graphs.
# #'
# #' @param match Graph matching result see \link[=graph_match_FW]{graph match methods}.
# #' @param A A matrix or an igraph object. Adjacency matrix of \eqn{G_1}.
# #' @param B A matrix or an igraph object. Adjacency matrix of \eqn{G_2}.
# #' @param true_label A vector. NULL if the true correspondence
# #'  between two graphs is unknown. A vector indicating the
# #'  true correspondence in the second graph if the true
# #'  correspondence is known.
# #' @param corr Correspondence data frame as given by match$corr
# #' @param directed Whether the graphs should be treated as directed
# #'  or undirected. NULL defaults to !isSymmetric(A).
# #'
# #' @rdname match_report
# #'
# #' @return \code{match_report} returns the match object
# #'  evaluation metrics including number of matches, true matches,
# #'  and a data frame with edge correctness information.
# #'  \code{edge_match_info} returns this data frame with columns
# #'  for number of common edges, missing edges, extra edges, and
# #'  common non-edges, and Frobenius norm.
# #'
# #'
# #' @details For multilayered graphs information is given per layer.
# #'  For weighted graphs the counts are based on non-zero entries.
# #'  Equality of weights is not tested.
# #'  If you want to ignore seeds in the edge match info you must
# #'  remove them from corr/match$corr.
# #'
# #' @examples
# #' graphs <- sample_correlated_gnp_pair(10, .5, .3)
# #' A <- graphs$graph1
# #' B <- graphs$graph2
# #' res <- graph_match_percolation(A, B, 1:4)
# #' match_report(res, A, B)
# #'
# #' gp_list <- replicate(3,
# #'    sample_correlated_gnp_pair(100, .8, .3),
# #'    simplify = FALSE)
# #' A <- lapply(gp_list, function(gp)gp[[1]])
# #' B <- lapply(gp_list, function(gp)gp[[2]])
# #' corr <- data.frame(corr_A = 1:100, corr_B = 1:100)
# #' edge_match_info(corr, A, B)
# #' @export
# #'
# match_report <- function(match, A, B, true_label = NULL, directed = NULL){
#   graph_pair <- check_graph(A, B)
#   A <- graph_pair[[1]]
#   B <- graph_pair[[2]]


#   # if(max(A)>1 || max(B)>1){
#   #   warning("Common egdes have positive weights but not necessarily have same weights.")
#   # }


#   cat("Call: \n")
#   print(match$call)

#   # Matched nodes
#   corr <- match$corr
#   match$n.match <- nrow(corr) - nrow(match$seeds)
#   cat("\n# Matches:", match$n.match)
#   if(!is.null(true_label)){
#     match$n.true.match <-
#       sum(true_label[corr$corr_A] == corr$corr_B) - nrow(match$seeds)
#     cat("\n# True Matches: ", match$n.true.match)
#   }
#   cat("\n")

#   # Matched edges
#   match$edge_match_info <-
#     edge_match_info(corr, A, B, directed)

#   ep <- as.data.frame(t(match$edge_match_info))
#   colnames(ep) <- NULL
#   print(ep)

#   # objective value: ||A-PBP^T||_F
#   cat("\n")
#   invisible(match)
# }


# #' @title document
# #' @description Return aligned versions of A and B according to
# #'  a result of match method
# #'
# #' @param match Result from a a graph matching method.
# #' @param A A matrix, igraph object, or list of either.
# #'  Likely used in the call for creating match.
# #' @param B A matrix, igraph object, or list of either.
# #'  Likely used in the call for creating match.
# #'
# #' @return A list of aligned graphs named \code{A_m} and \code{B_m}.
# #'
# #' @export
# matched_adjs <- function(match, A, B){
#   graph_pair <- check_graph(A, B)
#   A <- graph_pair[[1]]
#   B <- graph_pair[[2]]

#   list(A_m = A[match$corr$corr_A, match$corr$corr_A],
#     B_m = B[match$corr$corr_B, match$corr$corr_B])
# }



# #' @rdname match_report
# #'
# #' @section TODO: support weighted? loops? ...?
# #'
# #' @export

edge_match_info <- function(corr, A, B,
    directed = NULL) {
  graph_pair <- check_graph(A, B)
  A <- graph_pair$g1
  B <- graph_pair$g2
  nv <- min(graph_pair$totv1, graph_pair$totv2,
    nrow(corr))


  if (is.null(directed)) {
    directed <- !(Matrix::isSymmetric(A[[1]])) & !(Matrix::isSymmetric(B[[1]]))
  }

  nA <- nrow(A[[1]])
  nB <- nrow(B[[1]])
  # implement for non-equal

  # if(max(A)>1 || max(B)>1){
  # }
  layers <- names(A)
  if (is.null(layers)) {
    layers <- seq_along(A)
  }
  l <- length(A)
  corr_A <- corr$corr_A[seq(nv)]
  corr_B <- corr$corr_B[seq(nv)]


  Reduce(rbind, lapply(seq_along(layers), function(i) {
    A_m <- A[[i]][corr_A, corr_A]
    B_m <- B[[i]][corr_B, corr_B]
    # computations below are fast for sparse matrices

    res <- list()
    nzA <- which(A_m > 0)
    nzB <- which(B_m > 0)

    # non-zeros in both A and B
    res$common_edges <- length(intersect(nzA, nzB))
    # non-zero in A but not B
    res$missing_edges <- length(setdiff(nzA, nzB))
    # zero in A but not B
    res$extra_edges <- length(setdiff(nzB, nzA))
    # zero in A and B, ie the rest
    # res$common_non_edges <- nA^2 - nA - Reduce(sum, res)

    # if undirected, divide by 2
    if (!directed) {
      res <- lapply(res, function(x) x / 2)
    }
    res$fnorm <- Matrix::norm(A_m - B_m, "F")
    if(l > 1) {
      res <- c(layer = layers[i], res)
    }
    as.data.frame(res)
  }))

}

match_plot_igraph <- function(A, B, match,
  color = TRUE, linetype = TRUE, ...) {

  ch <- check_graph(A, B, same_order = FALSE, as_igraph = TRUE)

  nv <- min(ch$totv1, ch$totv2, nrow(match@corr))


  corr_A <- match@corr$corr_A[seq(nv)]
  corr_B <- match@corr$corr_B[seq(nv)]

  if(is.null(igraph::V(A)$name)){
    A <- igraph::set_vertex_attr(A, "name", corr_A, corr_A)
  }

  A <- igraph::permute(
    igraph::induced_subgraph(A, corr_A),
    rank(corr_A)
  )
  B <- igraph::permute(
    igraph::induced_subgraph(B, corr_B),
    rank(corr_B)
  )



  igraph::E(A)$in_A <- "A"
  igraph::E(B)$in_B <- "B"

  g <- igraph::union(A, B, byname = FALSE)
  igraph::E(g)$edge_match <-
    ifelse(is.na(igraph::E(g)$in_B), "Only A",
      ifelse(is.na(igraph::E(g)$in_A), "Only B", "Both"))
  g <- igraph::delete_edge_attr(g, "in_A")
  g <- igraph::delete_edge_attr(g, "in_B")

  pal <- c("#888888", "#4444AA", "#AA4444")
  igraph::E(g)$color <- pal[1]
  igraph::E(g)$lty <- 1
  if (color) {
    igraph::E(g)$color <-
      pal[as.numeric(factor(igraph::E(g)$edge_match,
        levels = c("Both", "Only A", "Only B")))]
  }
  if (linetype) {
    igraph::E(g)$lty <-
      as.numeric(factor(igraph::E(g)$edge_match,
        levels = c("Both", "Only A", "Only B")))
  }

  graphics::plot(g, ...)
  invisible(g)
}


match_plot_matrix <- function(A, B, match, col.regions = NULL, at = NULL, colorkey = NULL, ...) {
  ch <- check_graph(A, B, same_order = FALSE, as_list = FALSE)
  nv <- min(ch$totv1, ch$totv2, nrow(match@corr))

  corr_A <- match@corr$corr_A[seq(nv)]
  corr_B <- match@corr$corr_B[seq(nv)]

  A <- ch$g1[corr_A, corr_A]
  B <- ch$g2[corr_B, corr_B]

  m <- A - B
  m_max <- max(abs(m))
  if (is.null(at)) {
    at <- seq(-m_max * 1.0001, m_max * 1.0001, length.out = 16)
  }
  if (is.null(col.regions)) {
    col <- grDevices::colorRampPalette(
      c("#AA4444", "#888888", "#4444AA"))
    col.regions <- col(length(at) - 1)
  }
  if (is.null(colorkey)) {
    colorkey <- list(at = at)
  }
  print(image(m, col.regions = col.regions, at = at, colorkey = colorkey, ...))
  invisible(m)
}
