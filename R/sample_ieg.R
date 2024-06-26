#' @title Sample graphs from edge probability matrix and correlation matrix
#'
#' @description Sample a pair of graphs with specified edge probability and
#'   correlation between each pair of vertices.
#'
#' @param n An integer. Number of total vertices for the sampled graphs.
#' @param p_mat An \code{n-by-n} matrix. Edge probability matrix, each entry
#'   should be in the open (0,1) interval.
#' @param c_mat An \code{n-by-n} matrix. The target Pearson correlation matrix,
#'   each entry should be in the open (0,1) interval.
#' @param directed A logical. \code{TRUE} if the sampled graphs are directed.
#' @param X A matrix. Dot products matrix, each entry must be in open (0,1)
#'   interval.
#' @param corr A number. The target Pearson correlation between the adjacency
#'   matrices of the generated graphs. It must be in open (0,1) interval.
#' @param directed Logical scalar, whether to generate directed graphs.
#' @param loops Logical scalar, whether self-loops are allowed in the graph.
#' @param permutation A numeric vector,permute second graph.
#' @param ncore An integer. Number of core vertices.
#' @param ... Passed to \code{sample_correlated_ieg_pair}.
#'
#' @rdname sample_ieg
#' @return \code{sample_correlated_ieg_pair} returns two igraph objects named
#'   \code{graph1} and \code{graph2}. If sample two graphs with junk vertices,
#'   the first \code{ncore} vertices are core vertices and the rest are junk
#'   vertices.
#'
#' @examples
#' n <- 50
#' p_mat <- matrix(runif(n^2),n)
#' c_mat <- matrix(runif(n^2),n)
#' sample_correlated_ieg_pair(n,p_mat,c_mat,ncore=40)
#'
#' @references S. Young and E. Scheinerman (2007), \emph{Random Dot Product
#'   Graph Models for Social Networks}. Proceedings of the 5th International
#'   Conference on Algorithms and Models for the Web-graph, pages 138-149.
#' @references F. Fang and D. Sussman and V. Lyzinski (2018), \emph{Tractable
#'   Graph Matching via Soft Seeding}. \url{https://arxiv.org/abs/1807.09299}.
#'
#'
#' @seealso \code{\link{sample_correlated_gnp_pair}},
#'   \code{\link{sample_correlated_sbm_pair}}
#'
#' @export
sample_correlated_ieg_pair<- function(n, p_mat, c_mat, ncore=n, directed=FALSE, loops=FALSE, permutation=1:n){
  if (any(p_mat < 0 | p_mat > 1 | is.na(p_mat))) {
    stop("p_mat must have all entries between 0 and 1 and non-NA.")
  }
  if (any(c_mat < 0 | c_mat > 1 | is.na(c_mat))) {
    stop("c_mat must have all entries between 0 and 1 and non-NA.")
  }

  if(ncore > n){
    stop("ncore must be at most n.")
  }
  if(ncore != n){
    c_mat[((ncore+1):n),((ncore+1):n)] <- 0
  }

  if(nrow(p_mat) != n | ncol(p_mat) != n | nrow(c_mat) != n | ncol(c_mat) != n){
    stop("Edge probability matrix and Pearson correlation matrix should be square matrices of size n.")
  }
  g1 <- matrix(stats::rbinom(n^2,1,p_mat),n)
  z0 <- matrix(stats::rbinom(n^2,1,p_mat*(1-c_mat)),n)
  z1 <- matrix(stats::rbinom(n^2,1,p_mat*(1-c_mat)+c_mat),n)
  g2 <- z1*g1+z0*(1-g1)

  if(directed){
    mode <- "directed"
  } else{
    g1[row(g1)>=col(g1)] <- 0
    g1 <- g1 + t(g1)
    g2[row(g2)>=col(g2)] <- 0
    g2 <- g2 + t(g2)
    mode <- "undirected"
  }
  list(graph1 = graph_from_adjacency_matrix(g1,
         mode = mode, diag = loops),
       graph2 = igraph::permute(graph_from_adjacency_matrix(g2,
         mode = mode, diag = loops),permutation))
}

#' @rdname sample_ieg
#' @return \code{sample_correlated_rdpg_pair} returns two igraph objects named
#' \code{graph1} and \code{graph2} that are sampled from random dot product
#' graphs model. If sample two graphs with junk vertices, the first
#' \code{ncore} vertices are core vertices and the rest are junk vertices.
#'
#' @examples
#' ## sample a pair of igraph objects from random dot
#' ## product graphs model with dimension 3 and scale 8
#' n <- 50
#' xdim <- 3
#' scale <- 8
#' X <- matrix(rgamma(n*(xdim+1),scale,1),n,xdim+1)
#' X <- X/rowSums(X)
#' X <- X[,1:xdim]
#' sample_correlated_rdpg_pair(X,corr=0.5,ncore=40)
#'
#' @export
sample_correlated_rdpg_pair <- function(X, corr, ncore=nrow(X), ...){
  p_mat <- X %*% t(X)
  n <- nrow(X)
  if(length(corr) == 1){
    c_mat <- matrix(0, n, n)
    c_mat[1:ncore, 1:ncore] <- corr
  }else{
    c_mat <- corr
  }
  sample_correlated_ieg_pair(n, p_mat, c_mat, ncore, ...)
}
