#' @title Sample graphs pair from stochastic block model
#'
#' @description Sample a pair of random graphs from stochastic block model with
#'   correlation between two graphs being \code{corr} and edge probability being
#'   \code{p}.
#'
#' @param n An integer. Number of vertices in the graph.
#' @param pref.matrix The matrix giving the Bernoulli rates. This is a
#'   \code{K-by-K} matrix, where \code{k} is the number of groups. The
#'   probability of creating an edge between vertices from groups \code{i} and
#'   \code{j} is given by element \code{i,j}. For undirected graphs, this matrix
#'   must be symmetric.
#' @param block.sizes A numeric vector. Give the number of vertices in each
#'   group. The sum of the vector must match the number of vertices.
#' @param corr A number. The target Pearson correlation between the adjacency
#'   matrices of the generated graphs. It must be in open (0,1) interval.
#' @param permutation A numeric vector, permute second graph.
#' @param core.block.sizes A numeric vector. Give the number of core vertices in
#'   each group. Entries should be smaller than \code{block.sizes} and the
#'   vector length should be the same as \code{block.sizes}.
#' @param ... Passed to \code{sample_sbm}.
#'
#' @rdname sample_sbm
#' @return Returns a list of two igraph object, named \code{graph1} and
#'   \code{graph2}. If sample two graphs with junk vertices, in each
#'   corresponding block the first \code{core.block.sizes} vertices are core
#'   vertices and the rest are junk vertices.
#'
#' @references P. Holland and K. Laskey and S. Leinhardt (1983),
#'   \emph{Stochastic Blockmodels: First Steps}. Social Networks, pages 109-137.
#' @references F. Fang and D. Sussman and V. Lyzinski (2018), \emph{Tractable
#'   Graph Matching via Soft Seeding}. \url{https://arxiv.org/abs/1807.09299}.
#'
#'
#' @examples
#' pm <- cbind( c(.1, .001), c(.001, .05) )
#' sample_correlated_sbm_pair(n=1000, pref.matrix=pm, block.sizes=c(300,700), corr=0.5)
#' sample_correlated_sbm_pair(n=1000, pref.matrix=pm, block.sizes=c(300,700), corr=0.5,
#' core.block.sizes=c(200,500))
#'
#' @seealso \code{\link{sample_correlated_gnp_pair}},
#'   \code{\link{sample_correlated_rdpg_pair}}

#'
#' @export
#'

sample_correlated_sbm_pair <- function(
  n, pref.matrix, block.sizes, corr, core.block.sizes=NULL, permutation=1:n, ...){
  if (any(pref.matrix < 0 | pref.matrix > 1 | is.na(pref.matrix))) {
    stop("pref.matrix must have all entries between 0 and 1 and non-NA.")
  }
  if (any(corr < 0 | corr > 1 | is.na(corr))) {
    stop("corr must have all entries between 0 and 1 and non-NA.")
  }


  if(is.null(core.block.sizes)){
    sample_correlated_sbm_pair_no_junk(n, pref.matrix, block.sizes, corr, permutation, ...)
  } else if(sum(block.sizes >= core.block.sizes) == length(block.sizes)){
    sample_correlated_sbm_pair_w_junk(n, pref.matrix, block.sizes, corr, core.block.sizes, permutation, ...)
  } else{
    stop("Number of core vertices must be at most block size in each block.")
  }
}

sample_correlated_sbm_pair_no_junk <- function(n, pref.matrix, block.sizes, corr, permutation=1:n, ...){

  K <- length(block.sizes)
  # Make the first graph
  graph1 <- sample_sbm(n,pref.matrix,block.sizes,...)

  # Make two graphs which will be used to make the
  # second graph
  corr.matrix <- (1-corr)*pref.matrix
  Z0 <- sample_sbm(n,corr.matrix,block.sizes,...)
  Z1 <- sample_sbm(n,corr.matrix+corr,block.sizes,...)

  graph2 <- Z1 %s% graph1 %u% (Z0-graph1)

  list(graph1=graph1,graph2=igraph::permute(graph2,permutation))
}


sample_correlated_sbm_pair_w_junk <- function(
  n, pref.matrix, block.sizes, corr, core.block.sizes, permutation=1:n, ...){

  K <- length(block.sizes)
  ncore <- sum(core.block.sizes)
  core <- 1:ncore
  junk <- (ncore+1):n

  junk.block.sizes <- block.sizes - core.block.sizes
  all.block.sizes <- c(core.block.sizes,junk.block.sizes)
  all.pref.matrix <- kronecker(matrix(1,2,2),pref.matrix)
  # Make the first graph
  graph1 <- sample_sbm(n,all.pref.matrix,all.block.sizes,...)

  # Make two graphs which will be used to make the
  # second graph
  all.pref.matrix <- kronecker(matrix(c(1-corr,1,1,1),2),pref.matrix)
  Z0 <- sample_sbm(n,all.pref.matrix,all.block.sizes,...)

  all.pref.matrix[1:K,1:K] <- all.pref.matrix[1:K,1:K] + corr
  Z1 <- sample_sbm(n,all.pref.matrix,all.block.sizes,...)

  graph2 <- Z1 %s% graph1 %u% (Z0-graph1)

  list(graph1=graph1,graph2=igraph::permute(graph2,permutation))
}

