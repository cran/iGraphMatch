% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gm_convex.R, R/gm_indefinite.R, R/gm_path.R
\name{graph_match_convex}
\alias{graph_match_convex}
\alias{graph_match_indefinite}
\alias{graph_match_PATH}
\title{Frank-Wolfe Graph Matching Methods}
\usage{
graph_match_convex(
  A,
  B,
  seeds = NULL,
  similarity = NULL,
  start = "bari",
  max_iter = 100,
  tol = 1e-05,
  lap_method = NULL
)

graph_match_indefinite(
  A,
  B,
  seeds = NULL,
  similarity = NULL,
  start = "bari",
  max_iter = 20,
  lap_method = NULL
)

graph_match_PATH(
  A,
  B,
  seeds = NULL,
  similarity = NULL,
  epsilon = 1,
  tol = 1e-05,
  max_iter = 20,
  lap_method = NULL
)
}
\arguments{
\item{A}{A matrix, igraph object, or list of either.}

\item{B}{A matrix, igraph object, or list of either.}

\item{seeds}{A vector of integers or logicals, a matrix or a data frame. If
the seed pairs have the same indices in both graphs then seeds can be a
vector. If not, seeds must be  a matrix or a data frame, with the first
column being the indices of \eqn{G_1} and the second column being the
corresponding indices of \eqn{G_2}.}

\item{similarity}{A matrix. An \code{n-by-n} matrix containing vertex similarities.}

\item{start}{A matrix or a character. Any \code{nns-by-nns} matrix or
character value like "bari", "rds" or "convex" to initialize the starting matrix.}

\item{max_iter}{A number. Maximum number of replacing matches.}

\item{tol}{A number. Tolerance of edge disagreements.}

\item{lap_method}{Choice for lap method. One of "lapjv", "lapmod", or "clue".}

\item{epsilon}{A small number}
}
\value{
\code{graph_match_indefinite}, \code{graph_match_convex} and \code{graph_match_PATH}
  return an object of class "\code{\link{graphMatch}}" which is a list containing the following
  components:

  \describe{
    \item{corr_A}{matching correspondence in \eqn{G_1}}
    \item{corr_B}{matching correspondence in \eqn{G_2}}
    \item{soft}{the doubly stochastic matrix from the last iteration with which one can
          extract more than one matching candidates}
    \item{iter}{number of iterations until convergence or reaches the \code{max_iter}}
    \item{max_iter}{Maximum number of replacing matches}
    \item{lap_method}{Choice for solving the LAP}
    \item{seeds}{a vector of logicals indicating if the corresponding vertex is a seed}
  }
}
\description{
Match two given graphs, returns a list of graph matching
  results, including matching correspondence vector of \eqn{G_2} with respect
  to \eqn{G_1}, doubly stochastic matrix and permutation matrix.
}
\examples{
\donttest{
cgnp_pair <- sample_correlated_gnp_pair(n = 10, corr =  0.9, p =  0.5)
g1 <- cgnp_pair$graph1
g2 <- cgnp_pair$graph2
# match G_1 & G_2 with no seeds
gm(g1, g2, method = "convex", max_iter = 10)
seeds <- 1:10 <= 3
gm(g1, g2, seeds, method = "convex", max_iter = 10)
}



# match G_1 & G_2 with some known node pairs as seeds
cgnp_pair <- sample_correlated_gnp_pair(n = 10, corr =  0.3, p =  0.5)
g1 <- cgnp_pair$graph1
g2 <- cgnp_pair$graph2
seeds <- 1:10 <= 3
GM_bari <- gm(g1, g2, seeds, method = "indefinite", start = "bari")
GM_bari
GM_bari[!GM_bari$seeds] # matching correspondence for non-seeds

summary(GM_bari, g1, g2, true_label = 1:10)

# match G_1 & G_2 with some incorrect seeds
hard_seeds <- matrix(c(4,6,5,4),2)
seeds <- rbind(as.matrix(check_seeds(seeds, nv = 10)$seeds),hard_seeds)
GM_badseed <- gm(g1, g2, seeds, method = "indefinite")

GM_badseed[] # get the corresponding permutation matrix
GM_badseed \%*\% g2 # permute the second graph according to match result: PBP^T
GM_badseed$soft # doubly stochastic matrix from the last step of Frank-Wolfe iterations
GM_badseed$iter # number of iterations
GM_badseed$max_iter # preset maximum number of iterations: 20

# match two multi-layer graphs
gp_list <- replicate(3, sample_correlated_gnp_pair(20, .3, .5), simplify = FALSE)
A <- lapply(gp_list, function(gp)gp[[1]])
B <- lapply(gp_list, function(gp)gp[[2]])

match_multi_layer <- gm(A, B, seeds = 1:10, method = "indefinite", start = "bari", max_iter = 20)
summary(match_multi_layer, A, B)

# match G_1 & G_2 using PATH algorithm
gm(g1, g2, method = "PATH")


}
\references{
Y. Aflalo and A. Bronstein and R. Kimmel (2014), \emph{On convex
relaxation of graph isomorphism}. Proceedings of the National Academy of Sciences,
pages 2942-2947.

V. Lyzinski and D. E. Fishkind and M. Fiori and J. T. Vogelstein and C. E. Priebe
and G. Sapiro (2016), \emph{Graph Matching: Relax at Your Own Risk}. IEEE TPAMI, pages 60-73.

V. Lyzinski and D. E. Fishkind and C. E. Priebe (2014), \emph{Seeded Graph Matching
for Correlated Erdos-Renyi Graphs}.J. Mach. Learn. Res., pages 3513-3540.

M. Zaslavskiy, F. Bach and J. Vert (2009), \emph{A Path following
algorithm for the graph matching problem}. IEEE Trans Pattern Anal Mach Intell,
pages 2227-2242.
}
\keyword{internal}
