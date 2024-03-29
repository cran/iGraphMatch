% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lap.R
\name{do_lap}
\alias{do_lap}
\title{Linear (sum) assignment problem}
\usage{
do_lap(score, method = "clue")
}
\arguments{
\item{score}{matrix of pairwise scores}

\item{method}{One of "lapjv", "lapmod", or "clue"}
}
\value{
\code{do_lap} returns a vector which indicates the
 best matching column for each row.
}
\description{
Compute the best bipartite matching
using one of three methods. For an n x n score matrix it find
\eqn{\max_{v\in \Pi_n} \sum_{i=1}^n score_{i, v(i)}}
where \eqn{\Pi_n} denotes all permutations on n objects.
}
\details{
Solves a linear assignment using one of three methods.
 "clue" uses \code{solve_lsap} from the clue package.
 "lapjv" uses the Jonker-Volgenaut approach implemented in this package.
 "lapmod" use a modification of JV that exploits sparsity in the score matrix.
 
 Scores do not need to be non-negative. For "clue" the scores are pre-translated to be
 non-negative which preserves the LAP solution.
}
\examples{
set.seed(12345)
cost <- Matrix::rsparsematrix(10, 10, .5)
cbind(
 do_lap(cost, "lapjv"),
 do_lap(cost, "lapmod"),
 do_lap(cost, "clue")
)

}
\references{
R. Jonker, A. Volgenant (1987). \emph{A shortest augmenting path algorithm
for dense and sparse linear assignment problems}. Computing, pages 325-340.

A. Volgenant (1996). \emph{Linear and Semi-Assignment Problems: A
  Core Oriented Approach}. Computer Ops Res., pages 917-932.

C. H. Papadimitriou and K. Steiglitz (1998). \emph{Combinatorial Optimization:
Algorithms and Complexity}. Courier Corporation.
}
