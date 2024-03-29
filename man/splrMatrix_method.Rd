% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/splr_matrix_class.R
\name{show,splrMatrix-method}
\alias{show,splrMatrix-method}
\alias{splrMatrix_method}
\alias{print,splrMatrix-method}
\alias{\%*\%,splrMatrix,splrMatrix-method}
\alias{\%*\%,splrMatrix,matrix_list-method}
\alias{\%*\%,matrix_list,splrMatrix-method}
\alias{\%*\%,Matrix,splrMatrix-method}
\alias{\%*\%,matrix,splrMatrix-method}
\alias{\%*\%,numeric,splrMatrix-method}
\alias{\%*\%,numLike,splrMatrix-method}
\alias{\%*\%,ANY,splrMatrix-method}
\alias{dim,splrMatrix-method}
\alias{length,splrMatrix-method}
\alias{\%*\%,splrMatrix,Matrix-method}
\alias{\%*\%,splrMatrix,matrix-method}
\alias{\%*\%,splrMatrix,numeric-method}
\alias{\%*\%,splrMatrix,numLike-method}
\alias{\%*\%,splrMatrix,ANY-method}
\alias{*,splrMatrix,splrMatrix-method}
\alias{*,Matrix,splrMatrix-method}
\alias{*,splrMatrix,ddiMatrix-method}
\alias{*,ddiMatrix,splrMatrix-method}
\alias{*,matrix,splrMatrix-method}
\alias{*,numeric,splrMatrix-method}
\alias{*,ANY,splrMatrix-method}
\alias{*,splrMatrix,matrix-method}
\alias{*,splrMatrix,Matrix-method}
\alias{*,splrMatrix,numeric-method}
\alias{*,splrMatrix,ANY-method}
\alias{/,splrMatrix,matrix-method}
\alias{/,splrMatrix,Matrix-method}
\alias{/,splrMatrix,ANY-method}
\alias{+,splrMatrix,splrMatrix-method}
\alias{-,splrMatrix,splrMatrix-method}
\alias{+,splrMatrix,Matrix-method}
\alias{+,splrMatrix,numeric-method}
\alias{+,splrMatrix,ANY-method}
\alias{-,splrMatrix,missing-method}
\alias{-,splrMatrix,Matrix-method}
\alias{-,splrMatrix,ddiMatrix-method}
\alias{-,splrMatrix,numeric-method}
\alias{-,splrMatrix,ANY-method}
\alias{+,Matrix,splrMatrix-method}
\alias{+,numeric,splrMatrix-method}
\alias{+,ANY,splrMatrix-method}
\alias{-,Matrix,splrMatrix-method}
\alias{-,numeric,splrMatrix-method}
\alias{-,ANY,splrMatrix-method}
\alias{norm,splrMatrix,character-method}
\alias{rowSums,splrMatrix-method}
\alias{colSums,splrMatrix-method}
\alias{rowMeans,splrMatrix-method}
\alias{colMeans,splrMatrix-method}
\alias{sum,splrMatrix,ANY-method}
\alias{mean,splrMatrix-method}
\alias{[,splrMatrix,missing,missing,missing-method}
\alias{[,splrMatrix,numeric,numeric,logical-method}
\alias{[,splrMatrix,missing,numeric,logical-method}
\alias{[,splrMatrix,missing,numeric,missing-method}
\alias{[,splrMatrix,missing,logical,logical-method}
\alias{[,splrMatrix,missing,logical,missing-method}
\alias{[,splrMatrix,numeric,missing,logical-method}
\alias{[,splrMatrix,numeric,missing,missing-method}
\alias{[,splrMatrix,logical,missing,logical-method}
\alias{[,splrMatrix,logical,missing,missing-method}
\alias{[,splrMatrix,numeric,ANY,logical-method}
\alias{[,splrMatrix,numeric,logical,logical-method}
\alias{[,splrMatrix,numeric,ANY,missing-method}
\alias{[,splrMatrix,logical,ANY,ANY-method}
\alias{[,splrMatrix,logical,ANY,missing-method}
\alias{[,splrMatrix,logical,numeric,missing-method}
\alias{[,splrMatrix,logical,numeric,logical-method}
\alias{[,splrMatrix,matrix,missing,missing-method}
\alias{[<-,splrMatrix,numeric,numeric,ANY-method}
\alias{[<-,splrMatrix,numeric,missing,ANY-method}
\alias{[<-,splrMatrix,missing,numeric,ANY-method}
\alias{[<-,Matrix,ANY,ANY,splrMatrix-method}
\alias{str,splrMatrix-method}
\alias{t,splrMatrix-method}
\alias{diag,splrMatrix-method}
\title{"SPLR" Methods}
\usage{
\S4method{show}{splrMatrix}(object)

\S4method{print}{splrMatrix}(x)

\S4method{\%*\%}{splrMatrix,splrMatrix}(x, y)

\S4method{\%*\%}{splrMatrix,matrix_list}(x, y)

\S4method{\%*\%}{matrix_list,splrMatrix}(x, y)

\S4method{\%*\%}{Matrix,splrMatrix}(x, y)

\S4method{\%*\%}{matrix,splrMatrix}(x, y)

\S4method{\%*\%}{numeric,splrMatrix}(x, y)

\S4method{\%*\%}{numLike,splrMatrix}(x, y)

\S4method{\%*\%}{ANY,splrMatrix}(x, y)

\S4method{dim}{splrMatrix}(x)

\S4method{length}{splrMatrix}(x)

\S4method{\%*\%}{splrMatrix,Matrix}(x, y)

\S4method{\%*\%}{splrMatrix,matrix}(x, y)

\S4method{\%*\%}{splrMatrix,numeric}(x, y)

\S4method{\%*\%}{splrMatrix,numLike}(x, y)

\S4method{\%*\%}{splrMatrix,ANY}(x, y)

\S4method{*}{splrMatrix,splrMatrix}(e1, e2)

\S4method{*}{Matrix,splrMatrix}(e1, e2)

\S4method{*}{splrMatrix,ddiMatrix}(e1, e2)

\S4method{*}{ddiMatrix,splrMatrix}(e1, e2)

\S4method{*}{matrix,splrMatrix}(e1, e2)

\S4method{*}{numeric,splrMatrix}(e1, e2)

\S4method{*}{ANY,splrMatrix}(e1, e2)

\S4method{*}{splrMatrix,matrix}(e1, e2)

\S4method{*}{splrMatrix,Matrix}(e1, e2)

\S4method{*}{splrMatrix,numeric}(e1, e2)

\S4method{*}{splrMatrix,ANY}(e1, e2)

\S4method{/}{splrMatrix,matrix}(e1, e2)

\S4method{/}{splrMatrix,Matrix}(e1, e2)

\S4method{/}{splrMatrix,ANY}(e1, e2)

\S4method{+}{splrMatrix,splrMatrix}(e1, e2)

\S4method{-}{splrMatrix,splrMatrix}(e1, e2)

\S4method{+}{splrMatrix,Matrix}(e1, e2)

\S4method{+}{splrMatrix,numeric}(e1, e2)

\S4method{+}{splrMatrix,ANY}(e1, e2)

\S4method{-}{splrMatrix,missing}(e1, e2 = NULL)

\S4method{-}{splrMatrix,Matrix}(e1, e2)

\S4method{-}{splrMatrix,ddiMatrix}(e1, e2)

\S4method{-}{splrMatrix,numeric}(e1, e2)

\S4method{-}{splrMatrix,ANY}(e1, e2)

\S4method{+}{Matrix,splrMatrix}(e1, e2)

\S4method{+}{numeric,splrMatrix}(e1, e2)

\S4method{+}{ANY,splrMatrix}(e1, e2)

\S4method{-}{Matrix,splrMatrix}(e1, e2)

\S4method{-}{numeric,splrMatrix}(e1, e2)

\S4method{-}{ANY,splrMatrix}(e1, e2)

\S4method{norm}{splrMatrix,character}(x, type, ...)

\S4method{rowSums}{splrMatrix}(x, na.rm = FALSE, dims = 1, ...)

\S4method{colSums}{splrMatrix}(x, na.rm = FALSE, dims = 1, ...)

\S4method{rowMeans}{splrMatrix}(x, na.rm = FALSE, dims = 1, ...)

\S4method{colMeans}{splrMatrix}(x, na.rm = FALSE, dims = 1, ...)

\S4method{sum}{splrMatrix,ANY}(x, ..., na.rm = FALSE)

\S4method{mean}{splrMatrix}(x, ...)

\S4method{[}{splrMatrix,missing,missing,missing}(x, i = NULL, j = NULL, drop = NULL)

\S4method{[}{splrMatrix,numeric,numeric,logical}(x, i, j, ..., drop = TRUE)

\S4method{[}{splrMatrix,missing,numeric,logical}(x, i, j, ..., drop = TRUE)

\S4method{[}{splrMatrix,missing,numeric,missing}(x, i, j, ..., drop = TRUE)

\S4method{[}{splrMatrix,missing,logical,logical}(x, i, j, ..., drop = TRUE)

\S4method{[}{splrMatrix,missing,logical,missing}(x, i, j, ..., drop = TRUE)

\S4method{[}{splrMatrix,numeric,missing,logical}(x, i, j, ..., drop = TRUE)

\S4method{[}{splrMatrix,numeric,missing,missing}(x, i, j, ..., drop = TRUE)

\S4method{[}{splrMatrix,logical,missing,logical}(x, i, j, ..., drop = TRUE)

\S4method{[}{splrMatrix,logical,missing,missing}(x, i, j, ..., drop = TRUE)

\S4method{[}{splrMatrix,numeric,ANY,logical}(x, i, j, ..., drop = TRUE)

\S4method{[}{splrMatrix,numeric,logical,logical}(x, i, j, ..., drop = TRUE)

\S4method{[}{splrMatrix,numeric,ANY,missing}(x, i, j, ..., drop = FALSE)

\S4method{[}{splrMatrix,logical,ANY,ANY}(x, i, j, ..., drop = TRUE)

\S4method{[}{splrMatrix,logical,ANY,missing}(x, i, j, ..., drop = TRUE)

\S4method{[}{splrMatrix,logical,numeric,missing}(x, i, j, ..., drop = TRUE)

\S4method{[}{splrMatrix,logical,numeric,logical}(x, i, j, ..., drop = TRUE)

\S4method{[}{splrMatrix,missing,missing,missing}(x, i = NULL, j = NULL, drop = NULL)

\S4method{[}{splrMatrix,matrix,missing,missing}(x, i, j, ..., drop = TRUE)

\S4method{[}{splrMatrix,numeric,numeric,ANY}(x, i, j, ...) <- value

\S4method{[}{splrMatrix,numeric,missing,ANY}(x, i, j, ...) <- value

\S4method{[}{splrMatrix,missing,numeric,ANY}(x, i, j, ...) <- value

\S4method{[}{Matrix,ANY,ANY,splrMatrix}(x, i, j, ...) <- value

\S4method{dim}{splrMatrix}(x)

\S4method{str}{splrMatrix}(object)

\S4method{t}{splrMatrix}(x)

\S4method{diag}{splrMatrix}(x)
}
\arguments{
\item{object}{As in "Matrix"}

\item{x}{As in "Matrix"}

\item{y}{As in "Matrix"}

\item{e1}{As in "Matrix"}

\item{e2}{As in "Matrix"}

\item{type}{As in "Matrix"}

\item{...}{As in "Matrix"}

\item{na.rm}{As in "Matrix"}

\item{dims}{As in "Matrix"}

\item{i}{As in "Matrix"}

\item{j}{As in "Matrix"}

\item{drop}{As in "Matrix"}

\item{value}{As in "Matrix"}
}
\value{
Results of matrix operations for splrMatrix objects.
 Attempts are made such that the returned object is stored efficiently,
 either as a splrMatrix or sparse Matrix.
}
\description{
Methods for the splrMatrix class. Most behave like
Matrix methods though things like output show the
decomposition. Use as.matrix to see the computed
dense matrix.
}
\keyword{internal}
