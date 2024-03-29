# Measures for computing the goodness of matching for each vertex.

row_cor <- function(g1,g2){
  g1 <- g1[]
  g2 <- g2[]


  sapply(1:nrow(g1),
    function(v) suppressWarnings(1-stats::cor(g1[v,],g2[v,])))
}

row_diff <- function(g1,g2){
  g1 <- g1[]
  g2 <- g2[]
  g1 <- as.matrix(g1)
  g2 <- as.matrix(g2)
  rowSums(abs(g1-g2))
}

row_perm_stat <- function(g1,g2,exact=TRUE){
  g1 <- g1[]
  g2 <- g2[]

  if(exact){
    m <- mean_row_diff(g1,g2)
    v <- var_row_diff(g1,g2)
  } else {
    mv <- row_diff_perm(g1,g2)
    m <- mv$mean
    v <- mv$var
  }

  d <- rowSums(abs(g1-g2))

  (d-m)/sqrt(v)
}

row_diff_perm <- function(g1, g2, nmc = 1000, sym=FALSE){
  g1 <- g1[]
  g2 <- g2[]

  n <- nrow(g2)
  A <- Matrix(0,n,nmc)

  for(mc in 1:nmc){
    p <- sample(n)
    A[,mc] <- rowSums(abs(g1[,order(p)]-g2[p,]))
  }
  m <- rowMeans(A)
  v <- apply(A,1,stats::var)
  if(sym){
    mv <- row_diff_perm(g1,g2,nmc)
    m <- m+mv$mean
    v <- sqrt(v*mv$var)
  }
  list(mean=m,var=v)
}

mean_row_diff <- function(g1, g2, sym=FALSE){
  g1 <- g1[]
  g2 <- g2[]

  dg1 <- rowSums(g1)
  dg2 <- rowSums(g2)
  mdg2 <- mean(dg2)
  n <- nrow(g1)

  r1 <- mdg2
  r2 <- dg1 * (1 - 2 * mdg2 / (n - 1))
  ED <- r1 + r2
  if (sym){
    ED <- (ED + mean_row_diff(g2, g1)) / 2
  }
  ED
}

var_row_diff <- function(g1, g2, sym=FALSE){
  g1 <- g1[]
  g2 <- g2[]

  dg1 <- rowSums(g1)
  dg2 <- rowSums(g2)
  mdg2 <- mean(dg2)
  n <- nrow(g1)

  m2dg2 <- mean(dg2^2)
  vdg2 <- m2dg2-mdg2^2 # don't use var here

  v1 <- (1-2*dg1/(n-1))^2*vdg2
  v2 <- 4*dg1*(n-1-dg1)*((n-1)*mdg2-m2dg2)/((n-1)^2*(n-2))
  VD <- v1+v2

  if(sym){
    VD <- sqrt(VD*var_row_diff(g2,g1))
  }
  VD
}
