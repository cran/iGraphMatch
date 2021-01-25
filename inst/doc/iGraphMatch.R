## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(iGraphMatch)

## ---- warning=FALSE, message=FALSE--------------------------------------------
library(igraph)
library(iGraphMatch)

set.seed(8)
cgnp_pair <- sample_correlated_gnp_pair(n = 10, corr = .5, p = .5)
g1 <- cgnp_pair$graph1
g2 <- cgnp_pair$graph2

(seeds <- 1:10 <= 3)

## -----------------------------------------------------------------------------
bad_soft_seeds <- rbind(c(4,4),c(8,6))

## -----------------------------------------------------------------------------
nns <- 7
ns <- 3
(start_bari <- init_start(start = "bari", nns = nns, ns = ns, soft_seeds = bad_soft_seeds))
set.seed(123)
(start_rds <- init_start(start = "rds", nns = nns, ns = ns, soft_seeds = bad_soft_seeds))
(start_convex <- init_start(start = "convex", nns = nns, ns = ns, soft_seeds = bad_soft_seeds, 
                           A = g1, B = g2, seeds = seeds))

## -----------------------------------------------------------------------------
set.seed(123)
match_bari <- graph_match_FW(g1, g2, seeds, start = start_bari)
match_rds <- graph_match_FW(g1, g2, seeds, start = start_rds)
match_convex <- graph_match_FW(g1, g2, seeds, start = start_convex)
err_bari_bad <- mean(match_bari$corr$corr_A[!seeds] != match_bari$corr$corr_B[!seeds])
err_rds_bad <- mean(match_rds$corr$corr_A[!seeds] != match_rds$corr$corr_B[!seeds])
err_convex_bad <- mean(match_convex$corr$corr_A[!seeds] != match_convex$corr$corr_B[!seeds])

## -----------------------------------------------------------------------------
set.seed(123)
match_nss_bari <- graph_match_FW(g1, g2, seeds, start = "bari")
match_nss_rds <- graph_match_FW(g1, g2, seeds, start = "rds")
match_nss_convex <- graph_match_FW(g1, g2, seeds, start = "convex")

## ----echo=FALSE---------------------------------------------------------------
err_nss_bari <- mean(match_nss_bari$corr$corr_A[!seeds] != match_nss_bari$corr$corr_B[!seeds])
err_nss_rds <- mean(match_nss_rds$corr$corr_A[!seeds] != match_nss_rds$corr$corr_B[!seeds])
err_nss_convex <- mean(match_nss_convex$corr$corr_A[!seeds] != match_nss_convex$corr$corr_B[!seeds])

## -----------------------------------------------------------------------------
good_soft_seeds <- rbind(c(4,4),c(8,8))

## ----output="asis", echo=FALSE,warning=FALSE----------------------------------
nns <- 7
ns <- 3
start_bari <- init_start(start = "bari", nns = nns, ns = ns, soft_seeds = good_soft_seeds)
start_rds <- init_start(start = "rds", nns = nns, ns = ns, soft_seeds = good_soft_seeds)
start_convex <- init_start(start = "convex", nns = nns, ns = ns, soft_seeds = good_soft_seeds, 
                           A = g1, B = g2, seeds = seeds)

set.seed(123)
match_bari <- graph_match_FW(g1, g2, seeds, start = start_bari)
match_rds <- graph_match_FW(g1, g2, seeds, start = start_rds)
match_convex <- graph_match_FW(g1, g2, seeds, start = start_convex)
err_bari_good <- mean(match_bari$corr$corr_A[!seeds] != match_bari$corr$corr_B[!seeds])
err_rds_good <- mean(match_rds$corr$corr_A[!seeds] != match_rds$corr$corr_B[!seeds])
err_convex_good <- mean(match_convex$corr$corr_A[!seeds] != match_convex$corr$corr_B[!seeds])

bad_soft_seeds <- data.frame(bari=err_bari_bad,rds=err_rds_bad,convex=err_convex_bad)
good_soft_seeds <- data.frame(bari=err_bari_good,rds=err_rds_good,convex=err_convex_good)
non_soft_seeds <- data.frame(bari=err_nss_bari,rds=err_nss_rds,convex=err_nss_convex)
result <- rbind(good_soft_seeds,bad_soft_seeds,non_soft_seeds)
row.names(result) <- c("good soft seeds","bad soft seeds","non soft seeds")

knitr::kable(result, caption = "Matching Errors With Various Initialization Methods")


## -----------------------------------------------------------------------------
set.seed(5)
cgnp_pair <- sample_correlated_gnp_pair_w_junk(n = 50, corr = .5, p = .5, ncore = 45)
g1 <- cgnp_pair$graph1
g2 <- cgnp_pair$graph2

seeds <- 1 : 50 <= 10
core <- 1 : 50 <= 45
junk <- !core
non_seeds <- !seeds

match <- graph_match_FW(g1, g2, seeds, start = "rds")

## -----------------------------------------------------------------------------
r <- best_matches(A = g1, B = g2, match = match,
                  measure = "row_perm_stat", num = sum(non_seeds))$A_best
r

## ----echo=FALSE, output="asis", warning=FALSE---------------------------------
r <- r-10
junk <- junk[non_seeds]
core <- core[non_seeds]
nc <- sum(core)
nj <- sum(junk)
core_precision <- sapply(seq(nc), function(n) mean(core[r <= n ]))
junk_precision <- sapply(seq(nj), function(n) mean(junk[r > nc + nj- n ]))

core <- data.frame(k1=max(core_precision),k5=quantile(core_precision,6/7),k10=quantile(core_precision,5/7),k25=quantile(core_precision,2/7),k35=min(core_precision))
result <- rbind(core)
row.names(result) <- c("core identification precision")
knitr::kable(result, caption = "Summarization Table for Core Identification Precisions")

junk <- data.frame(k1=max(junk_precision),k2=quantile(junk_precision,2/5),k3=quantile(junk_precision,3/5),k4=quantile(junk_precision,4/5),k5=min(junk_precision))
result <- rbind(junk)
row.names(result) <- c("junk identification precision")
knitr::kable(result, caption = "Summarization Table for Junk Identification Precisions")

## ---- warning=FALSE, message=FALSE--------------------------------------------
set.seed(6)
pm <- cbind( c(.3, .5), c(.5, .7) )
sbm_pair <- sample_correlated_sbm_pair_w_junk(n = 50, pref.matrix = pm, rho = 0.5,
                                              block.sizes = c(15,35), 
                                              core.block.sizes = c(10,30))
g1 <- sbm_pair$graph1
g2 <- sbm_pair$graph2

seeds <- 1 : 50 <= 5
seeds[16:20] <- TRUE

## ---- warning=FALSE, message=FALSE--------------------------------------------
match <- graph_match_FW(g1, g2, seeds = seeds, start = "convex")
err <- mean(match$corr$corr_A[!seeds] != match$corr$corr_B[!seeds])

## ----warning=FALSE, message=FALSE---------------------------------------------
seeds_adp <- best_matches(A = g1, B = g2, match = match, 
                            measure = "row_perm_stat", num = 3)
seeds_adp

## -----------------------------------------------------------------------------
seeds <- rbind(as.matrix(check_seeds(seeds, nv = 50)$seeds), as.matrix(seeds_adp)[,1:2])
match_adp <- graph_match_FW(g1, g2, seeds=seeds, start = "convex")
seeds <- 1 : 50 <= 10
err_adp <- mean(match_adp$corr$corr_A[!seeds] != match_adp$corr$corr_B[!seeds])

## ----output="asis", echo=FALSE------------------------------------------------
knitr::kable(data.frame(err_orig=err, err_adp=err_adp), caption = "Matching Errors for Adaptive Seeding")

