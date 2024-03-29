# Codes below test Expand When Stuck

set.seed(12)
G <- sample_correlated_gnp_pair(n = 10, corr = .5, p = .5)
A <- G$graph1
B <- G$graph2
seeds <- seq(4)

test_that("matching correspondence between graph1 and graph2", {
  tt <- gm(A, B, seeds, method = "percolation", r = 2, ExpandWhenStuck = TRUE)
  expect_snapshot_output(print(tt))
})

test_that("test number of seeds", {
  expect_equal(
    sum(gm(A, B,seeds, method = "percolation", r = 2, ExpandWhenStuck = TRUE)$seeds),
    length(seeds)
  )
})

# with similarity score
sim <- matrix(rnorm(100), 10)
test_that("exp w. similarity score", {
  m <- gm(A, B, seeds = seeds, similarity = sim, method = "percolation", ExpandWhenStuck = TRUE)
  expect_equal(
    m[m$seeds],
    data.frame(corr_A = 1:4, corr_B = 1:4),
    ignore_attr = TRUE
  )
})

test_that("exp w. similarity score & no seeds", {
  m <- gm(A, B, seeds = NULL, similarity = sim, method = "percolation", ExpandWhenStuck = TRUE)
  expect_equal(sum(m$seeds), 0)
})

# directed graphs
set.seed(123)
g <- sample_correlated_gnp_pair(n = 10, corr = 0.5, p = 0.5, directed = TRUE)
A <- g$graph1
B <- g$graph2
test_that("perco w. similarity score", {
  m <- gm(A, B, seeds = seeds, similarity = sim, method = "percolation", ExpandWhenStuck = TRUE)
  expect_equal(
    m[m$seeds],
    data.frame(corr_A = 1:4, corr_B = 1:4),
    ignore_attr = TRUE
  )
})


# multiple candidate matches with same score
set.seed(12)
gp_list <- replicate(2, sample_correlated_gnp_pair(10, .5, .5), simplify = FALSE)
A <- lapply(gp_list, function(gp)gp[[1]])
B <- lapply(gp_list, function(gp)gp[[2]])
seeds <- 1:3
test_that("perco w. similarity score", {
  m <- gm(A, B, seeds = seeds, method = "percolation", ExpandWhenStuck = TRUE)
  expect_equal(
    m[m$seeds],
    data.frame(corr_A = 1:3, corr_B = 1:3),
    ignore_attr = TRUE
  )
})



