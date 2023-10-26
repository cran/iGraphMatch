## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(warning=FALSE, message=FALSE)
options(scipen = 1, digits = 2)

## ----sample_graph-------------------------------------------------------------
library(iGraphMatch)
set.seed(1)
cgnp_pair <- sample_correlated_gnp_pair(n = 5, corr = 0.7, p = 0.5)
(cgnp_g1 <- cgnp_pair$graph1)
cgnp_g1[]
cgnp_g2 <- cgnp_pair$graph2

## ----sample_ieg---------------------------------------------------------------
set.seed(123)
p <- matrix(runif(5^2, .5, .8),5)
c <- matrix(runif(5^2, .5, .8),5)
ieg_pair <- sample_correlated_ieg_pair(n = 5, p_mat = p, c_mat = c)

## ----sample_sbm---------------------------------------------------------------
pm <- cbind(c(.7, .001), c(.001, .5))
sbm_pair <- sample_correlated_sbm_pair(n = 5, pref.matrix = pm,
                                       block.sizes = c(2,3), corr = 0.5)

## ----center_graph-------------------------------------------------------------
center_graph(cgnp_g1, scheme = "center", use_splr = FALSE)
center_graph(cgnp_g1, scheme = 2, use_splr = FALSE)

## ----self_defined_method------------------------------------------------------
graph_match_rand <- function(A, B, seeds = NULL, 
                             similarity = NULL, rand_seed){
  totv1 <- nrow(A[[1]])
  totv2 <- nrow(B[[1]])
  nv <- max(totv1, totv2)

  set.seed(rand_seed)
  corr <- data.frame(corr_A = 1:nv, 
                     corr_B = c(1:nv)[sample(nv)])

  graphMatch(
    corr = corr,
    nnodes = c(totv1, totv2),
    detail = list(
      rand_seed = rand_seed
    )
  )
}

match_rand <- gm(cgnp_g1, cgnp_g2, 
                 method = graph_match_rand, rand_seed = 123)

## ----seeds--------------------------------------------------------------------
hard_seeds <- 1:5 <= 2
soft_seeds <- data.frame(seed_A = 3, seed_B = 4)

## ----start_bari---------------------------------------------------------------
as.matrix(start_bari <- init_start(start = "bari", nns = 3,
      ns = 2, soft_seeds = soft_seeds))

## ----start_rds----------------------------------------------------------------
set.seed(1)
as.matrix(start_rds <- init_start(start = "rds", nns = 3,
      ns = 2, soft_seeds = soft_seeds))

## ----match_rds----------------------------------------------------------------
match_rds <- gm(cgnp_g1, cgnp_g2, seeds = hard_seeds,
                method = "indefinite", start = "rds")

## ----match_convex-------------------------------------------------------------
set.seed(123)
match_convex <- gm(cgnp_g1, cgnp_g2, seeds = hard_seeds,
                   method = "indefinite", start = "convex")

## ----match_perco--------------------------------------------------------------
sbm_g1 <- sbm_pair$graph1
sbm_g2 <- sbm_pair$graph2
match_perco <- gm(sbm_g1, sbm_g1, seeds = hard_seeds, 
                  method = "percolation", r = 2)
match_perco

## ----match_perco_multi, warning=FALSE-----------------------------------------
matrix_lA <- list(sbm_g1, ieg_pair$graph1, cgnp_g1)
matrix_lB <- list(sbm_g2, ieg_pair$graph2, cgnp_g2)
match_perco_list <- gm(A = matrix_lA, B = matrix_lB, seeds = hard_seeds, 
                       method = "percolation", r = 2)
match_perco_list

## ----match_IsoRank------------------------------------------------------------
set.seed(1)
sim <- as.matrix(init_start(start = "bari", nns = 5, 
                            soft_seeds = hard_seeds))
match_IsoRank <- gm(A = matrix_lA, B = matrix_lB, 
                    seeds = hard_seeds, similarity = sim, 
                    method = "IsoRank", lap_method = "LAP")

## ----graphMatch_class---------------------------------------------------------
match_convex@corr
match_convex@call
match_convex@nnodes

## ----match_nonseeds-----------------------------------------------------------
match_convex[!match_convex$seeds]

## ----permutation_matrix-------------------------------------------------------
match_convex[]

## ----graphMatch_multilicity---------------------------------------------------
match_convex %*% cgnp_g2

## ----summary_convex-----------------------------------------------------------
summary(match_convex, cgnp_g1, cgnp_g2, true_label = 1:5)

## ----summary_IsoRank----------------------------------------------------------
summary(match_IsoRank, matrix_lA, matrix_lB)

## ----best_matches-------------------------------------------------------------
best_matches(cgnp_g1, cgnp_g2, match = match_convex, 
             measure = "row_perm_stat", num = 3, 
             true_label = 1:igraph::vcount(cgnp_g1))

## ----visualization, out.width="47.5%", fig.show="hold", fig.cap=" Match visualizations. Grey, blue, and red colors indicate common edges, missing edges present only in the first network, and extra edges present only in the second network, respectively."----
plot(cgnp_g1, cgnp_g2, match_convex)
plot(cgnp_g1[], cgnp_g2[], match_convex)

## ----dataset-overview, echo=FALSE, results="asis"-----------------------------
summarize_graph <- function(gp, name) {
  nn <- nrow(gp[[1]][])
  if (nn != nrow(gp[[2]][])) {
    nn <- paste(nn, "/", nrow(gp[[2]][]))
  }
  ne <- paste(igraph::ecount(gp[[1]]), "/", igraph::ecount(gp[[2]]))

  data.frame(
    `Dataset` = name,
    `Nodes` = as.character(nn),
    `Edges` = ne,
    `Correlation` = cor(as.numeric(gp[[1]][]), as.numeric(gp[[2]][])),
    `Weighted` = 
      ifelse(igraph::is_weighted(gp[[1]]), "Yes", "No")
    ,
    `Directed` = 
      ifelse(igraph::is_directed(gp[[1]]), "Yes", "No")
    ,
    `Loop` = paste(
      ifelse(any(igraph::is.loop(gp[[1]])), "Yes", "No"),
      "/",
      ifelse(any(igraph::is.loop(gp[[2]])), "Yes", "No")
    )
  )
}

knitr::kable(
  rbind(
    summarize_graph(Enron, "Enron"),
    summarize_graph(C.Elegans, "C. Elegans")
  ),
  col.names = c("Dataset", "\\# Nodes", "\\# Edges", "Correlation", "Weighted", "Directed", "Loop"),
  booktabs = TRUE,
  caption =  "Overview of the Enron and C. Elegans graphs.",
  escape = FALSE
)

## ----edge-summary, echo=FALSE, results="asis"---------------------------------
summarize_edge <- function(gp1, gp2, name) {
  nn <- nrow(gp1[])
  corr <- graphMatch(data.frame(corr_A = 1:nn, corr_B = 1:nn), nn)
  edge_info <- summary(corr, gp1, gp2)$edge_match_info

  data.frame(
    `Dataset` = name,
    `Common Edges` = edge_info$common_edges,
    `Missing Edges` = edge_info$missing_edges,
    `Extra Edges` = edge_info$extra_edges
  )
}

knitr::kable(
  rbind(
    summarize_edge(Enron[[1]], Enron[[2]], "Enron"),
    summarize_edge(C.Elegans[[1]], C.Elegans[[2]], "C. Elegans")
  ),
  col.names = c("Dataset", "Common edges", "Missing edges", "Extra edges"),
  booktabs = TRUE,
  caption =  "Edge summary under the true alignments of the Enron and C. Elegans graphs. The columns indicate the number of common edges, missing edges in $G_1$, and extra edges in $G_2$. For weighted graphs, we define a pair of corresponding edges as a common edge as long as they both have positive weights. ",
  escape = FALSE
)

## ----load_packages, message=FALSE---------------------------------------------
library(igraph)
library(iGraphMatch)
library(purrr)
library(dplyr)

## ----Enron-graph, fig.show="hold", fig.cap="Asymmetric adjacency matrices of aligned Enron Corporation communication networks. The vertices are sorted by a community detection algorithm (@community_detection) and degree."----
g <- igraph::as.undirected(Enron[[1]])
com <- igraph::membership(igraph::cluster_fast_greedy(g))
deg <- rowSums(as.matrix(g[]))
ord <- order(max(deg)*com+deg)
plot(Enron[[1]][][ord,ord], Enron[[2]][][ord,ord])

## ----enron_lcc----------------------------------------------------------------
vid1 <- which(largest_cc(Enron[[1]])$keep)
vid2 <- which(largest_cc(Enron[[2]])$keep)

vinsct <- intersect(vid1, vid2) 
v1 <- setdiff(vid1, vid2)
v2 <- setdiff(vid2, vid1)
A <- Enron[[1]][][c(vinsct, v1), c(vinsct, v1)]
B <- Enron[[2]][][c(vinsct, v2), c(vinsct, v2)]

## ----enron_match_fw-----------------------------------------------------------
set.seed(1)
match_FW <- gm(A = A, B = B, start = "bari", max_iter = 200)
head(match_FW)

## ----enron_summary------------------------------------------------------------
summary(match_FW, A, B)

## ----compute, echo=FALSE------------------------------------------------------
emi <- summary(match_FW, A, B)$edge_match_info

## ----enron_center, warning=TRUE, message=TRUE---------------------------------
A_center <- center_graph(A = A, scheme = "naive", use_splr = TRUE)
B_center <- center_graph(A = B, scheme = "center", use_splr = TRUE)
set.seed(1)
match_FW_center <- gm(A = A_center, B = B_center, 
                           start = "bari", max_iter = 200)
summary(match_FW_center, A, B)

## ----enron_best_matches-------------------------------------------------------
bm <- best_matches(A = A, B = B, match = match_FW_center, 
             measure = "row_perm_stat")
head(bm)

## ----enron_match_w_hard_seeds, warning=FALSE, message=FALSE, results='hide'----
match_w_hard_seeds <- function(ns){
  seeds_bm <- head(bm, ns)
  precision <- mean(seeds_bm$A_best == seeds_bm$B_best)
  match_FW_center_seeds <- gm(A = A_center, B = B_center,
                           seeds = seeds_bm, similarity = NULL,
                           start = "bari", max_iter = 100)
  edge_info <- summary(match_FW_center_seeds, A, B)$edge_match_info
  cbind(ns, precision, edge_info)
}
set.seed(12345)
map_dfr(seq(from = 0, to = 80, by = 20), match_w_hard_seeds)

## ----enron_match_w_hard_seeds_table, warning=FALSE, message=FALSE, echo=FALSE----
set.seed(12345)
map_dfr(seq(from = 0, to = 80, by = 20), match_w_hard_seeds) %>% 
  knitr::kable(col.names = c("ns", "precision","common", "missing", "extra", "fnorm"), 
               booktabs = TRUE, digits = 2)

## ----enron_match_w_soft_seeds, warning=FALSE, message=FALSE, results='hide'----
match_w_soft_seeds <- function(ns){
  seeds_bm <- head(bm, ns)
  precision <- mean(seeds_bm$A_best == seeds_bm$B_best)
  start_soft <- init_start(start = "bari", 
                           nns = max(dim(A)[1], dim(B)[1]), 
                           soft_seeds = seeds_bm)
  match_FW_center_soft_seeds <- gm(A = A_center, B = B_center, 
                           start = start_soft, max_iter = 100)
  edge_info <- summary(match_FW_center_soft_seeds, A, B)$edge_match_info
  cbind(ns, precision, edge_info)
}
set.seed(12345)
map_dfr(seq(from = 0, to = 80, by = 20), match_w_soft_seeds)

## ----enron_match_w_soft_seeds_table, warning=FALSE, message=FALSE, echo=FALSE----
set.seed(12345)
map_dfr(seq(from = 0, to = 80, by = 20), match_w_soft_seeds) %>% 
  knitr::kable(col.names = c("ns", "precision","common", "missing", "extra", "fnorm"), 
               booktabs = TRUE, digits = 2)

## ----enron_core---------------------------------------------------------------
nc <- length(vinsct)
nj <- max(length(v1), length(v2))
core_precision <- 1:nc %>% map_dbl(~mean(bm$A_best[1:.x]<=nc))
junk_precision <- 1:nj %>% map_dbl(~mean(bm$A_best[(nc+.x):(nc+nj)]>nc))

## ----echo=FALSE, out.width="0.5\\textwidth", out.height="0.24\\textheight", fig.cap="\\label{fig:core}Mean precision for identifying core and junk vertices for the Enron networks by using the row permutation test. The vertical lines separate the performance of identifying core vertices with low ranks from junk vertices with high ranks. The horizontal lines indicate the performance of a random classifier."----
library(ggplot2)

map_df_core <- data.frame(rank = 1:nc, which = "core", 
                          mean_precision = core_precision)
map_df_junk <- data.frame(rank = (nc+1):(nc+nj), which = "junk", 
                          mean_precision = junk_precision)
map_df <- rbind(map_df_core, map_df_junk)
map_df <- map_df %>% 
  mutate(random = ifelse(which == "junk", nj/(nc+nj), nc/(nc+nj))) %>%
  mutate(start=ifelse(which=="core",0,nc-5),
         end=ifelse(which=="core",nc+5,nc+nj))

map_df  %>% 
  ggplot(aes(x=rank,y=mean_precision)) +
  geom_line(alpha=.5,size=.3) +
  geom_segment(aes(y=random,yend=random, x=start,xend=end),
               data=map_df,color="black",alpha=1,linetype=2,size=.1)+
  geom_vline(aes(xintercept=nc),size=.1,linetype=2)+
  ylab("mean precision")+
  theme_minimal()+
  theme(text=element_text(family="Palatino"),
        axis.line=element_line(size=.3,color="black"),
        panel.grid.major=element_line(color="grey95",size=0),
        panel.grid.minor=element_line(color="grey98",size=0))

## ----C-Elegans-edge, warning = FALSE, fig.show="hold", fig.cap="Edge discrepancies for the matched graphs with the true correspondence (left) and FW algorithm starting at the true correspondence (right). Green pixels represents an edge in the chemical graph while no edge in the electrical graph. Red pixels represent only an edge in the electrical graph. Grey pixels represent there is an edge in both graphs and white represents no edge in both graphs."----
C1 <- C.Elegans[[1]][] > 0
C2 <- C.Elegans[[2]][] > 0
plot(C1[], C2[])
match <- gm(C1, C2, start = Matrix::Diagonal(nrow(C1)))
plot(C1[], C2[], match)

## ----c.el_emi-----------------------------------------------------------------
nv <- nrow(C1)
id_match <- graphMatch(data.frame(corr_A = 1:nv, corr_B = 1:nv), nv)
i_sum <- summary(id_match, C.Elegans[[1]], C.Elegans[[2]])
m_sum <- summary(match, C.Elegans[[1]], C.Elegans[[2]], id_match)
i_emi <- i_sum$edge_match_info
m_emi <- m_sum$edge_match_info

## ----celegans_soft, warning=FALSE, message=FALSE------------------------------
seeds <- sample(nrow(C1), 20)
sim <- init_start(start = "bari", nns = nrow(C1), soft_seeds = seeds)

## ----celegans_match, warning=FALSE--------------------------------------------
set.seed(123)
m_FW <- gm(A = C1, B = C2, seeds = seeds, 
           similarity = sim, method = "indefinite",
           start = "bari", max_iter = 100)
m_PATH <- gm(A = C1, B = C2, seeds = seeds,
             similarity = NULL, method = "PATH",
             epsilon = 1, tol = 1e-05)
m_Iso <- gm(A = C1, B = C2, seeds = seeds,
             similarity = as.matrix(sim), method = "IsoRank",
             max_iter = 50, lap_method = "LAP")

## ----celegans_match_eval------------------------------------------------------
match_eval <- function(match){
  precision <- mean(match$corr_A == match$corr_B) 
  order <- apply(match$soft, MARGIN = 1, FUN = order, decreasing = TRUE)
  top3 <- t(order[1:3,]) - 1:ncol(order) 
  MAP3 <- mean(apply(top3, MARGIN = 1, FUN = function(v){0 %in% v}))
  
  round(data.frame(precision, MAP3),4)
}

sapply(list(m_FW, m_PATH, m_Iso), match_eval) %>% 
  knitr::kable(col.names = c("Frank Wolfe", "PATH", "IsoRank"), 
               booktabs = TRUE, digits = 2)

## ----transpo_info, echo=FALSE-------------------------------------------------
library(Matrix)

tm <- Transportation[[1]]
cm <- Transportation[[2]]
candidate <- Transportation[[3]]

tn <- nrow(tm[[1]])
wn <- nrow(cm[[1]])
mc <- candidate %>%
  with(Matrix::sparseMatrix(i = tem, j = wor, x = 1,
                            dims = c(tn,wn)))
c_count <- rowSums(mc)

## ----transpo_set_up_graphs, eval=FALSE----------------------------------------
#  tm <- Transportation[[1]]
#  cm <- Transportation[[2]]
#  candidate <- Transportation[[3]]

## ----edge-summary-trans, echo=FALSE, results="asis"---------------------------
summarize_graph_trans <- function(gp1, gp2, name) {
  nn <- nrow(gp1)
  if (nn != nrow(gp2)) {
    nn <- paste(nn, "/", nrow(gp2))
  }
  ne <- paste(sum(gp1), "/", sum(gp2))

  data.frame(
    `Layer` = name,
    `Nodes` = as.character(nn),
    `Edges` = ne,
    `Correlation` = cor(as.numeric(gp1), as.numeric(gp2[1:nrow(gp1), 1:nrow(gp1)]))
  )
}

 
summary_table <- 
  rbind(
    summarize_graph_trans(tm[[1]], cm[[1]], "Ferry"),
    summarize_graph_trans(tm[[2]], cm[[2]], "Rail"),
    summarize_graph_trans(tm[[3]], cm[[3]], "Metro"),
    summarize_graph_trans(tm[[4]], cm[[4]], "Coach"),
    summarize_graph_trans(tm[[5]], cm[[5]], "Bus")
  )
a <- tm
b <- cm %>% purrr::map(~.x[][1:tn, 1:tn])

edge_layer <- summary(graphMatch(data.frame(corr_A = 1:tn, corr_B = 1:tn), as.integer(tn)), a, b)$edge_match_info
# edge_layer$layer <- c("Ferry", "Rail", "Metro", "Coach", "Bus")

knitr::kable(
  cbind(summary_table, edge_layer[,2:4]),
  col.names = c("Layer",  "\\# Nodes", "\\# Edges", "Correlation","Common", "Missing", "Extra"),
  booktabs = TRUE,
  caption =  "Overview of the Britain Transportation Network layers. Correlation is calculted using the template graph and the aligned induced subgraph of the world graph. The final three columns indicate the number of common edges, missing edges, and extra edges in the aligned subgraph of the world graph.",
  escape = FALSE
)
# Edge summary of the Britain Transportation Network layers.  

## ----out.width="47.5%", fig.show="hold", echo=FALSE, fig.cap="\\label{Fig:trans_net} Visualization of the template graph (left) and the world graph (right) with corresponding vertices, both derived from the Britain Transportation network with five layers: ferry, rail, metro, coach, and bus. Edges represent transportation transactions and each color indicates a different means of transportation from a different layer of network."----

l <- 1:length(tm)

tg_sub <- l %>% map(~ tm[[.x]] %>% graph_from_adjacency_matrix)
wg_sub <- l %>% map(~ cm[[.x]][1:tn,1:tn] %>% graph_from_adjacency_matrix)

for(l in 1:5){
  tg_sub[[l]] <- set_edge_attr(tg_sub[[l]], "layer", value = l)
  wg_sub[[l]] <- set_edge_attr(wg_sub[[l]], "layer", value = l)
}

tg_colored <- igraph::union(tg_sub[[1]], tg_sub[[2]], tg_sub[[3]], tg_sub[[4]], tg_sub[[5]])
wg_colored <- igraph::union(wg_sub[[1]], wg_sub[[2]], wg_sub[[3]], wg_sub[[4]], wg_sub[[5]])

tg_layer <- rep(0, ecount(tg_colored))
wg_layer <- rep(0, ecount(wg_colored))
for(l in 1:5){
  tg_layer <- tg_layer + ifelse(is.na(edge_attr(tg_colored)[[l]]) | tg_layer != 0, 0, edge_attr(tg_colored)[[l]])
  wg_layer <- wg_layer + ifelse(is.na(edge_attr(wg_colored)[[l]]) | wg_layer != 0, 0, edge_attr(wg_colored)[[l]])
}

layout <- layout_nicely(wg_colored)
plot(tg_colored, edge.color = factor(tg_layer),
     vertex.size = 5, vertex.label.cex = .666,
     edge.curved = TRUE,
     edge.arrow.size = 0.2, layout = layout)

plot(wg_colored, edge.color = factor(wg_layer),
     vertex.size = 5, vertex.label.cex = .666,
     edge.curved = TRUE,
     edge.arrow.size = 0.2, layout = layout)

## ----set_up_similarity_matrix, echo=FALSE-------------------------------------
start <- mc %>%
  rbind(Matrix(0, nrow = wn - tn, ncol = wn))
start[1:tn, ] <- diag(1 / rowSums(start[1:tn, ])) %*% start[1:tn, ]
similarity <- start * 1e5


## ----transpo_match------------------------------------------------------------
match <- gm(A = tm, B = cm, similarity = similarity, 
            method = "percolation", r = 4)
summary(match, tm, cm)

