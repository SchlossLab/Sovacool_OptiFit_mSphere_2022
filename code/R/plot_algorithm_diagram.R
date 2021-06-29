# Adapted from Pat's code
# source: https://github.com/SchlossLab/Westcott_OptiClust_mSphere_2017/blob/a8bc26855423bba85acc0b8e7cca075e5c94f533/submission/supplemental_text.Rmd#L26-L28
library(DiagrammeR)
library(diagram)
library(Matrix)
library(tidyverse)
set.seed(20200308)

# blank placeholder plot for now
ggplot() + theme_void() + annotate('text', x=1, y=1, 
                                   label='insert algorithm diagram here', 
                                   size=8, angle=45)
dims <- eval(parse(text=snakemake@params[['dim']]))
ggsave(snakemake@output[['tiff']],
       device = 'tiff', dpi=300,
       width=dims[1], height=dims[2], units='in')

# Pat's MCC function
# source: https://github.com/SchlossLab/Westcott_OptiClust_mSphere_2017/blob/a8bc26855423bba85acc0b8e7cca075e5c94f533/submission/supplemental_text.Rmd#L26-L28
mcc <- function(tp, tn, fp, fn) {
  format(round((tp * tn - fp * fn) / 
                 sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)), 
               digits = 2),
  digits = 2,
  nsmall = 2L
  )
}

calc_mcc_from_conf_mat <- function(conf_mat) {
  tp <- conf_mat[1, 1]
  tn <- conf_mat[2, 2]
  fp <- conf_mat[1, 2]
  fn <- conf_mat[2, 1]
  return(mcc(tp, tn, fp, fn))
}

get_seqs_to_otus <- function(otu_list) {
  seqs_to_otus <- list()
  for (i in 1:length(otu_list)) {
    otu <- otu_list[[i]]
    for (seq in otu) {
      seqs_to_otus[seq] <- i
    }
  }
  return(seqs_to_otus)
}


get_otu_label <- function(otus, idx) {
  x <- otus[[idx]]
  paste0(x, collapse = ",")
}

# assignments of sequences in OTUs as matrix
get_otu_mat_from_list <- function(otu_list, seqs) {
  otu_mat <- matrix(
    nrow = length(seqs), ncol = length(seqs),
    dimnames = list(seqs, seqs), data = 0
  )
  for (i in 1:length(otu_list)) {
    otu <- otu_list[[i]]
    for (seq1 in otu) {
      for (seq2 in otu) {
        otu_mat[seq1, seq2] <- 1
        otu_mat[seq2, seq1] <- 1
      }
    }
  }
  return(otu_mat)
}

get_condition <- function(pair) {
  seq1 <- pair[1]
  seq2 <- pair[2]
  in_dist <- dist_mat[seq1, seq2]
  in_otu <- otu_mat[seq1, seq2]
  return(c(actual = in_dist, pred = in_otu))
}

# compare otu_mat to dist_mat to get tp/tn/fp/fn values
get_conf_mat <- function(pairs, dist_mat, otu_mat) {
  conditions <- apply(seq_pairs, 2, get_condition)
  conf_mat <- caret::confusionMatrix(table(
    conditions["pred", ],
    conditions["actual", ]
  )) %>%
    as.matrix()
  return(conf_mat)
}

pat_example <- function() {

# example data
seqs <- LETTERS[1:17]
seq_pairs <- combn(seqs, 2)
# distances within threshold
dat <- data.frame(
  Seq1 = c("D", "F", "G", "H", "I", "I", "J", "J", "N", "O", "P", "P", "P", "Q", "Q"),
  Seq2 = c("B", "E", "C", "A", "B", "D", "A", "H", "M", "L", "K", "L", "O", "E", "F"),
  Distance = c(0.024, 0.028, 0.028, 0.027, 0.016, 0.024, 0.028, 0.024, 0.024, 0.024, 0.016, 0.027, 0.027, 0.024, 0.028)
)

# matrix of pairwise sequences within distance threshold
dist_mat <- matrix(
  nrow = length(seqs), ncol = length(seqs),
  dimnames = list(seqs, seqs), data = 0
)
diag(dist_mat) <- 1
for (i in 1:nrow(dat)) {
  r <- dat[i, ]
  seq1 <- r[["Seq1"]]
  seq2 <- r[["Seq2"]]
  dist_mat[seq1, seq2] <- 1
  dist_mat[seq2, seq1] <- 1
}

# example otu assignment
otu_list <- list(
  c("B", "D", "I"),
  c("E", "F"),
  c("C", "G"),
  c("A", "H", "J"),
  c("M", "N"),
  c("L", "O", "P"),
  c("K"),
  c("Q")
)
seqs_to_otus <- get_seqs_to_otus(otu_list)
otu_mat <- get_otu_mat_from_list(otu_list, seqs)

current_mcc <- calc_mcc_from_conf_mat(get_conf_mat(seq_pairs, dist_mat, otu_mat))

# Pat's method for plotting network graph
tp <- 0
tn <- 1210
fp <- 0
fn <- 15
names <- c("B,D,I", "E,F", "C,G", "A,H,J", "L", "M", "N", "O", "K,P", "Q", "...")
n_seqs <- length(names)
M <- matrix(nrow = n_seqs, ncol = n_seqs, byrow = TRUE, data = 0)
rownames(M) <- colnames(M) <- names
C <- matrix(nrow = n_seqs, ncol = n_seqs, byrow = TRUE, data = 0)
rownames(C) <- colnames(C) <- names
M["L", "L"] <- mcc(tp = tp, tn = tn, fp = fp, fn = fn)
M["O", "L"] <- mcc(tp = tp + 1, tn = tn, fp = fp, fn = fn - 1)
M["K,P", "L"] <- mcc(tp = tp + 1, tn = tn - 1, fp = fp + 1, fn = fn - 1)
C["O", "L"] <- -0.4
C["K,P", "L"] <- 0.4
par(mar = c(1, 2, 1, 2), xpd = T)
plotmat(M, pos = n_seqs, name = names, lwd = 1, curve = C, box.lwd = 2, 
        cex.txt = 0.6, box.size = 0.03, box.type = "circle", shadow.size = 0, 
        arr.type = "triangle", dtext = 0.5, self.shiftx = -0.02, 
        self.shifty = -0.04, xpd = T, box.cex = 0.6, arr.length = 0.3)
}

# example diagrammeR graphs
create_graph(nodes_df=create_node_df(4), 
             edges_df = create_edge_df(from = 1:3, to = (2:4)), 
             attr_theme = 'lr') %>% 
  render_graph()

from_adj_matrix(otu_mat, use_diag = FALSE) %>% render_graph()

from_adj_matrix(dist_mat, use_diag = FALSE) %>% render_graph()

from_adj_matrix(otu_mat, use_diag = FALSE) %>% 
  add_global_graph_attrs('rankdir', 'TB', attr_type = 'graph') %>% 
  add_global_graph_attrs('layout', 'dot', attr_type = 'graph') %>% 
  render_graph()

# can I make part of a node label bold?
library(ggraph)
library(tidygraph)
graph <- tbl_graph(nodes = data.frame(name = c("abc", "def", "g**h**i", "jkl")), 
               edges = data.frame(from = c(1, 1, 1, 2, 3, 3, 4, 4, 4),
                             to = c(2, 3, 4, 1, 1, 2, 1, 2, 3)))
ggraph(graph, layout = 'fr') + 
  geom_edge_link() + 
  geom_node_label(aes(label = name)) + 
  #geom_richtext(aes(label = name)) +
  labs(title = 'an *example* plot') +
  theme(plot.title = element_markdown())

# have to use geom_richtext() to get markdown in labels

nodes <- data.frame(
  otu_label = c("abc", "def", "g**h**i", "jkl"),
  x = c(1, 2, 3, 4),
  y = c(1, 1, 1, 1)
) 
edges <- data.frame(from = c(2, 2),
                    to   = c(1, 3)
                    ) %>% 
  mutate(x1 = nodes[['x']][from],
         y1 = nodes[['y']][from],
         x2 = nodes[['x']][to],
         y2 = nodes[['y']][to])


ggplot() +
  geom_richtext(aes(x, y, label = otu_label), data = nodes) +
  geom_curve(
    aes(x = x1, y = y1, xend = x2, yend = y2),
    data = edges,
    arrow = arrow(length = unit(0.03, "npc"))
  )
ggsave('figures/tmp.tiff', device = 'tiff', dpi=300, 
       width = 10, heigh = 2, units = 'in')
