library(ape)
library(Biostrings)
library(ggplot2)
library(ggtree)
library(tidytree)
library(stringr)
library(dplyr)
library(bio3d)
args = commandArgs(trailingOnly=TRUE)
uniprot_id <- args[4]

aa_to_num <- function(aa) {
  amino_acids <- c("G", "A", "L", "M", "F", "W", "K", "Q", "E", "S", "P", "V", "I", "C", "Y", "H", "R", "N", "D", "T")
  num <- sapply(aa, function(a){ifelse(sum(amino_acids %in% a) == 1, as.numeric(which(amino_acids %in% a)), 21)})
  # num <- ifelse(sum(amino_acids %in% aa) == 1, as.numeric(which(amino_acids %in% aa)), 21)
  return(num)
}

num_to_aa <- function(num) {
  amino_acids <- c("G", "A", "L", "M", "F", "W", "K", "Q", "E", "S", "P", "V", "I", "C", "Y", "H", "R", "N", "D", "T")
  aa <- ifelse(num == 21, 21, amino_acids[num])
  return(aa)
}

compute_score <- function(file_nwk, file_fasta, output_name, human_id, pos_chosen) {
  
  # Read tree file
  tr_org <- read.tree(file_nwk)
  
  # Tree_info: node-node, node-leaf connections
  tree_info <- as.data.frame(as_tibble(tr_org))
  
  # Read fasta file, MSA
  fasta <- read.fasta(file = file_fasta)
  msa <- fasta$ali
  
  # connections_1: Parent node, connections_2: connected node/leaf
  connections_1 <- tree_info$parent
  connections_2 <- tree_info$node
  
  # Names of leaves
  names_all <- tr_org[["tip.label"]]
  msa <- msa[names_all, ]
  # Number of total leaves&nodes
  num_leaves <- length(tr_org[["tip.label"]])
  num_nodes <- tr_org[["Nnode"]]
  
  # Distance between leaves & nodes
  dd_node <- dist.nodes(tr_org)
  dist_leaf <- dd_node[1:num_leaves, 1:num_leaves]
  dist_node <- dd_node[(num_leaves+1):(num_leaves + num_nodes), (num_leaves+1):(num_leaves + num_nodes)]
  
  # Human position (leaf & node)
  h_name <- human_id
  human_codeml <- names_all[grep(pattern = h_name, x = names_all, fixed = TRUE)]
  leaf_human <- tree_info[which(tree_info$label == human_codeml), "node"]
  human_plc <- leaf_human
  node_human <- tree_info[which(tree_info$label == human_codeml), "parent"]
  nodes_raxml <- as.numeric(gsub(pattern = "Node", replacement = "", x = tree_info[num_leaves+1:num_nodes, "label"])) #Node or Branch
  names(nodes_raxml) <- tree_info[num_leaves+1:num_nodes, "node"]
  
  # Total number of positions from ancestralProbs file
  total_pos <- length(msa[1,])
  
  # Chosen positions (all or some)
  if (pos_chosen[1] == "all"){
    positions <- 1:total_pos
    score_all <- matrix(0, total_pos, 21)
  } else {
    positions <- pos_chosen
    score_all <- matrix(0, length(positions), 21)
  }
  
  # Connections between leaves & nodes
  chosen_leaves <- tree_info[1:num_leaves,c("parent", "node")]
  # Connections between nodes & nodes
  chosen_nodes <- tree_info[(num_leaves+2):(num_leaves +num_nodes),c("parent", "node")]
  leaf_names <- tree_info$label
  
  human_leaf_len <- as.double(tree_info[human_plc, "branch.length"])
  
  d_n <- dist_node[as.character(node_human),] + human_leaf_len
  d_l <- dist_leaf[leaf_human,]
  
  # chosen_nodes2: ordered connections (for probability differences)
  chosen_nodes2 <- matrix(0, num_nodes-1, 2)
  
  n1 <- as.numeric(chosen_nodes$parent)
  n2 <- as.numeric(chosen_nodes$node)
  dist_f <- d_n[as.character(n1)]
  dist_s <- d_n[as.character(n2)]
  
  # chosen_nodes2: ordered connections (for probability differences)
  chosen_nodes2[which(dist_f < dist_s), 1] <- n2[which(dist_f < dist_s)]
  chosen_nodes2[which(dist_f < dist_s), 2] <- n1[which(dist_f < dist_s)]
  
  chosen_nodes2[which(dist_f >= dist_s), 1] <- n1[which(dist_f >= dist_s)]
  chosen_nodes2[which(dist_f >= dist_s), 2] <- n2[which(dist_f >= dist_s)]
  
  score_norm <- t(mapply(function(ps){position_score(ps, msa, num_nodes, num_leaves, total_pos, human_plc,
                                                     node_human, nodes_raxml, human_leaf_len, dist_node, dist_leaf,
                                                     chosen_leaves, chosen_nodes2, d_n, d_l)}, rep(positions)))
  
  approaches <- c("_Blosum", "_Diversity", "_Frequency", "_WeightedFrequency",
                  "_Freq_WithoutGap", "_WeightedFreq_WithoutGap")
  
  print_wl <- lapply(1:length(approaches), function(p){
    score <- matrix(unlist(score_norm[ ,p]), nrow = length(positions) , ncol = 20, byrow = TRUE)
    score <- cbind(rep(positions), score)
    colnames(score) <- c("Pos/AA", num_to_aa(1:20))
    filename <- approaches[p]
    write.csv(score, sprintf("%s.csv", paste(output_name, filename, sep = "")), row.names = FALSE, quote = FALSE)
  })
  
}

position_score <- function(ps, msa, num_nodes, num_leaves, total_pos, human_plc, node_human, nodes_raxml, human_leaf_len, dist_node, dist_leaf, chosen_leaves, chosen_nodes2, d_n, d_l) {
  position <- ps
  
  position_vec <- msa[, ps]
  position_num <- aa_to_num(position_vec)
  
  ################## Weights per Node
  weights2 <- weight_fnc(d_n, d_l, human_plc, "Max05")
  weight_leaf2 <- weights2[1:num_leaves]
  weight_node2 <- tail(weights2,num_nodes)
  ####################
  
  # Blosum62
  logodd2prob <- function(A) {
    A <- exp(A)
    A <- A / (1 + A) 
    return(A)
  }
  amino_acids <- c("G", "A", "L", "M", "F", "W", "K", "Q", "E", "S", "P", "V", "I", "C", "Y", "H", "R", "N", "D", "T")
  data("BLOSUM62")
  blosum_mat <- BLOSUM62[amino_acids, amino_acids]
  Z_prob <- logodd2prob(blosum_mat)
  ref_aa <- position_num[human_plc]
  score_blosum <- Z_prob[ref_aa, ]
  
  # Diversity
  unq <- length(unique(position_num)) - 1*is.element(21,position_num)
  score_div <- matrix(unq, 1, 20)
  
  # Frequency, Frequency Without Gaps
  score <- matrix(0,1,20)
  for (i in 1:20){
    score[i] <- sum(position_num == i)
  }
  score_freq <- score/length(position_num)
  score_freq_wogap <- score/(length(position_num)-length(which(position_num>20)))
  
  # Weighted Frequency, Weighted Frequency Without Gaps
  score <- matrix(0,1,20)
  for (i in 1:num_leaves){
    wh_aa <- position_num[i]
    if (wh_aa<21){
      score[wh_aa] <- score[wh_aa] + weight_leaf2[i]*1
    }
  }
  score_weightedfreq <- score/length(position_num)
  score_weightedfreq_wogap <- score/(length(position_num)-length(which(position_num>20)))
  
  scores <- list()
  scores$blosum <- score_blosum
  scores$diversity <- score_div
  scores$frequency <- score_freq
  scores$weightedfreq <- score_weightedfreq
  scores$freqnogap <- score_freq_wogap
  scores$weightedfreqnogap <- score_weightedfreq_wogap
  
  return(scores)
}

weight_fnc <- function(d_n, d_l, human_plc, parameter) {
  if (parameter=="Max05"){
    d_l_n <- d_l[-human_plc]
    min_l <- min(d_l_n)
    d_l2 <- d_l/min_l
    d_n2 <- d_n/min_l
    
    d_l2 <- d_l2 +1
    d_n2 <- d_n2 +1
    weight_node <- 1/d_n2
    weight_leaf <- 1/d_l2
    
  }
  weights = c(weight_leaf, weight_node)
  
  return(weights)
}

csv_file <- compute_score(file_nwk=args[1],file_fasta=args[2], output_name=args[3],human_id=args[4],'all')
