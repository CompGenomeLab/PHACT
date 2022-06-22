library(ape)
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

compute_score <- function(file_nwk, file_rst, file_fasta, output_name, human_id, pos_chosen, parameters) {
  
  # Read tree file
  tr_org <- read.tree(file_nwk)
  x <- read.table(file = file_rst, sep = '\t', header = TRUE, fill = TRUE)
  colnames(x)[4:ncol(x)] <- gsub("p_", replacement = "", x = colnames(x)[4:ncol(x)], fixed = TRUE )
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
  total_pos <- max(x$Site)
  
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
  
  # Number of nodes between nodes & leaf of human
  nodes_conn <- numeric(num_nodes)
  nodes_conn[node_human-num_leaves] <- 1
  names(nodes_conn) <- names(d_n)
  chs <- c()
  chs2 <- c()

  inds <- chosen_nodes2[chosen_nodes2[,2]==(node_human),1]
  nodes_conn[as.character(inds)] <- 2
  chs <- inds
  
  s0 <- sapply(3:num_leaves, function(i){
    for (j in chs){
      inds <- chosen_nodes2[chosen_nodes2[,2]==j,1]
      if (length(inds)!=0){
        nodes_conn[as.character(inds)] <<- i
        chs2 <- c(chs2, inds)
      }
    }
    chs <<- chs2
    chs2 <- c()
  })
  
  # Number of nodes between leaves & leaf of human
  leaves_conn <- nodes_conn[as.character(chosen_leaves[,1])]
  
  #########################################
  
  score_norm <- t(mapply(function(ps, parameter){position_score(ps, x, msa, num_nodes, num_leaves, total_pos, human_plc, node_human, nodes_raxml, human_leaf_len, dist_node, dist_leaf, leaves_conn, nodes_conn, chosen_leaves, chosen_nodes2, d_n, d_l)}, rep(positions)))
  approaches <- c("_PHACT", "_WithoutLeaves", "_EqualWeight", "_AllDifferences", "_AllProbability_Equal", "_AllProbability", "_OnlyScore",
                  "_OnlyDiversity", "_EqualImportance", "_HighDiversity")
  
  print_wl <- lapply(1:length(approaches), function(p){
    score <- matrix(unlist(score_norm[ ,p]), nrow = length(positions) , ncol = 20, byrow = TRUE)
    score <- cbind(rep(positions), score)
    colnames(score) <- c("Pos/AA", num_to_aa(1:20))
    filename <- approaches[p]
    write.csv(score, sprintf("%s.csv", paste(output_name, filename, sep = "")), row.names = FALSE, quote = FALSE)
  })
  
}

position_score <- function(ps, x, msa, num_nodes, num_leaves, total_pos, human_plc, node_human, nodes_raxml, human_leaf_len, dist_node, dist_leaf, leaves_conn, nodes_conn, chosen_leaves, chosen_nodes2, d_n, d_l) {
  position <- ps
  
  b1 <- position + total_pos*(0:(num_nodes-1))
  TT <- x[b1,]
  matrix_prob <- matrix(0, num_nodes, 20)
  
  probs <- data.matrix((TT[, (4:ncol(TT))]))
  rownames(probs) <- NULL
  rr <- aa_to_num(colnames(x)[4:ncol(TT)])
  matrix_prob[,rr] <- probs
  matrix_prob <- matrix_prob[nodes_raxml,]
  
  position_vec <- msa[, ps]
  
  position_num <- aa_to_num(position_vec)
  gaps <- which(position_num == 21)
  
  prob_leaves <- matrix(0, num_leaves, 20)
  prob_leaves[cbind(which(position_num <= 20), position_num[which(position_num <= 20)])] <- 1
  diff_leaves <- matrix(0, num_leaves, 20)
  diff_leaves <- prob_leaves - matrix_prob[(chosen_leaves$parent - num_leaves), ]
  diff_leaves[human_plc,]<- -diff_leaves[human_plc,]
  
  diff_leaves_all <- prob_leaves
  diff_nodes <- matrix(0, num_nodes-1, 20)
  diff_nodes <- matrix_prob[chosen_nodes2[,1] - num_leaves, ] - matrix_prob[chosen_nodes2[,2] - num_leaves, ]
  diff_nodes_all <- matrix_prob[chosen_nodes2[,1] - num_leaves, ]
  
  ################## weights
  weights <- weight_fnc(d_n, d_l, human_plc, "Max05", leaves_conn, nodes_conn)
  weight_leaf <- weights[1:num_leaves]
  weight_node <- tail(weights,num_nodes)
  
  weights_eq <- weight_fnc(d_n, d_l, human_plc, "Equal", leaves_conn, nodes_conn)
  weight_leaf_eq <- weights_eq[1:num_leaves]
  weight_node_eq <- tail(weights_eq,num_nodes)
  ####################
  
  score <- matrix(0,1,20)
  score_eq <- matrix(0,1,20)
  score_alldiff <- matrix(0,1,20)
  score_allprob <- matrix(0,1,20)
  score_allprob_eq <- matrix(0,1,20)
  
  s1 <- sapply(1:20, function(ii){
    dif_pr <- diff_nodes[1:(num_nodes-1),ii]
    dif_pr[dif_pr<0] <- 0
    sel_node <- chosen_nodes2[1:length(dif_pr), 1] - num_leaves
    
    score[ii] <<- score[ii] + sum(weight_node[sel_node] * dif_pr)
    score_eq[ii] <<- score_eq[ii] + sum(weight_node_eq[sel_node] * dif_pr)
    
    dif_pr2 <- diff_nodes[1:(num_nodes-1),ii]
    score_alldiff[ii] <<- score_alldiff[ii] + sum(weight_node[sel_node] * dif_pr2)
    
    dif_pr3 <- diff_nodes_all[1:(num_nodes-1),ii]
    score_allprob[ii] <<- score_allprob[ii] + sum(weight_node[sel_node] * dif_pr3)
    score_allprob_eq[ii] <<- score_allprob_eq[ii] + sum(weight_node_eq[sel_node] * dif_pr3)
    
  })
  
  aa_f <- position_num[human_plc]
  vect_human <- matrix_prob[node_human - num_leaves,]
  vect_human[aa_f]<-0
  score_wol <- score + weight_node[(node_human-num_leaves)]*vect_human
  score <- score + weight_node[(node_human-num_leaves)]*vect_human
  score_wol_eq <- score_eq + weight_node_eq[(node_human-num_leaves)]*vect_human
  score_eq <- score_eq + weight_node_eq[(node_human-num_leaves)]*vect_human
  
  vect_human <- matrix_prob[node_human - num_leaves,]
  vect_human[aa_f]<-1-vect_human[aa_f]
  score_wol_alldiff <- score_alldiff + weight_node[(node_human-num_leaves)]*vect_human
  score_alldiff <- score_alldiff + weight_node[(node_human-num_leaves)]*vect_human
  
  vect_human <- matrix_prob[node_human - num_leaves,]
  score_wol_allprob <- score_allprob + weight_node[(node_human-num_leaves)]*vect_human
  score_allprob <- score_allprob + weight_node[(node_human-num_leaves)]*vect_human
  score_wol_allprob_eq <- score_allprob_eq + weight_node_eq[(node_human-num_leaves)]*vect_human
  score_allprob_eq <- score_allprob_eq + weight_node_eq[(node_human-num_leaves)]*vect_human
  
  
  s2 <- sapply(1:20, function(ii){
    diff_lf <- diff_leaves[1:num_leaves,ii]
    diff_lf[gaps] <-  0
    diff_lf[diff_lf<0] <- 0
    
    s1 <- sum(weight_leaf[((1:length(diff_lf)) != human_plc )] * diff_lf[((1:length(diff_lf)) != human_plc )])
    score[ii] <<- score[ii] + s1
    s1_eq <- sum(weight_leaf_eq[((1:length(diff_lf)) != human_plc )] * diff_lf[((1:length(diff_lf)) != human_plc )])
    score_eq[ii] <<- score_eq[ii] + s1_eq
    
    diff_lf <- diff_leaves[1:num_leaves,ii]
    diff_lf[gaps] <-  0
    
    s1_ad <- sum(weight_leaf[((1:length(diff_lf)) != human_plc )] * diff_lf[((1:length(diff_lf)) != human_plc )])
    score_alldiff[ii] <<- score_alldiff[ii] + s1_ad
    
    diff_lf_ap <- diff_leaves_all[1:num_leaves,ii]
    diff_lf_ap[gaps] <-  0
    
    s1 <- sum(weight_leaf[((1:length(diff_lf)) != human_plc )] * diff_lf_ap[((1:length(diff_lf)) != human_plc )])
    score[ii] <<- score[ii] + s1
    s1_eq <- sum(weight_leaf_eq[((1:length(diff_lf_ap)) != human_plc )] * diff_lf_ap[((1:length(diff_lf)) != human_plc )])
    score_eq[ii] <<- score_eq[ii] + s1_eq
    
  })
  
  aa_f <- position_num[human_plc]
  if (aa_f != 21){
    score[aa_f] <- score[aa_f] + weight_leaf[human_plc]*1
    score_wol[aa_f] <- score_wol[aa_f] + weight_leaf[human_plc]*1
    
    score_eq[aa_f] <- score_eq[aa_f] + weight_leaf_eq[human_plc]*1
    score_wol_eq[aa_f] <- score_wol_eq[aa_f] + weight_leaf_eq[human_plc]*1
    
    score_alldiff[aa_f] <- score_alldiff[aa_f] + weight_leaf[human_plc]*1
    score_wol_alldiff[aa_f] <- score_wol_alldiff[aa_f] + weight_leaf[human_plc]*1
    
    score_allprob[aa_f] <- score_allprob[aa_f] + weight_leaf[human_plc]*1
    score_wol_allprob[aa_f] <- score_wol_allprob[aa_f] + weight_leaf[human_plc]*1
    
    score_allprob_eq[aa_f] <- score_eq[aa_f] + weight_leaf_eq[human_plc]*1
    score_wol_allprob_eq[aa_f] <- score_wol_allprob_eq[aa_f] + weight_leaf_eq[human_plc]*1
  }
  
  sums <- sum(score_wol)-max(score_wol)
  sums_eq <- sum(score_wol_eq)-max(score_wol_eq)
  
  sums_alldiff <- sum(score_wol_alldiff)-max(score_wol_alldiff)
  sums_allprob_eq <- sum(score_wol_allprob_eq)-max(score_wol_allprob_eq)
  sums_allprob <- sum(score_wol_allprob)-max(score_wol_allprob)
  
  eps <- 10^(-15)
  sums_upd <- (-(length(which(score_wol<0.0001))*0.1)/20+0.1)*(sums)
  sums_upd_eq <- (-(length(which(score_wol_eq<0.0001))*0.1)/20+0.1)*(sums_eq)
  
  sums_upd_alldiff <- (-(length(which(score_alldiff<0.0001))*0.1)/20+0.1)*(sums_alldiff)
  sums_upd_allprob_eq <- (-(length(which(score_allprob_eq<0.0001))*0.1)/20+0.1)*(sums_allprob_eq)
  sums_upd_allprob <- (-(length(which(score_allprob<0.0001))*0.1)/20+0.1)*(sums_allprob)
  
  scores <- list()
  scores$PHACT <- 1-log((score*0.9 + sums_upd)/(num_nodes+num_leaves)+eps)/log(eps)
  scores$score_wol <- 1-log((score_wol*0.9 + sums_upd)/(num_nodes)+eps)/log(eps)
  scores$equalweight <- 1-log((score_eq*0.9 + sums_upd_eq)/(num_nodes+num_leaves)+eps)/log(eps)
  
  # Since the negative values are possible, we do not apply the scaling that we apply to
  # other PHACT-related approaches
  scores$alldiff <- (score_alldiff*0.9 + sums_upd_alldiff)/(num_nodes+num_leaves)
  
  scores$allprob_eq <- 1-log((score_allprob_eq*0.9 + sums_upd_allprob_eq)/(num_nodes+num_leaves)+eps)/log(eps)
  scores$allprob_max05 <- 1-log((score_allprob*0.9 + sums_upd_allprob)/(num_nodes+num_leaves)+eps)/log(eps)
  
  scores$onlyscore <- 1-log(score*1/(num_nodes+num_leaves)+eps)/log(eps)
  scores$onlydiv <- matrix(1-log(sums_upd/(num_nodes+num_leaves)+eps)/log(eps),1,20)
  sums_upd_eq <- (-(length(which(score_wol<0.0001))*0.5)/20+0.5)*(sums)
  scores$equalimport <- 1-log((score*0.5 + sums_upd_eq)/(num_nodes+num_leaves)+eps)/log(eps)
  sums_upd_div <- (-(length(which(score_wol<0.0001))*0.5)/20+0.5)*(sums)
  scores$higherdiv <- 1-log((score*0.1 +sums_upd_div)/(num_nodes+num_leaves)+eps)/log(eps)
  return(scores)
}

weight_fnc <- function(d_n, d_l, human_plc, parameter, leaves_conn, nodes_conn) {
  if (parameter=="Max05"){
    d_l_n <- d_l[-human_plc]
    min_l <- min(d_l_n)
    d_l2 <- d_l/min_l
    d_n2 <- d_n/min_l
    
    d_l2 <- d_l2 +1
    d_n2 <- d_n2 +1
    weight_node <- 1/d_n2
    weight_leaf <- 1/d_l2
 
  } else if (parameter=="Equal"){
    weight_node <- matrix(1,1,length(d_n))
    weight_leaf <- matrix(1,1,length(d_l))

  }
  weights = c(weight_leaf, weight_node)
  return(weights)
}

csv_file <- compute_score(file_nwk=args[1],file_rst=args[2],file_fasta=args[3], output_name=args[4],human_id=args[5],'all')
