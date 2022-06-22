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
  
  # Length of the leaf of human
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
  
  parameters <- unlist(str_split(parameters, pattern = ","))
  
  score_norm <- mapply(function(ps, parameter){position_score(ps, x, msa, num_nodes, num_leaves, total_pos, human_plc, node_human, nodes_raxml, human_leaf_len, dist_node, dist_leaf, parameter, leaves_conn, nodes_conn, chosen_leaves, chosen_nodes2, d_n, d_l)}, rep(positions, length(parameters)), rep(parameters, each = length(positions)))
  score_norm_with_leaf <- matrix(unlist(score_norm), nrow = length(positions) * length(parameters), ncol = 20, byrow = TRUE)
  score_norm_with_leaf <- cbind(rep(positions, length(parameters)), score_norm_with_leaf)

  colnames(score_norm_with_leaf) <- c("Pos/AA", num_to_aa(1:20))
  
  print_wl <- lapply(1:length(parameters), function(p){
    score_to_print_wl <- score_norm_with_leaf[(positions + length(positions)*(p - 1)), ]
    filename <- parameters[p]
    write.csv(score_to_print_wl, sprintf("%s/%s/%s.csv", save_path, uniprot_id, paste(output_name, "_weight_", filename, sep = "")), row.names = FALSE, quote = FALSE)
  })
  
}

position_score <- function(ps, x, msa, num_nodes, num_leaves, total_pos, human_plc, node_human, nodes_raxml, human_leaf_len, dist_node, dist_leaf, parameter, leaves_conn, nodes_conn, chosen_leaves, chosen_nodes2, d_n, d_l) {
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
  prob_leaves <- matrix(0, num_leaves, 20)
  prob_leaves[cbind(which(position_num <= 20), position_num[which(position_num <= 20)])] <- 1
  
  gaps <- which(position_num == 21)
  
  diff_leaves <- matrix(0, num_leaves, 20)
  diff_leaves <- prob_leaves - matrix_prob[(chosen_leaves$parent - num_leaves), ]
  diff_leaves[human_plc,]<- -diff_leaves[human_plc,]
  
  diff_nodes <- matrix(0, num_nodes-1, 20)
  diff_nodes <- matrix_prob[chosen_nodes2[,1] - num_leaves, ] - matrix_prob[chosen_nodes2[,2] - num_leaves, ]

  
  ################## Weight per node
  weights <- weight_fnc(d_n, d_l, human_plc, parameter, leaves_conn, nodes_conn, mxx) # Updated with related parameters
  weight_leaf <- weights[1:num_leaves]
  weight_node <- tail(weights,num_nodes)
  ####################
  
  score <- matrix(0,1,20)
  
  s1 <- sapply(1:20, function(ii){
    dif_pr <- diff_nodes[1:(num_nodes-1),ii]
    dif_pr[dif_pr<0] <- 0
    sel_node <- chosen_nodes2[1:length(dif_pr), 1] - num_leaves
    score[ii] <<- score[ii] + sum(weight_node[sel_node] * dif_pr)
  })
  
  aa_f <- position_num[human_plc]
  vect_human <- matrix_prob[node_human - num_leaves,]
  vect_human[aa_f]<-0
  
  score_without_leaf <- score
  score_without_leaf <- score_without_leaf + weight_node[(node_human-num_leaves)]*vect_human
  score <- score + weight_node[(node_human-num_leaves)]*vect_human
  
  s2 <- sapply(1:20, function(ii){
    diff_lf <- diff_leaves[1:num_leaves,ii]
    diff_lf[gaps] <-  0
    diff_lf[diff_lf<0] <- 0
    
    s1 <- sum(weight_leaf[((1:length(diff_lf)) != human_plc )] * diff_lf[((1:length(diff_lf)) != human_plc )])
    score[ii] <<- score[ii] + s1
  })
  
  aa_f <- position_num[human_plc]
  if (aa_f != 21){
    score[aa_f] <- score[aa_f] + weight_leaf[human_plc]*1
    score_without_leaf[aa_f] <- score_without_leaf[aa_f] + weight_leaf[human_plc]*1
  }
  
  score_norm <- score
  score_wol_norm <- score_without_leaf
  
  kk1 <- sum(score_norm)
  kk2 <- sum(score_wol_norm)
  
  sums1 <- kk1-max(score_norm)
  sums2 <- kk2-max(score_wol_norm)
  
  eps <- 10^(-15)
  sums1_upd <- (-(length(which(score_norm<0.0001))*0.1)/20+0.1)*(sums1)
  sums2_upd <- (-(length(which(score_wol_norm<0.0001))*0.1)/20+0.1)*(sums2)
  scores <- list()
  scores$score_with_leaf <- 1-log((score_norm*0.9 + sums2_upd)/(num_nodes+num_leaves)+eps)/log(eps)

  return(scores)
}

weight_fnc <- function(d_n, d_l, human_plc, parameter, leaves_conn, nodes_conn, mxx) {
  if (parameter=="Max05"){
    d_l_n <- d_l[-human_plc]
    min_l <- min(d_l_n)
    d_l2 <- d_l/min_l
    d_n2 <- d_n/min_l
    
    d_l2 <- d_l2 +1
    d_n2 <- d_n2 +1
    weight_node <- 1/d_n2
    weight_leaf <- 1/d_l2
 
  } else if (parameter=="Max05_Node"){
    d_l_n <- d_l[-human_plc]
    min_l <- min(d_n)
    d_l2 <- d_l/min_l
    d_n2 <- d_n/min_l
    
    d_l2 <- d_l2 +1
    d_n2 <- d_n2 +1
    weight_node <- 1/d_n2
    weight_leaf <- 1/d_l2
    
  } else if (parameter=="Mix_Gaussian_Max05"){
    d_l_n <- d_l[-human_plc]
    min_l <- min(d_n)
    d_l2 <- d_l/min_l
    d_n2 <- d_n/min_l
    
    d_l2 <- d_l2 +1
    d_n2 <- d_n2 +1
    weight_node1 <- 1/d_n2
    weight_leaf1 <- 1/d_l2
    
    param <- mean(c(d_n, d_l))
    weight_leaf2 <- exp(-d_l^2/param^2)/2
    weight_node2 <- exp(-d_n^2/param^2)/2
    
    weight_node<-sqrt(weight_node1*weight_node2)
    weight_leaf<-sqrt(weight_leaf1*weight_leaf2)
    
  } else if (parameter == "Gaussian_mean"){
    param <- mean(c(d_n, d_l))
    weight_leaf <- exp(-d_l^2/param^2)
    weight_node <- exp(-d_n^2/param^2)
    
  } else if (parameter == "Gaussian_median"){
    param <- median(c(d_n, d_l))
    weight_leaf <- exp(-d_l^2/param^2)
    weight_node <- exp(-d_n^2/param^2)

  } else if (parameter == "Mix_Gaussian_NodeCount"){
    param <- mean(c(d_n, d_l))
    weight_node = sqrt(exp(-d_n^2/param^2)*1/nodes_conn)
    weight_leaf = sqrt(exp(-d_l^2/param^2)*1/leaves_conn)
    
  }
  weights = c(weight_leaf, weight_node)
  return(weights)
}

csv_file <- compute_score(file_nwk=args[1],file_rst=args[2],file_fasta=args[3], output_name=args[4],human_id=args[5],'all', parameters = args[6])
