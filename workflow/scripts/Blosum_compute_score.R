#!/usr/bin/env Rscript
library(ape)
library(Biostrings)
library(ggplot2)
library(ggtree)
library(tidytree)
library(stringr)
library(dplyr)
#source("./position_score.R")
args = commandArgs(trailingOnly=TRUE)

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

compute_score <- function(file_nwk, file_rst, output_name, species_id, pos_chosen, parameters) {
  
  tr_org = read.tree(file_nwk)
  x <- readLines(file_rst)
  
  # ggtree(tr_org)
  tree_info <- as_tibble(tr_org)
  
  connections <- x[12]
  
  names_all <- tr_org[["tip.label"]]
  
  num_leaves <- length(tr_org[["tip.label"]])
  num_nodes <- tr_org[["Nnode"]]
  
  max_nodes <- x[17]
  num_nodes_codeml <- as.numeric(unlist(regmatches(max_nodes, gregexpr("[[:digit:]]+", max_nodes)))[2])
  
  if (num_nodes_codeml-num_leaves != num_nodes){
    print("Number of nodes is less than expected: CODEML")
  }
  
  dd_node <- dist.nodes(tr_org)
  dist_leaf <- dd_node[1:num_leaves, 1:num_leaves]
  dist_node <- dd_node[(num_leaves+1):(num_leaves + num_nodes), (num_leaves+1):(num_leaves + num_nodes)]
  connects <- str_split(connections, pattern = "[ ]", simplify = TRUE)
  connects <- connects[connects !=""]
  
  h_name <- species_id
  human_codeml <- names_all[sapply(X = names_all, function(n){ grepl(n, h_name, fixed = TRUE)})]
  
  leaf_human <- which(tree_info$label == human_codeml)
  
  index1 <- which(grepl(pattern = "Prob distribution at", x = x[30:length(x)], fixed = TRUE) == TRUE)[1] + 29
  
  seq_len_vect <- x[index1-2]
  total_pos <- str_split(seq_len_vect, pattern = "[ ]", simplify = TRUE)
  total_pos <- as.numeric(total_pos[total_pos !=""][1])
  if (pos_chosen[1] == "all"){
    positions <- 1:total_pos
    score_all <- matrix(0, total_pos, 21)
  } else {
    positions <- pos_chosen
    score_all <- matrix(0, length(positions), 21)
  }
  
  parameters <- unlist(str_split(parameters, pattern = ","))
  # parameters <- c(parameters, "mean", "median", 0)
  
  score_norm <- t(mapply(function(ps, parameter){position_score(ps, x, num_nodes, num_leaves, total_pos, leaf_human, tree_info, dist_node, dist_leaf, parameter)}, rep(positions, length(parameters)), rep(parameters, each = length(positions))))
  score_norm_with_leaf <- matrix(unlist(score_norm[ ,1]), nrow = length(positions) * length(parameters), ncol = 20, byrow = TRUE)
  score_norm_without_leaf <- matrix(unlist(score_norm[ ,2]), nrow = length(positions) * length(parameters), ncol = 20, byrow = TRUE)
  
  score_norm_with_leaf <- cbind(rep(positions, length(parameters)), score_norm_with_leaf)
  score_norm_without_leaf <- cbind(rep(positions, length(parameters)), score_norm_without_leaf)
  
  colnames(score_norm_with_leaf) <- c("Pos/AA", num_to_aa(1:20))
  colnames(score_norm_without_leaf) <- c("Pos/AA", num_to_aa(1:20))
  
  print_wl <- lapply(1:length(parameters), function(p){
    score_to_print_wl <- score_norm_with_leaf[(positions + length(positions)*(p - 1)), ]
    score_to_print_wol <- score_norm_without_leaf[(positions + length(positions)*(p - 1)), ]
    filename <- ifelse(parameters[p] == "0", "max05", parameters[p])
    filename <- ifelse(parameters[p] == "X", "max05_Gauss", parameters[p]) ### New Line (11.11)
    write.csv(score_to_print_wl, sprintf("%s.csv", paste(output_name, "_wl_param_", filename, sep = "")), row.names = FALSE, quote = FALSE)
    write.csv(score_to_print_wol, sprintf("%s.csv", paste(output_name, "_wol_param_", filename, sep = "")), row.names = FALSE, quote = FALSE)
    })
  
}

position_score <- function(ps, x, num_nodes, num_leaves, total_pos, leaf_human, tree_info, dist_node, dist_leaf, parameter) {
  position <- ps
  b1 <- matrix(0,1,num_nodes)
  b1[1] <- 29 + position -1 
  b1[2:length(b1)] <- (5 + total_pos) * (c(2:length(b1)) -1) + b1[1]
  
  TT <- x[b1]
  matrix_prob <- matrix(0, num_nodes, 20)
  
  splits_1 <- apply(as.matrix(TT), c(1), function(t){unlist(str_split(t, pattern = ":"))[1]})
  splits_2 <- apply(as.matrix(TT), c(1), function(t){unlist(str_split(t, pattern = ":"))[2]})
  
  probs <- matrix(as.numeric(unlist(regmatches(splits_2, gregexpr("(?<=\\().*?(?=\\))", splits_2, perl=T)))), nrow = num_nodes, 20, byrow = TRUE)
  pos_index <- unlist(gregexpr(pattern ="[(]",splits_2[1]))
  rr <- aa_to_num(str_sub(splits_2[1], pos_index - 1, pos_index -1))
  matrix_prob[,rr] <- probs
  
  ###########################################################################
  ###############      BLOSUM_BASED AA SUBST. - 30DEC (Lines 113:149) #######
  ###########################################################################
  logodd2prob <- function(A) {
    A <- exp(A)
    A <- A / (1 + A) 
    return(A)
  }
  logodd2dist <- function(A) {
    A <- exp(-A)
    A <- A / (1 + A) 
    return(A)
  }
  amino_acids <- c("G", "A", "L", "M", "F", "W", "K", "Q", "E", "S", "P", "V", "I", "C", "Y", "H", "R", "N", "D", "T")
  data("BLOSUM62")
  blosum_mat <- BLOSUM62[amino_acids, amino_acids]
  Z_dist <- logodd2dist(blosum_mat)
  Z_prob <- logodd2prob(blosum_mat)
  hc <- hclust(as.dist(Z_dist), method = "ward.D2")
  # hc <- hclust(as.dist(Z), method = "complete")
  # hc <- hclust(as.dist(Z), method = "single")
  aa_clusters <- cutree(hc, k = 4) #change number of clusters k if needed
  
  observed_aa <- (matrix_prob > 0) * 1
  colnames(observed_aa) <- amino_acids
  
  observed_aa_clusters <- sapply(1:nrow(matrix_prob), function(i){aa_clusters[names(which(observed_aa[i, ] == 1))]})
  added_aa <- lapply(1:nrow(matrix_prob), function(i){ names(aa_clusters[aa_clusters %in% unique(observed_aa_clusters[[i]])]) })
  
  matrix_prob_new <- matrix(0, nrow(matrix_prob), 20)
  colnames(matrix_prob_new) <- amino_acids
  s1 <- sapply(1:nrow(matrix_prob), function(i){
    Z_added <- Z_prob[, added_aa[[i]]]
    expected_prob <- Z_added / matrix(rowSums(Z_added), nrow(Z_added), ncol(Z_added), byrow = FALSE)
    added_prob <- matrix_prob[i, ] %*% expected_prob
    matrix_prob_new[i, colnames(added_prob)] <<- added_prob 
  })
  
  weight_freq <- sapply(1:nrow(matrix_prob), function(i){sqrt((1/length(unique(observed_aa_clusters[[i]]))) * max(matrix_prob[i, ])) })
  
  matrix_prob <- matrix_prob * matrix(weight_freq, nrow(matrix_prob), ncol(matrix_prob), byrow = FALSE) + 
    matrix_prob_new * (1 - matrix(weight_freq, nrow(matrix_prob), ncol(matrix_prob), byrow = FALSE))
  ##############################
  ##############################
  
  position_vec <- str_sub(splits_1[1], nchar(splits_1[1]) - num_leaves + 1, nchar(splits_1[1]))
  
  human_plc <- leaf_human
  human_leaf_len <- as.double(tree_info[human_plc, 3])
  human_node <- as.double(tree_info[human_plc, 1]) - num_leaves
  
  position_num <- aa_to_num(unlist(strsplit(position_vec, "")))
  
  prob_leaves <- matrix(0, num_leaves, 20)
  prob_leaves[cbind(which(position_num <= 20), position_num[which(position_num <= 20)])] <- 1
  
  gaps <- which(position_num == 21)
  
  chosen_leaves <- tree_info[1:num_leaves,1:2]
  chosen_nodes <- tree_info[num_leaves+2:num_nodes,1:2]
  
  d_n <- dist_node[human_node,1:num_nodes] + human_leaf_len
  d_l <- dist_leaf[human_plc,1:num_leaves]
  
  chosen_nodes2 <- matrix(0, num_nodes-1, 2)
  n1 <- as.numeric(chosen_nodes$parent)
  n11 <- n1-num_leaves
  n2 <- as.numeric(chosen_nodes$node)
  n22 <- n2-num_leaves
  dist_f <- d_n[n11] 
  dist_s <- d_n[n22] 
  
  chosen_nodes2[which(dist_f < dist_s), 1] <- n2[which(dist_f < dist_s)]
  chosen_nodes2[which(dist_f < dist_s), 2] <- n1[which(dist_f < dist_s)]
  
  chosen_nodes2[which(dist_f >= dist_s), 1] <- n1[which(dist_f >= dist_s)]
  chosen_nodes2[which(dist_f >= dist_s), 2] <- n2[which(dist_f >= dist_s)]
  
  ###########################################################################
  ###############      NEW PART - 15NOV                      ################
  ###########################################################################
  
  # Number of nodes between nodes & leaf of human
  nodes_conn <- numeric(num_nodes)
  nodes_conn[human_node-num_leaves] <- 1
  chs <- c()
  chs2 <- c()
  ##########################
  inds <- chosen_nodes2[chosen_nodes2[1:dim(chosen_nodes2)[1],2]==(human_node+num_leaves),1]
  nodes_conn[inds-num_leaves] <- 2
  chs <- inds
  
  s0 <- sapply(3:num_leaves, function(i){
    for (j in chs){
      inds <- chosen_nodes2[chosen_nodes2[1:dim(chosen_nodes2)[1],2]==j,1]
      if (length(inds)!=0){
        nodes_conn[inds-num_leaves] <<- i
        chs2 <- c(chs2, inds)
      }
    }
    chs <<- chs2
    chs2 <- c()
  })
  
  # Number of nodes between leaves & leaf of human
  leaves_conn <- nodes_conn[pull(chosen_leaves[1:dim(chosen_leaves)[1],1])-num_leaves] 
  
  diff_leaves <- matrix(0, num_leaves, 20)
  diff_leaves <- prob_leaves - matrix_prob[as.numeric(chosen_leaves$parent) - num_leaves, ]
  
  diff_nodes <- matrix(0, num_nodes-1, 20)
  diff_nodes <- matrix_prob[chosen_nodes2[,1] - num_leaves, ] - matrix_prob[chosen_nodes2[,2] - num_leaves, ]
  
  ################## weights
  weights <- weight_fnc(d_n, d_l, human_plc, parameter, leaves_conn, nodes_conn) # Updated with related parameters
  weight_leaf <- weights[1:num_leaves]
  weight_node <- tail(weights,num_nodes)
  ####################
  
  score <- matrix(0,1,20)
  
  
  
  s1 <- sapply(1:20, function(ii){
    dif_pr <- diff_nodes[1:num_nodes-1,ii]
    v1 <- matrix_prob[1:num_nodes,ii]
    if (dif_pr[1]>0) {
      sel_node <- chosen_nodes2[1, 1] - num_leaves
      score[ii] <- weight_node[sel_node] * dif_pr[1]
    }
    
    sel_node <- chosen_nodes2[2:length(dif_pr), 1] - num_leaves
    score[ii] <<- score[ii] + sum(weight_node[sel_node] * dif_pr[2:length(dif_pr)] * (dif_pr[2:length(dif_pr)] > 0))
  })
  
  score_without_leaf <- score
  
  s2 <- sapply(1:20, function(ii){ 
    diff_lf <- diff_leaves[1:num_leaves,ii]
    diff_lf[gaps] <-  0
    
    s1 <- sum(weight_leaf[((1:length(diff_lf)) != human_plc ) & diff_lf > 0] * diff_lf[((1:length(diff_lf)) != human_plc ) & diff_lf > 0])
    score[ii] <<- score[ii] + s1
  })
  
  aa_f <- position_num[human_plc]
  if (aa_f != 21){
    score[aa_f] <- score[aa_f] + weight_leaf[human_plc]*1
    score_without_leaf[aa_f] <- score_without_leaf[aa_f] + weight_leaf[human_plc]*1
  }
  
  
  score_norm <- score*(1/sum(score))
  score_wol_norm <- score_without_leaf*(1/sum(score_without_leaf))
  scores <- list()
  scores$score_with_leaf <- score_norm
  scores$score_without_leaf <- score_wol_norm
  return(scores)
}

weight_fnc <- function(d_n, d_l, human_plc, parameter, leaves_conn, nodes_conn) {
  if (parameter=="0"){
    d_l_n <- d_l[-human_plc]
    min_l <- min(d_l_n)
    d_l2 <- d_l/min_l
    d_n2 <- d_n/min_l
    
    d_l2 <- d_l2 +1
    d_n2 <- d_n2 +1
    weight_node <- 1/d_n2
    weight_leaf <- 1/d_l2
  } else if (parameter == "mean"){
    param <- mean(c(d_n, d_l))
    weight_leaf <- exp(-d_l^2/param^2)
    weight_node <- exp(-d_n^2/param^2)
  } else if (parameter == "median"){
    param <- median(c(d_n, d_l))
    weight_leaf <- exp(-d_l^2/param^2)
    weight_node <- exp(-d_n^2/param^2)
  } else if (parameter == "X"){
    d_l_n <- d_l[-human_plc]
    min_l <- min(d_l_n)
    min_n <- min(d_n)
    min_ch <- min(min_l, min_n)
    d_l2 <- d_l-min_ch
    d_n2 <- d_n-min_ch
    weight_leaf <- exp(-d_l^2)/2
    weight_node <- exp(-d_n^2)/2
    weight_leaf[human_plc] <- 1
  } else if (parameter == "CountNodes_1"){      ### NewFunction 1 (15Nov)
    weight_node = (exp(-d_n^2) + 1/nodes_conn)/2
    weight_leaf = (exp(-d_l^2) + 1/leaves_conn)/2 
  } else if (parameter == "CountNodes_2"){      ### NewFunction 2 (15Nov)
    weight_node = (exp(-d_n^2) + exp(-nodes_conn^2))/2  
    weight_leaf = (exp(-d_l^2) + exp(-leaves_conn^2))/2  
  } else if (parameter == "CountNodes_3"){      ### NewFunction 3 (15Nov)
    weight_node = sqrt(exp(-d_n^2)*1/nodes_conn)  
    weight_leaf = sqrt(exp(-d_l^2)*1/leaves_conn)
  } else if (parameter == "CountNodes_4"){      ### NewFunction 4 (15Nov)
    weight_node = exp(-(sqrt(d_n*nodes_conn))^2) 
    weight_leaf = exp(-(sqrt(d_l*leaves_conn))^2)
  } else {
    param <- as.double(parameter)
    weight_leaf <- exp(-d_l^2/param^2)
    weight_node <- exp(-d_n^2/param^2)
  }
  weights = c(weight_leaf, weight_node)
  return(weights)
}

csv_file <- compute_score(file_nwk=args[1],file_rst=args[2],output_name=args[3],species_id=args[4],'all', parameters = args[5])