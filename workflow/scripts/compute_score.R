#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
aa_to_num <- function(aa) {
  if (aa == "G") {
    num <- 1
  } else if (aa == "A"){
    num <- 2
  } else if (aa == "L"){
    num <- 3
  } else if (aa == "M"){
    num <- 4
  } else if (aa == "F"){
    num <- 5
  } else if (aa == "W"){
    num <- 6
  } else if (aa == "K"){
    num <- 7
  } else if (aa == "Q"){
    num <- 8
  } else if (aa == "E"){
    num <- 9
  } else if (aa == "S"){
    num <- 10
  } else if (aa == "P"){
    num <- 11
  } else if (aa == "V"){
    num <- 12
  } else if (aa == "I"){
    num <- 13
  } else if (aa == "C"){
    num <- 14
  } else if (aa == "Y"){
    num <- 15
  } else if (aa == "H"){
    num <- 16
  } else if (aa == "R"){
    num <- 17
  } else if (aa == "N"){
    num <- 18
  } else if (aa == "D"){
    num <- 19
  } else if (aa == "T"){
    num <- 20
  } else {
    num <- 21
  }
  return(num)
}
num_to_aa <- function(num) {
  if (num == 1) {
    aa <- "G"
  } else if (num == 2){
    aa <- "A"
  } else if (num == 3){
    aa <- "L"
  } else if (num == 4){
    aa <- "M"
  } else if (num == 5){
    aa <- "F"
  } else if (num == 6){
    aa <- "W"
  } else if (num == 7){
    aa <- "K"
  } else if (num == 8){
    aa <- "Q"
  } else if (num == 9){
    aa <- "E"
  } else if (num == 10){
    aa <- "S"
  } else if (num == 11){
    aa <- "P"
  } else if (num == 12){
    aa <- "V"
  } else if (num == 13){
    aa <- "I"
  } else if (num == 14){
    aa <- "C"
  } else if (num == 15){
    aa <- "Y"
  } else if (num == 16){
    aa <- "H"
  } else if (num == 17){
    aa <- "R"
  } else if (num == 18){
    aa <- "N"
  } else if (num == 19){
    aa <- "D"
  } else if (num == 20){
    aa <- "T"
  } else {
    aa <- 21
  }
  return(aa)
}
compute_score <- function(file_nwk, file_rst, output_name, human_id, pos_chosen) {
  
  library("ape")
  library("Biostrings")
  library("ggplot2")
  library("ggtree")
  library(tidytree)
  
  tr_org = read.tree(file_nwk)
  x <- readLines(file_rst)
  
  
  ggtree(tr_org)
  tree_info <- as_tibble(tr_org)
  
  connections <- x[12]
  
  names_all = tr_org[["tip.label"]]
  
  num_leaves = length(tr_org[["tip.label"]])
  num_nodes = tr_org[["Nnode"]]
  
  max_nodes <- x[17]
  num_nodes_codeml <- 0
  for (i in 1:nchar(max_nodes)){
    if (substr(max_nodes,i,i+1) == "to"){
      for (j in i+4:nchar(max_nodes)){
        if (substr(max_nodes,j,j+2) == "are"){
          num_nodes_codeml <- as.double(substr(max_nodes,i+3,j-2))
          break
        }
      }
    }
    if (num_nodes_codeml != 0){
      break
    }
  }  
  
  if (num_nodes_codeml-num_leaves != num_nodes){
    print("Number of nodes is less than expected: CODEML")
  }
  
  
  dd_node = dist.nodes(tr_org)
  dist_leaf <- dd_node[1:num_leaves,1:num_leaves]
  dist_node <- dd_node[num_leaves+1:num_nodes,num_leaves+1:num_nodes]
    
  connects <- 0
  i <- 3
  for (val in 1:28)
  {
    s <- substr(connections, i, i+5)
    i <- i+9
    connects[val] <- s
  }
  
  h_name <- human_id
  
  for (i in 1:num_leaves)
  {
    vect <- names_all[i]
    if (nchar(vect)>=nchar(h_name)){
      if (substr(vect,1,nchar(h_name)) == h_name){
        human_codeml = vect;
        break
      } else if (substr(vect,nchar(vect)-nchar(h_name)+1,nchar(vect)) == h_name) {
        human_codeml = vect;
        break
      }
    }
  }
  
  for (i in 1:num_leaves){
    if (tree_info[i,4] == human_codeml) {
      leaf_human <- i
      break
    }
  }
    
  vect2 <- 1
  
  for (i in 30:length(x)){
    vect <- x[i]
    if (substr(vect, 1, 20) == "Prob distribution at") {
      seq_len_vect <- x[i-2]
      ind <- nchar(seq_len_vect)-2
      for (j in 1:ind){
        if (substr(seq_len_vect,j,j) != " " & substr(seq_len_vect,j+1,j+1) == " "){
          vect2 <- as.double(substr(seq_len_vect,1,j))
          break
        }
      }
    }
    if (vect2 > 1){
      break
    }
  }
  
  total_pos <- vect2
  if (pos_chosen[1] == "all"){
    positions <- 1:total_pos
    score_all <- matrix(0, vect2,21)
  } else {
    positions <- pos_chosen
    score_all <- matrix(0, length(positions),21)
  }
  
  for (ps in 1:length(positions)){
    
    position <- positions[ps]
    b1 <- matrix(0,1,num_nodes)
    
    for (i in 1:num_nodes){
      if (i==1){
        b1[1] <- 29 + position -1
      } else {
        b1[i] <- 5 + total_pos + b1[i-1]
      }
    }
    TT <- x[b1]
    
    matrix_prob <- matrix(0, num_nodes, 20)
    ent <- 0
    for (num_i in 1:num_nodes){
      T_vec <- TT[num_i]
      for (num_j in 1:nchar(T_vec)){
        if (substr(T_vec,num_j,num_j)=="("){
          ent <- ent + 1
          prb <- substr(T_vec,num_j+1,num_j+5)
          aa <- substr(T_vec,num_j-1,num_j-1)
          if (aa !="-"){
            rr <- aa_to_num(aa)
            matrix_prob[num_i,rr] <- as.double(prb)  ### CORRECT!!!
          }
          if (ent == 1){
            position_vec <- substr(T_vec, num_j-4-num_leaves+1,num_j-4)
          }
        }
      }
    }
    
    human_plc <- leaf_human
    human_leaf_len <- as.double(tree_info[human_plc, 3])
    human_node <- as.double(tree_info[human_plc, 1]) - num_leaves
    
    position_num <- matrix(0, 1, num_leaves)
    
    for (i in 1:nchar(position_vec)){
      position_num[i] <- aa_to_num(substr(position_vec,i,i))
    }
    
    prob_leaves <- matrix(0, num_leaves, 20)
    for (i in 1:length(position_num)){
      if (position_num[i]<=20){
        prob_leaves[i,position_num[i]] <- 1
      }
    }
    
    gaps <- which(position_num %in% 21)
    
    chosen_leaves <- tree_info[1:num_leaves,1:2]
    chosen_nodes <- tree_info[num_leaves+2:num_nodes,1:2]
  
    d_n <- dist_node[human_node,1:num_nodes] + human_leaf_len
    d_l <- dist_leaf[4,1:num_leaves]
    
    chosen_nodes2 = matrix(0, num_nodes-1, 2)
    dim_cn <- dim(chosen_nodes)
    for (i in 1:dim_cn[1]){
      n1 <- as.double(chosen_nodes[i,1])
      n11 <- n1-num_leaves
      n2 <- as.double(chosen_nodes[i,2])
      n22 <- n2-num_leaves
      dist_f <- d_n[n11] 
      dist_s <- d_n[n22] 
      
      if (dist_f < dist_s){
        chosen_nodes2[i,1] <- n2
        chosen_nodes2[i,2] <- n1
      } 
      else {
        chosen_nodes2[i,1] <- n1
        chosen_nodes2[i,2] <- n2
      }
      
    }
    
    diff_leaves <- matrix(0, num_leaves, 20)
    for (i in 1:num_leaves){
      c_n <- as.double(chosen_leaves[i,1]) - num_leaves
      val <- matrix_prob[c_n,1:20]
      for (j in 1:20){
        diff_leaves[i,j] <- prob_leaves[i,j]-val[j]
      }
    }
  
    diff_nodes <- matrix(0, num_nodes-1, 20)
    for (i in 1:dim_cn[1]){
      c_n1 <- chosen_nodes2[i,1]-num_leaves
      c_n2 <- chosen_nodes2[i,2]-num_leaves
      val1 <- matrix_prob[c_n1,1:20]
      val2 <- matrix_prob[c_n2,1:20]
      for (j in 1:20){
        diff_nodes[i,j] <- val1[j]-val2[j]
      }
    }
  
    d_l_n <- d_l[-4]
    min_l <- min(d_l_n)
    d_l2 <- d_l/min_l
    d_n2 <- d_n/min_l
    
    d_l2 <- d_l2 +1
    d_n2 <- d_n2 +1
    weight_node <- 1/d_n2
    weight_leaf <- 1/d_l2
    
    score <- matrix(0,1,20)
    score2 <- matrix(0,1,20)
    
    for (ii in 1:20){
      dif_pr <- diff_nodes[1:num_nodes-1,ii]
      v1 <- matrix_prob[1:num_nodes,ii]
      
      for (i in 1:length(dif_pr)){
        if (i==1 & dif_pr[1]>0){
          selected_pair <- chosen_nodes2[i,1:2]
          sel_node <- selected_pair[1]-num_leaves
          score[ii] <- score[ii] +  weight_node[sel_node]*dif_pr[1]
          score2[ii] <- score2[ii] + weight_node[sel_node]*v1[2]
        }
        else if (i>1 && dif_pr[i] > 0){
          selected_pair <- chosen_nodes2[i,1:2]
          sel_node = selected_pair[1]-num_leaves
          score[ii] = score[ii] + weight_node[sel_node]*dif_pr[i]
          score2[ii] = score2[ii] + weight_node[sel_node]*v1[sel_node]
        }
        
      }
      
    }
    
    for (ii in 1:20){
      dif_lf <- diff_leaves[1:num_leaves,ii]
      dif_lf[gaps] = 0
      
      for (s in 1:length(dif_lf)){
        if (s != human_plc & dif_lf[s]>0){
          score[ii] <- score[ii] + weight_leaf[s]*dif_lf[s]
          score2[ii] <- score2[ii] + weight_leaf[s]*dif_lf[s]
        }
        else if (s == human_plc & dif_lf[s]<0){
          score[ii] = score[ii] - weight_leaf[s]*dif_lf[s]
          score2[ii] = score2[ii] - weight_leaf[s]*dif_lf[s]
        }
      }
      
    }
  
    aa_f <- position_num[human_plc]
    if (aa_f != 21){
      score[aa_f] <- score[aa_f] + weight_leaf[human_plc]*1
    }
    
    
    score3 <- score*(1/sum(score))
    
    if (length(positions)>1){
      score_all[ps,1] <- ps
      score_all[ps,2:21] <- score3[1:20]
    } else {
      score_all[1] <- ps
      score_all[2:21] <- score3[1:20]
    }
    
  }
  
  header_aa <- vector()
  header_aa[1] <- "Pos/AA"
  for (i in 1:20){
    header_aa[i+1]<- num_to_aa(i)
  }
  
  total_length = 1+length(positions)
  score_to_print = matrix(0,total_length,21)
  score_to_print[1,1:21] = header_aa
  score_to_print[2:total_length,1:21] = score_all
  
  write.csv(score_to_print, paste(output_name,".csv",sep=""), row.names = FALSE, col.names = FALSE, quote = FALSE)
}

csv_file<-compute_score(file_nwk=args[1],file_rst=args[2],output_name=args[3],human_id=args[4],'all')

