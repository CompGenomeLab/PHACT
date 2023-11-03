#!/usr/bin/env Rscript
library(bio3d)
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

compute_score <- function(file_fasta, output_name, folder_name) {
  # Read fasta file, MSA
  fasta <- read.fasta(file = file_fasta)
  msa <- fasta$ali

  # Total number of positions from ancestralProbs file
  total_pos <- length(msa[1, ])
  
  names_all <- row.names(msa)
  
  h_name <- output_name
  human_codeml <- names_all[grep(pattern = h_name, x = names_all, fixed = TRUE)]
  
  # Chosen positions (all or some)
  positions <- 1:total_pos
  
  num_leaves <- length(msa[,1])
  msa_upd <- msa
  
  trim_msa <- msa_masking(msa, positions)
  trim_msa1 <- trim_msa$gaps
  trim_msa2 <- trim_msa$incomp
  trim_msa3 <- rownames(msa)[intersect(which(msa[,1]!="M"),which(msa[,1]!="-"))]
  trim_msa4 <- trim_msa$blocks
  ####################################################
  for (i in positions){
    vec <- c(trim_msa2[which(trim_msa2[,i]!=0),i], trim_msa4[which(trim_msa4[,i]!=0),i])
    if (i==1){
      vec <- c(vec, trim_msa3)
    }
    if (is.element(human_codeml, vec)){
      vec <- vec[-which(vec==human_codeml)]
    }
    if (length(vec)>0){
      vec <- unique(vec)
      msa_upd[t(mapply(function(jj){which(rownames(msa_upd)==vec[jj])}, rep(1:length(vec)))),i] <- "-"
    }
  }

  vect <- matrix(0, (2*num_leaves), 1)
  names <- rownames(msa_upd)
  k <- 1
  for (i in 1:length(names)){
    vect[k, ] <- paste(">", names[i], sep = "")
    vect[(k+1), ] <- paste(msa_upd[i,],collapse="")
    k <- k + 2
  }
  
  write.table(vect, quote = F, col.names = F, row.names = F, folder_name)
  
}

msa_masking <- function(msa, positions) {
  trim_msa_gaps <- matrix(0, length(msa[,1]), length(msa[1,]))
  trim_msa_incomp <- matrix(0, length(msa[,1]), length(msa[1,]))
  trim_msa_blocks <- matrix(0, length(msa[,1]), length(msa[1,]))
  
  total_pos <- length(msa[1,])
  num_leaves <- length(msa[,1])
  
  thr_gap <- 5
  cand <- matrix(0, num_leaves, total_pos)
  cand2 <- matrix(0, num_leaves, total_pos)
  
  leng <- 10
  
  # ONLY GAPS
  for (i in positions) {
    names <- rownames(msa)
    trs <- which(msa[,i]=="-") 
    if (length(trs)>0) {
      trim_msa_gaps[1:length(trs),i] <- names[trs]
    }
    
    # INCOMP CODE
    ind_inc <- 1
    pos_inc <- aa_to_num(msa[,i])
    if (i>5 && i<=(total_pos-5)){
      xx_inc <- setdiff(c((i-5):(i+5)),i)
      pos_bef_inc <- matrix(aa_to_num(msa[,xx_inc]),num_leaves,length(xx_inc))
    } else if (i<=5) {
      xx_inc <- setdiff(c(1:(i+5)),i)
      pos_bef_inc <- matrix(aa_to_num(msa[,xx_inc]),num_leaves,length(xx_inc))
    } else if (i>=(total_pos-5)) {
      xx_inc <- setdiff(c((i-5):total_pos),i)
      pos_bef_inc <- matrix(aa_to_num(msa[,xx_inc]),num_leaves,length(xx_inc))
    }
    
    sc <- matrix(0,num_leaves,length(xx_inc))
    for (j in 1:num_leaves) {
      # INCOMP
      if (length(which(pos_inc==pos_inc[j]))==1 && pos_inc[j]!=21){
        # Check 5 before, 5 after j
        for (k in 1:length(xx_inc)){
          sc[j,k] <- length(which(pos_bef_inc[,k]==pos_bef_inc[j,k]))
        }
        if (length(which(sc[j,]==1))>=length(xx_inc)/2){
          trim_msa_incomp[ind_inc,i] <- names(msa[j,i])
          ind_inc <- ind_inc + 1
        }
      }
      
      # BLOCKS
      if (msa[j,i]!="-") {
        pos_bef <- seq((i-1),(i-leng),-1)
        elim <- c()
        for (ii in 1:leng){
          if (pos_bef[ii]<=0){
            elim <- c(elim, ii)
          }
        }
        if (length(elim)>0){
          pos_bef <- pos_bef[-elim]
        }
        if (length(pos_bef)>=1){
          obs_bef <- msa[j, pos_bef]
          inds1 <- which(obs_bef=="-")
          if (length(inds1)>=1) {
            if (obs_bef[1]!="-"){
              obs_bef2 <- obs_bef[-c(1:(inds1[1]-1))]
              inds2 <- which(obs_bef2!="-")
              if (length(inds2)>=1){
                inds <- inds2[1] + inds1[1] - 1
                pos_bef <- pos_bef[-c(inds:length(obs_bef))]
              }
              
            } else {
              obs_bef2 <- obs_bef
              inds2 <- which(obs_bef2!="-")
              if (length(inds2)>=1){
                inds <- inds2[1] + inds1[1] - 1
                pos_bef <- pos_bef[-c(inds:length(obs_bef))]
              }
            }
            
          }
        }
        pos_aft <- seq((i+1),(i+leng),1)
        elim <- c()
        for (ii in 1:leng){
          if (pos_aft[ii]>total_pos){
            elim <- c(elim, ii)
          }
        }
        if (length(elim)>0){
          pos_aft <- pos_aft[-elim]
        }
        if (length(pos_aft)>=1){
          obs_aft <- msa[j, pos_aft]
          inds1 <- which(obs_aft=="-")
          if (length(inds1)>=1) {
            if (obs_aft[1]!="-"){
              obs_aft2 <- obs_aft[-c(1:(inds1[1]-1))]
              inds2 <- which(obs_aft2!="-")
              if (length(inds2)>=1){
                inds <- inds2[1] + inds1[1] - 1
                pos_aft <- pos_aft[-c(inds:length(obs_aft))]
              }
              
            } else {
              obs_aft2 <- obs_aft
              inds2 <- which(obs_aft2!="-")
              if (length(inds2)>=1){
                inds <- inds2[1] + inds1[1] - 1
                pos_aft <- pos_aft[-c(inds:length(obs_aft))]
              }
            }
            
          }
        }
        if (length(pos_aft)==0) {
          gp <- length(which(msa[j,pos_bef]=="-"))
          if (gp>=thr_gap) {
            cand[j,i] <- 1
          }
        } else if (length(pos_bef)==0) {
          gp <- length(which(msa[j,pos_aft]=="-"))
          if (gp>=thr_gap) {
            cand[j,i] <- 1
          }
        } else {
          gp <- length(which(msa[j,pos_bef]=="-"))
          gp2 <- length(which(msa[j,pos_aft]=="-"))
          if (gp>=thr_gap || gp2>=thr_gap) {
            cand[j,i] <- 1
          }
        }
        
        if (cand[j,i]==1) {
          # # print(c(i,j))
          char_ij <- msa[j,i]
          if (length(which(msa[,i] == char_ij))==1) {
            # Check the initial and final positions???
            # If they are correct do not delete
            # If one correct other not, check the overall ratio
            pos_bef <- seq((i-1),(i-leng),-1)
            elim <- c()
            for (ii in 1:leng){
              if (pos_bef[ii]<=0){
                elim <- c(elim, ii)
              }
            }
            if (length(elim)>0){
              pos_bef <- pos_bef[-elim]
            }
            if (length(pos_bef)>=1){
              obs_bef <- msa[j, pos_bef]
              inds1 <- which(obs_bef=="-")
              if (length(inds1)>=1) {
                if (obs_bef[1]!="-"){
                  obs_bef2 <- obs_bef[-c(1:(inds1[1]-1))]
                  inds2 <- which(obs_bef2!="-")
                  if (length(inds2)>=1){
                    inds <- inds2[1] + inds1[1] - 1
                    pos_bef <- pos_bef[-c(inds:length(obs_bef))]
                  }
                  
                } else {
                  obs_bef2 <- obs_bef
                  inds2 <- which(obs_bef2!="-")
                  if (length(inds2)>=1){
                    inds <- inds2[1] + inds1[1] - 1
                    pos_bef <- pos_bef[-c(inds:length(obs_bef))]
                  }
                }
                
              }
            }
            left <- 0
            obs_bef <- msa[j, pos_bef]
            plc <- which(obs_bef=="-")
            if (length(plc)>1){
              if (length(which(msa[,pos_bef[(min(plc)-1)]] == msa[j, pos_bef[(min(plc)-1)]]))==1){
                left <- 1
              }
            }
            
            pos_aft <- seq((i+1),(i+leng),1)
            elim <- c()
            for (ii in 1:leng){
              if (pos_aft[ii]>total_pos){
                elim <- c(elim, ii)
              }
            }
            if (length(elim)>0){
              pos_aft <- pos_aft[-elim]
            }
            if (length(pos_aft)>=1){
              obs_aft <- msa[j, pos_aft]
              inds1 <- which(obs_aft=="-")
              if (length(inds1)>=1) {
                if (obs_aft[1]!="-"){
                  obs_aft2 <- obs_aft[-c(1:(inds1[1]-1))]
                  inds2 <- which(obs_aft2!="-")
                  if (length(inds2)>=1){
                    inds <- inds2[1] + inds1[1] - 1
                    pos_aft <- pos_aft[-c(inds:length(obs_aft))]
                  }
                  
                } else {
                  obs_aft2 <- obs_aft
                  inds2 <- which(obs_aft2!="-")
                  if (length(inds2)>=1){
                    inds <- inds2[1] + inds1[1] - 1
                    pos_aft <- pos_aft[-c(inds:length(obs_aft))]
                  }
                }
                
              }
            }
            obs_aft <- msa[j, pos_aft]
            right <- 0
            plc <- which(obs_aft=="-")
            if (length(plc)>1){
              if (length(which(msa[,pos_aft[(min(plc)-1)]] == msa[j, pos_aft[(min(plc)-1)]]))==1){
                right <- 1
              }
            }
            if (left+right >= 1){
              cand2[j,i] <- 2
              # print(c(i,j))
            } else {
              cand2[j,i] <- 1
            }
          } else if (length(which(msa[,i] == char_ij))>1) {
            cand2[j,i] <- 1
          }
        }
      }
    }
    
  }
  
  
  locs <- which(cand2==2, arr.ind=TRUE)
  if (length(locs)>=2){
    ks <- matrix(1,total_pos,1)
    for (i in 1:length(locs[,2])){
      ps <- locs[i,2]
      sps <- locs[i,1]
      trim_msa_blocks[ks[ps],ps] <- rownames(msa)[sps]
      ks[ps] <- ks[ps]+1
    }
  }
  
  
  trim_msa <- list()
  trim_msa$gaps <- trim_msa_gaps
  trim_msa$incomp <- trim_msa_incomp
  trim_msa$blocks <- trim_msa_blocks
  
  return(trim_msa)
  
}

compute_score(file_fasta=args[1], output_name=args[2], folder_name = args[3])
