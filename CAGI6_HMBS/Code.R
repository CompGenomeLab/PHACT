### Read the experimental data file
file_name <- "urn_mavedb_00000108-a-3_scores.csv"
main <- read.csv(file_name)

### Three letters code for each amino acid
thr_let <- read.table("AA_3Let.txt")

### Convert experimental results to a compatible format to match with PHACT scores
vect <- c()
for (i in 1:length(main$accession)){
  var <- main$hgvs_pro[i]
  var <- substr(var, 3, nchar(var))
  ref <- substr(var, 1, 3)
  ref <- thr_let[which(thr_let[,2]==ref), 1]
  var <- substr(var, 4, nchar(var))
  score1 <- main$score[i]
  score2 <- main$se[i]
  
  if (substr(var, nchar(var), nchar(var))=="="){
    alt <- ref
    pos <- as.numeric(substr(var, 1, (nchar(var)-1)))
  } else if (substr(var, (nchar(var)-2), nchar(var))!="Ter") {
    alt <- substr(var, (nchar(var)-2), nchar(var))
    alt <- thr_let[which(thr_let[,2]==alt), 1]
    pos <- as.numeric(substr(var, 1, (nchar(var)-3)))
  } else {
    alt <- "Ter"
    pos <- as.numeric(substr(var, 1, (nchar(var)-3)))
  }
  vect <- rbind(vect, c("P08397", pos, ref, alt, score1, score2))
  if (length(c("P08397", pos, ref, alt, score1, score2))!=6){
    print(i)
  }
}

vect <- as.data.frame(vect)
colnames(vect) <- c("Id", "Pos", "Ref", "Alt", "Score", "StdErr")

### Our submission to CAGI6 including draft version of PHACT approach
submitted <- "/Users/nurdankuru/Downloads/EvoPower_model_1.txt"
submitted_data <- read.table(submitted)
vect_sub <- c()
for (i in 1:length(submitted_data$aa_substitution)){
  var <- submitted_data$aa_substitution[i]
  var <- substr(var, 3, nchar(var))
  ref <- substr(var, 1, 3)
  ref <- thr_let[which(thr_let[,2]==ref), 1]
  var <- substr(var, 4, nchar(var))
  score1 <- submitted_data$score[i]
  score2 <- submitted_data$sd[i]
  
  if (substr(var, nchar(var), nchar(var))=="="){
    alt <- ref
    pos <- as.numeric(substr(var, 1, (nchar(var)-1)))
  } else if (substr(var, (nchar(var)-2), nchar(var))!="Ter") {
    alt <- substr(var, (nchar(var)-2), nchar(var))
    alt <- thr_let[which(thr_let[,2]==alt), 1]
    pos <- as.numeric(substr(var, 1, (nchar(var)-3)))
  } else {
    alt <- "Ter"
    pos <- as.numeric(substr(var, 1, (nchar(var)-3)))
  }
  vect_sub <- rbind(vect_sub, c("P08397", pos, ref, alt, score1, score2))
  if (length(c("P08397", pos, ref, alt, score1, score2))!=6){
    print(i)
  }
}
vect_sub <- as.data.frame(vect_sub)
colnames(vect_sub) <- c("Id", "Pos", "Ref", "Alt", "Score", "StdErr")

### Keep only missense mutations

vect <- vect[-which(vect$Alt=="Ter"),]
vect_sub <- vect_sub[-which(vect_sub$Alt=="Ter"),]
vect2 <- vect[-which(vect$Ref==vect$Alt),]
vect2_sub <- vect_sub[-which(vect_sub$Ref==vect_sub$Alt),]

### Remove the variants with >1.36 experimental score as mentioned in Zhang et al., 2024
vect2 <- vect2[-which(as.numeric(vect2$Score)>1.36),]

### Combine experimental data with the submitted predictions
vect2$var <- paste(vect2$Pos, vect2$Ref, vect2$Alt, sep = "")
vect2_sub$var <- paste(vect2_sub$Pos, vect2_sub$Ref, vect2_sub$Alt, sep = "")
com <- intersect(vect2$var, vect2_sub$var)
vect2_sub <- vect2_sub[match(vect2$var, vect2_sub$var),]

vect <- vect2
vect$SUBMITTED <- vect_sub$Score

### Get PHACT results
result_PHACT <- read.csv("P08397_wl_param_0.csv")
result_PHACT <- as.data.frame(result_PHACT)

phact <- c()
for (i in 1:length(vect$Id)){
  i1 <- which(result_PHACT[,1]==vect$Pos[i])
  i3 <- which(colnames(result_PHACT)==vect$Alt[i])
  
  sc <- result_PHACT[i1, i3]
  
  phact <- c(phact, sc)
}

vect$PHACT <- phact




















