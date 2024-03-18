library(PRROC)
library(AUC)

### Read the experimental data file
file_name <- "Data/urn_mavedb_00000108-a-3_scores.csv"
main <- read.csv(file_name)

### Three letters code for each amino acid
thr_let <- read.table("Data/AA_3Let.txt")

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
submitted <- "Data/EvoPower_model_1.txt"
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
result_PHACT <- read.csv("Data/P08397_wl_param_0.csv")
result_PHACT <- as.data.frame(result_PHACT)

phact <- c()
for (i in 1:length(vect$Id)){
  i1 <- which(result_PHACT[,1]==vect$Pos[i])
  i3 <- which(colnames(result_PHACT)==vect$Alt[i])
  
  sc <- result_PHACT[i1, i3]
  
  phact <- c(phact, sc)
}

eps <- 10^(-15)
vect$PHACT <- 1-log(phact+eps)/log(eps)
vect$SUBMITTED <- vect_submitted$Score

result_table <- matrix(0, 2, 4)
result_table <- as.data.frame(result_table)

rownames(result_table) <- c("DraftApp", "PHACT")
colnames(result_table) <- c("Tau", "Spearman", "Dele_Roc", "Wild_Roc")

### KEANDALL TAU

tau_sub <- cor.test(as.numeric(vect$Score), as.numeric(vect$SUBMITTED), method="kendall")
tau_sub <- tau_sub$estimate
tau_phact <- cor.test(as.numeric(vect$Score), as.numeric(vect$PHACT), method="kendall")
tau_phact <- tau_phact$estimate

result_table$Tau <- c(tau_sub, tau_phact)

### SPEARMAN

spearman_sub <- cor.test(as.numeric(vect$Score), as.numeric(vect$SUBMITTED), method="spearman")
spearman_sub <- spearman_sub$estimate

spearman_phact <- cor.test(as.numeric(vect$Score), as.numeric(vect$PHACT), method="spearman")
spearman_phact <- spearman_phact$estimate

result_table$Spearman <- c(spearman_sub, spearman_phact)

### ROC Comparison

i1 <- which(as.numeric(vect$Score)<0.3)
i2 <- intersect(which(as.numeric(vect$Score)>=0.3), which(as.numeric(vect$Score)<0.8))
i3 <- which(as.numeric(vect$Score)>=0.8)

##### DELE_ROC

vect3 <- rbind(vect[i1,], vect[i2,], vect[i3,])
y <- c(matrix(1, 1, length(i1)), matrix(-1, 1, length(i2)), matrix(-1, 1, length(i3)))

roc1 <- roc(1-as.numeric(vect3$SUBMITTED), as.factor(1 * (y == 1)))
auc_sub <- auc(roc1)
roc2 <- roc(1-as.numeric(vect3$PHACT), as.factor(1 * (y == 1)))
auc_phact <- auc(roc2)

result_table$Dele_Roc <- c(auc_sub, auc_phact)

### WILD_ROC

vect3 <- rbind(vect[i1,], vect[i2,], vect[i3,])
y <- c(matrix(1, 1, length(i1)), matrix(1, 1, length(i2)), matrix(-1, 1, length(i3)))
  
roc1 <- roc(1-as.numeric(vect3$SUBMITTED), as.factor(1 * (y == 1)))
auc_sub <- auc(roc1)
roc2 <- roc(1-as.numeric(vect3$PHACT), as.factor(1 * (y == 1)))
auc_phact <- auc(roc2)

result_table$Wild_Roc <- c(auc_sub, auc_phact)

print(result_table)

### EXTRA MEASURES

table_extra <- matrix(0, 2, 3)
table_extra <- as.data.frame(table_extra)

rownames(table_extra) <- c("DraftApp", "PHACT")
colnames(table_extra) <- c("RMSD", "Pearson", "Value_Diff")

##### RMSD

rmsd_sub <- sqrt(sum((as.numeric(vect$PHACT) - as.numeric(vect$Score))^2)/length(vect$Id))
rmsd_phact <- sqrt(sum((as.numeric(vect$SUBMITTED) - as.numeric(vect$Score))^2)/length(vect$Id))

table_extra$RMSD <- c(rmsd_sub, rmsd_phact)

##### PEARSON

pearson_sub <- cor.test(as.numeric(vect$Score), as.numeric(vect$SUBMITTED), method="pearson")
pearson_phact <- cor.test(as.numeric(vect$Score), as.numeric(vect$PHACT), method="pearson")

table_extra$Pearson <- c(pearson_sub, pearson_phact)

##### VALUE_DIFF

calculate_AUC <- function(mutants, cutoff) {
  percentage <- sum(mutants <= cutoff) / length(mutants)
  return(percentage)
}

cutoffs <- seq(0, 1, by = 0.01)
mutants <- abs(as.numeric(vect$Score)-as.numeric(vect$SUBMITTED))
AUC_values <- sapply(cutoffs, function(cut) calculate_AUC(mutants, cut))
value_diff_sub <- sum(AUC_values)/length(AUC_values)

mutants <- abs(as.numeric(vect$Score)-as.numeric(vect$PHACT))
AUC_values <- sapply(cutoffs, function(cut) calculate_AUC(mutants, cut))
value_diff_phact <- sum(AUC_values)/length(AUC_values)

table_extra$Value_Diff <- c(value_diff_sub, value_diff_phact)
print(table_extra)


