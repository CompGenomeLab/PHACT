
library(readxl)
library(AUC)
library(PRROC)
library(RColorBrewer)

data <- read_xlsx("Dataset_DS5.xlsx")
data <- as.data.frame(data)
y <- data$Variant_type
colnames(data)[15] <- "PHACT"

selected_methods_HCG <- c("PHACT", "SIFT_phact_msa", "SIFT_converted_rankscore", "SIFT4G_converted_rankscore", 
                          "PROVEAN_converted_rankscore", "LIST-S2_rankscore", "GenoCanyon_rankscore", 
                          "integrated_fitCons_rankscore", "GM12878_fitCons_rankscore", 
                          "H1-hESC_fitCons_rankscore", "HUVEC_fitCons_rankscore", "GERP++_RS_rankscore", 
                          "phyloP100way_vertebrate_rankscore", "phyloP30way_mammalian_rankscore", 
                          "phyloP17way_primate_rankscore", "phastCons100way_vertebrate_rankscore", 
                          "phastCons30way_mammalian_rankscore", "phastCons17way_primate_rankscore", 
                          "SiPhy_29way_logOdds_rankscore", "bStatistic_converted_rankscore")

auroc_results <- rep(NA, length(selected_methods_HCG))
prauc_results <- rep(NA, length(selected_methods_HCG))

names(auroc_results) <- selected_methods_HCG
names(prauc_results) <- selected_methods_HCG

for(method in selected_methods_HCG) {
  if(method %in% c("PHACT", "SIFT_phact_msa") ) {
    roc_val <- roc(1 - data[,method], as.factor(1 * (y == 1)))
    auroc_results[method] <- auc(roc_val)
    
    prauc_val <- pr.curve(scores.class0 = 1 - data[,method], weights.class0 = (1 * (y == 1)), curve = TRUE)
    prauc_results[method] <- prauc_val$auc.integral
  }else{
    roc_val <- roc(data[,method], as.factor(1 * (y == 1)))
    auroc_results[method] <- auc(roc_val)
    
    prauc_val <- pr.curve(scores.class0 = data[,method], weights.class0 = (1 * (y == 1)), curve = TRUE)
    prauc_results[method] <- prauc_val$auc.integral
  }
}

names(auroc_results) <- gsub(pattern = "_rankscore", replacement = "", x =names(auroc_results))
names(auroc_results)[c(1,2)] <- c("PHACT", "SIFT - PHACT's MSA")
names(auroc_results) <- gsub(pattern = "_converted", replacement = "", x = names(auroc_results))

names(prauc_results) <- gsub(pattern = "_rankscore", replacement = "", x =names(prauc_results))
names(prauc_results)[c(1,2)] <- c("PHACT", "SIFT - PHACT's MSA")
names(prauc_results) <- gsub(pattern = "_converted", replacement = "", x = names(prauc_results))

best_names <- names(sort(auroc_results, decreasing = T))[1:20]

cols <- c(colorRampPalette(colors = c("#21457A", "#f7fbff"))((20)))

pdf(file = "Figure_7A.pdf", width = 9, height = 4, colormodel = "cmyk")
par(oma = c(2, 1, 1, 0),  mar = c(1,4,1,0), cex = 11/12, lwd = 1, cex.axis = 1.0, cex.main = 1)
x1 <- barplot(auroc_results[best_names], col = cols, ylim = c(0,1.1), horiz = F, yaxt = "n",
              las = 2, xaxt="n")
text(x1, (auroc_results[best_names]*0.5), labels = best_names, srt = 90, cex = 1, col = c(rep("white", 6), rep("black", (length(best_names)-6))))

axis(side = 2, at = seq(0,1,0.1), labels = sprintf("%.1f", seq(0,1,0.1), las = 2))

text(x1,(auroc_results[best_names] + 0.03), labels = sprintf("%.2f",auroc_results[best_names]), cex = 1)
mtext(side = 2, at = 0.5 , text = "AUC", line = 2.5, cex = 1)
dev.off()

pdf(file = "Figure_7B.pdf", width = 9, height = 4, colormodel = "cmyk")
par(oma = c(2, 1, 1, 0),  mar = c(1,4,1,0), cex = 11/12, lwd = 1, cex.axis = 1.0, cex.main = 1)
x2 <- barplot(sort(prauc_results[best_names], decreasing = T), col= cols,#[order((prauc_results[best_names]), decreasing = T)],
              yaxt = "n", xaxt = "n", ylim = c(0,1.1), horiz = F, las = 2)
text(x2, (sort(prauc_results[best_names], decreasing = T)*0.5), labels = names(sort(prauc_results[best_names], decreasing = T)), srt = 90, cex = 1, col = c(rep("white", 6), rep("black", (length(best_names)-6))))

axis(side = 2, at = seq(0,1,0.1), labels = sprintf("%.1f", seq(0,1,0.1), las = 1))
text(x2,(sort(prauc_results[best_names], decreasing = T) + 0.03), labels = sprintf("%.2f",sort(prauc_results[best_names], decreasing = T)), cex = 1)
mtext(side = 2, at = 0.5 , text = "AUPR", line = 2.5, cex = 1)
dev.off()
