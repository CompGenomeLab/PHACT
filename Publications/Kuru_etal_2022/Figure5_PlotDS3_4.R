library(AUC)
library(PRROC)

# Update Lines 6&7&9 depending on whether DS3 or DS4 is used.

output_roc <- "DS4_RocCurves.pdf"
output_pr <- "DS4_PRCurves.pdf"

data <- readxl::read_xlsx("Dataset_DS3.xlsx")
# data <- readxl::read_xlsx("Dataset_DS4.xlsx")
y <- as.numeric(data$Variant_type)

cols <- c("#ca0020", "#4393c3", "#2166ac", "#ffa056")

auc_val <- list()
pdf(file = sprintf(output_roc), width = 4, height = 4)
par(oma = c(2, 1, 1, 1), mar = c(3,3,2,1), cex = (11/12), lwd = 1, cex.axis = 1.0, cex.main = 1)
auc_values <- c()
ind <- 1
ind <-  ind + 1
sift_ourmsa_roc <- roc(1 - data$SIFT_phact_msa, as.factor(1 * (y == 1)))
plot(sift_ourmsa_roc, col = cols[ind], lwd = 2)
# lines(sift_ourmsa_roc$fpr, sift_ourmsa_roc$tpr, col = cols[ind], lwd = 2)
text(0.5, 0.5 - (ind - 1) * 0.07, label = sprintf("SIFT - PHACT's MSA = %.3f", auc(sift_ourmsa_roc)), col = cols[ind], cex = 0.9)

ind <-  ind + 1
sift <- roc(data$SIFT_converted_rankscore, as.factor(1 * (y == 1)))
lines(sift$fpr, sift$tpr, col = cols[ind], lwd = 2)
text(0.5, 0.5 - (ind - 1) * 0.07, label = sprintf("SIFT = %.3f", auc(sift)), col = cols[ind], cex = 1)

ind <-  ind + 1
pph2 <- roc(data$Polyphen2_HVAR_rankscore, as.factor(1 * (y == 1)))
lines(pph2$fpr, pph2$tpr, col = cols[ind], lwd = 2)
text(0.5, 0.5 - (ind - 1) * 0.07, label = sprintf("PPH2_HVAR = %.3f", auc(pph2)), col = cols[ind], cex = 1)


ind <- 1
phylas_roc <- roc(1 - data$`PHACT Score`, as.factor(1 * (y == 1)))
lines(phylas_roc$fpr, phylas_roc$tpr, col = cols[ind], lwd = 2)
text(0.5, 0.5 - (ind - 1) * 0.07, label = sprintf("PHACT = %.3f", auc(phylas_roc)), col = cols[1], cex = 1)

mtext(side = 2, at = 0.5 , text = "TPR (sensitivity)", line = 2.5, cex = 1)

mtext(side = 1, at = 0.5 , text = "FPR (1 - specificity)", line = 2.5, cex = 1)

dev.off()  


auc_val <- list()
pdf(file = sprintf(output_pr), width = 4, height = 4)
par(oma = c(2, 1, 1, 1), mar = c(3,3,2,1), cex = (11/12), lwd = 1, cex.axis = 1.0, cex.main = 1)

auc_values <- c()
ind <- 1

ind <-  ind + 1

sift_ourmsa_prauc <- pr.curve(scores.class0 = 1 - data$SIFT_phact_msa, weights.class0 = (1 * (y == 1)), curve = TRUE)
plot(sift_ourmsa_prauc, col = cols[ind], lwd = 2, main = "", auc.main = FALSE)
text(0.5, 0.3 - (ind - 1) * 0.07, label = sprintf("SIFT - PHACT's MSA = %.3f", (sift_ourmsa_prauc$auc.integral)), col = cols[ind], cex = 1)

ind <-  ind + 1

sift <- pr.curve(scores.class0 = data$SIFT_converted_rankscore, weights.class0 = (1 * (y == 1)), curve = TRUE)
lines(sift$curve[,1], sift$curve[,2], col = cols[ind], lwd = 2)
text(0.5, 0.3 - (ind - 1) * 0.07, label = sprintf("SIFT = %.3f", (sift$auc.integral)), col = cols[ind], cex = 1)

ind <-  ind + 1
pph2 <- pr.curve(scores.class0 = data$Polyphen2_HVAR_rankscore, weights.class0 = (1 * (y == 1)), curve = TRUE)
lines(pph2$curve[,1], pph2$curve[,2], col = cols[ind], lwd = 2)
text(0.5, 0.3 - (ind - 1) * 0.07, label = sprintf("PPH2_HVAR = %.3f", (pph2$auc.integral)), col = cols[ind], cex = 1)

ind <- 1
phylas_prauc <- pr.curve(scores.class0 = 1 - data$`PHACT Score`, weights.class0 = (1 * (y == 1)), curve = TRUE)
lines(phylas_prauc$curve[,1], phylas_prauc$curve[,2], col = cols[ind], lwd = 2)
text(0.5, 0.3 - (ind - 1) * 0.07, label = sprintf("PHACT = %.3f", phylas_prauc$auc.integral), col = cols[1], cex = 1)


mtext(side = 2, at = 0.5 , text = "Precision", line = 2.5, cex = 1)

mtext(side = 1, at = 0.5 , text = "Recall", line = 2.5, cex = 1)

dev.off() 
