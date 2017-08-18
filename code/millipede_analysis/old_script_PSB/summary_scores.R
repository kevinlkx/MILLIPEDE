rm(list = ls())

# source("~/DNaseData/script/lr_functions.R")
require(ROCR)
dir_results_human = paste('~/DNaseData/Results_forPSB/summary/human/', sep ='')
dir_results_yeast = paste('~/DNaseData/Results_forPSB/summary/yeast/', sep ='')
dir.create(dir_results_yeast, showWarnings = FALSE, recursive = TRUE)
dir.create(paste(dir_results_human, "site/", sep =''), showWarnings = FALSE, recursive = TRUE)
dir.create(paste(dir_results_human, "peak/", sep =''), showWarnings = FALSE, recursive = TRUE)


# weights = "Trained"
# k = 5
# humanTFs = c("REST", "MAX", "SRF")
humanTFs = c("REST", "CTCF", "MAX", "SRF", "GABPA", "JUND")
yeastTFs = c("reb1", "rap1",  "gal4", "abf1",
             "cbf1", "mcm1", "fkh1", "fkh2",
             "pho4", "swi4", "sfp1", "fhl1",
             "swi5", "pho2", "leu3", "stb2",
                           "phd1", 
             "ace2", "gcn4", "mbp1", "met32"
)


# Capitalizing - toupper every first letter of the TF names
simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2), sep="", collapse=" ")
}

yeast_scores = function(k, weights){
  if (weights == "Trained"){
    yeastTFs = yeastTFs[-1]
  }
  AUROC = AUPRC = Sens = Prec = matrix(0, ncol = length(yeastTFs), nrow = 10)
  
  for (i in 1:length(yeastTFs)){
    tf = yeastTFs[i]
    tfName = paste(tf,"_macisaac_score4",sep = "")
    scores = read.table(paste("~/DNaseData/Results_forPSB/", tfName, "/", simpleCap(tf), "_fold", k,"_weights", weights,"_scores.txt", sep = ''),sep = "\t", header = T)
    AUROC[,i] = scores$AUROC
    AUPRC[,i] = scores$AUPRC
    Sens[,i] = scores$Sens
    Prec[,i] = scores$Prec
    
  }
  colnames(AUROC) = colnames(AUPRC) = colnames(Sens) = colnames(Prec) = yeastTFs
  AUROC = data.frame(TF = scores$Models, AUROC, mean = round(rowMeans(AUROC),3), sd = round(apply(AUROC,1,sd),3))
  AUPRC = data.frame(TF = scores$Models, AUPRC, mean = round(rowMeans(AUPRC),3), sd = round(apply(AUPRC,1,sd),3))
  Sens = data.frame(TF = scores$Models, Sens, mean = round(rowMeans(Sens),3), sd = round(apply(Sens,1,sd),3))
  Prec = data.frame(TF = scores$Models, Prec, mean = round(rowMeans(Prec),3), sd = round(apply(Prec,1,sd),3))
  
  fileName = paste(dir_results_yeast, "AUROC", "_fold", k, "_weights", weights, ".txt", sep = '')
  write.table(t(AUROC), file = fileName, sep = "\t", quote = F, col.names = F, row.names = T)
  fileName = paste(dir_results_yeast, "AUPRC", "_fold", k, "_weights", weights, ".txt", sep = '')
  write.table(t(AUPRC), file = fileName, sep = "\t", quote = F, col.names = F, row.names = T)
  fileName = paste(dir_results_yeast, "Sens", "_fold", k, "_weights", weights, ".txt", sep = '')
  write.table(t(Sens), file = fileName, sep = "\t", quote = F, col.names = F, row.names = T)
  fileName = paste(dir_results_yeast, "Prec", "_fold", k, "_weights", weights, ".txt", sep = '')
  write.table(t(Prec), file = fileName, sep = "\t", quote = F, col.names = F, row.names = T)
}

human_scores = function(goldstandard, k, weights){
  if (weights == "Trained"){
    humanTFs = humanTFs[-1]
  }
  AUROC = AUPRC = Sens = Prec = matrix(0, ncol = length(humanTFs), nrow = 11)
  
  for (i in 1:length(humanTFs)){
    tfName = humanTFs[i]
    scores = read.table(paste("~/DNaseData/Results_forPSB/", tfName, "/", tfName,"_", goldstandard, "_fold", k,"_weights", weights,"_scores.txt", sep = ''),sep = "\t", header = T)
    AUROC[,i] = scores$AUROC
    AUPRC[,i] = scores$AUPRC
    Sens[,i] = scores$Sens
    Prec[,i] = scores$Prec
    
  }
  colnames(AUROC) = colnames(AUPRC) = colnames(Sens) = colnames(Prec) = humanTFs
  AUROC = data.frame(TF = scores$Models, AUROC, mean = round(rowMeans(AUROC),3), sd = round(apply(AUROC,1,sd),3))
  AUPRC = data.frame(TF = scores$Models, AUPRC, mean = round(rowMeans(AUPRC),3), sd = round(apply(AUPRC,1,sd),3))
  Sens = data.frame(TF = scores$Models, Sens, mean = round(rowMeans(Sens),3), sd = round(apply(Sens,1,sd),3))
  Prec = data.frame(TF = scores$Models, Prec, mean = round(rowMeans(Prec),3), sd = round(apply(Prec,1,sd),3))
  
  fileName = paste(dir_results_human, goldstandard, "/", "AUROC","_fold", k, "_weights", weights, ".txt", sep = '')
  write.table(t(AUROC), file = fileName, sep = "\t", quote = F, col.names = F, row.names = T)
  fileName = paste(dir_results_human, goldstandard, "/", "AUPRC", "_fold", k, "_weights", weights, ".txt", sep = '')
  write.table(t(AUPRC), file = fileName, sep = "\t", quote = F, col.names = F, row.names = T)
  fileName = paste(dir_results_human,  goldstandard, "/", "Sens", "_fold", k, "_weights", weights, ".txt", sep = '')
  write.table(t(Sens), file = fileName, sep = "\t", quote = F, col.names = F, row.names = T)
  fileName = paste(dir_results_human,  goldstandard, "/", "Prec", "_fold", k, "_weights", weights, ".txt", sep = '')
  write.table(t(Prec), file = fileName, sep = "\t", quote = F, col.names = F, row.names = T)
}

yeast_scores(k = 5, weights = "LR")
# yeast_scores(k = 1, weights = "Trained")
# yeast_scores(k = 1, weights = "Equal")

human_scores("site", k = 5, weights = "LR")
human_scores("site", k = 1, weights = "Trained")
human_scores("site", k = 1, weights = "Equal")

human_scores("peak", k = 5, weights = "LR")
# human_scores("peak", k = 1, weights = "Trained")
# human_scores("peak", k = 1, weights = "Equal")