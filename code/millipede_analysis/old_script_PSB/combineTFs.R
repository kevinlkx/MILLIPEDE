# source("~/DNaseData/script/lr_functions.R")
require(ROCR)
dir_results_human = paste('~/DNaseData/Results_forPSB/summary/human/', sep ='')
dir_results_yeast = paste('~/DNaseData/Results_forPSB/summary/yeast/', sep ='')
dir.create(dir_results_human, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_results_yeast, showWarnings = FALSE, recursive = TRUE)
dir.create(paste(dir_results_human, "site/", sep =''), showWarnings = FALSE, recursive = TRUE)
dir.create(paste(dir_results_human, "peak/", sep =''), showWarnings = FALSE, recursive = TRUE)

weights = "LR"
k = 5

humanTFs = c("REST", "CTCF", "MAX", "SRF", "GABPA", "JUND")
yeastTFs = c("reb1", "rap1",  "gal4", "abf1",
             "cbf1", "mcm1", "fkh1", "fkh2",
             "phd1", "swi4", "sfp1", "fhl1",
             "swi5", "pho2", "leu3", "stb2",
             "ace2", "gcn4", "mbp1", "met32"
)

# MpdColors = c("pink", "orange", "red", "purple", "magenta", "darkred")
# CpdColors = c("blue", "lightblue", "darkblue", "cyan", "green", "lightgreen")

human_curves = function(goldstandard,shrinkage, Mpd){
  
  human_results.l= list()
  
  for (i in 1:length(humanTFs)){
    tfName = humanTFs[i]
    results.l = get(load(paste("~/DNaseData/Results_forPSB/", tfName, "/", goldstandard, "/weightsLRresults.l", sep = '')))
    results.df = results.l$results
    human_results.l[[i]] = list(
      tfName = tfName,
      PostPr_NoShrinkage = results.df$PostPr_NoShrinkage, 
      PostPr_Shrinkage = results.df$PostPr_Shrinkage, 
      M5 = results.df$M5, 
      M3 = results.df$M3,
      M11 = results.df$M11,
      label = results.l$label)
  }
  
  # fileName = paste(dir_results, "yeast_results", "_fold", k, "_weights", weights, ".l", sep = '')
  # save(yeast_results.l, file = fileName)
  # yeast_results.l = get(load(fileName))
  
  fileName = paste(dir_results_human, goldstandard, "/fold", k, "_weights", weights, "_shrinkage", shrinkage, "_", Mpd, "_ROCs", sep = "")
  pdf(file= paste(fileName, ".pdf", sep = "")) 
  par(xaxs = "i")
  par(yaxs = "i")
  plot(x = c(0,1), y = c(0,1), type = "n", ylab = "Sensitivity", xlab = "False positive rate")
  mtext(paste("ROC curves for human TFs"),3,line=1, cex = 1.5)
  multiROC(human_results.l, shrinkage, Mpd)
  modelNames = c("MILLIPEDE M5", "CENTIPEDE w/ shrinkage")
  legend("bottomright", legend = modelNames, lty = 1, lwd=2, col= c("red", "blue"), bty="n", cex = 0.8)
# #   dev.off()
#   
# #   fileName = paste(dir_results, "humanTFs_", "fold", k, "_weights", weights, "_shrinkage", shrinkage, "_", Mpd, "_Sens", sep = "")
# #   pdf(file= paste(fileName, ".pdf", sep = "")) 
# #   par(xaxs = "i")
# #   par(yaxs = "i")
#   plot(x = c(0,0.01), y = c(0,1), type = "n", ylab = "Sensitivity", xlab = "False positive rate")
#   mtext(paste("ROC curves for human TFs"),3,line=1,font=2, cex = 1)
#   multiROC(human_results.l, shrinkage, Mpd)
#   modelNames = c("MILLIPEDE M5", "CENTIPEDE w/ shrinkage")
#   legend("topleft", legend = modelNames, lty = 1, lwd=2, col= c("red", "blue"), bty="n", cex = 0.8)
# #   dev.off()
#   
# #   fileName = paste(dir_results, "humanTFs_", "fold", k, "_weights", weights, "_shrinkage", shrinkage, "_", Mpd, "_PRCs", sep = "")
# #   pdf(file= paste(fileName, ".pdf", sep = "")) 
# #   par(xaxs = "i")
# #   par(yaxs = "i")
#   plot(x = c(0,1), y = c(0,1), type = "n", ylab = "Precision", xlab = "Recall")
#   mtext(paste("PR curves for human TFs"),3,line=1,font=2, cex = 1)
#   multiPRC(human_results.l, shrinkage, Mpd)
#   modelNames = c("MILLIPEDE M5", "CENTIPEDE w/ shrinkage")
#   legend("bottomleft", legend = modelNames, lty = 1, lwd=2, col= c("red", "blue"), bty="n", cex = 0.8)
  
  dev.off()
  
  cat("Human curves finished ...\n")
  
  
}


yeast_curves = function(shrinkage, Mpd){
  
  yeast_results.l= list()
  
  for (i in 1:length(yeastTFs)){
    tf = yeastTFs[i]
    tfName = paste(tf,"_macisaac_score4",sep = "")
    results.l = get(load(paste("~/DNaseData/Results_forPSB/", tfName, "/weightsLR_results.l", sep = '')))
    results.df = results.l$results
    yeast_results.l[[i]] = list(
      tfName = tf,
      PostPr_NoShrinkage = results.df$PostPr_NoShrinkage, 
      PostPr_Shrinkage = results.df$PostPr_Shrinkage, 
      M5 = results.df$M5, 
      M3 = results.df$M3,
      M11 = results.df$M11,
      label = results.l$label)
  }
  
  # fileName = paste(dir_results, "yeast_results", "_fold", k, "_weights", weights, ".l", sep = '')
  # save(yeast_results.l, file = fileName)
  # yeast_results.l = get(load(fileName))
  
  fileName = paste(dir_results_yeast, "fold", k, "_weights", weights, "_shrinkage", shrinkage, "_", Mpd, "_ROCs", sep = "")
  pdf(file= paste(fileName, ".pdf", sep = "")) 
  par(xaxs = "i")
  par(yaxs = "i")
  plot(x = c(0,1), y = c(0,1), type = "n", ylab = "Sensitivity", xlab = "False positive rate")
  mtext(paste("ROC curves for yeast TFs"),3,line=1, cex = 2)
  multiROC(yeast_results.l, shrinkage, Mpd)
  modelNames = c("MILLIPEDE M5", "CENTIPEDE w/ shrinkage")
  legend("bottomright", legend = modelNames, lty = 1, lwd=2, col= c("red", "blue"), bty="n", cex = 0.8)
  
  #   dev.off()
  
  #   fileName = paste(dir_results, "humanTFs_", "fold", k, "_weights", weights, "_shrinkage", shrinkage, "_", Mpd, "_Sens", sep = "")
  #   pdf(file= paste(fileName, ".pdf", sep = "")) 
  #   par(xaxs = "i")
  #   par(yaxs = "i")
#   plot(x = c(0,0.01), y = c(0,1), type = "n", ylab = "Sensitivity", xlab = "False positive rate")
#   mtext(paste("ROC curves for yeast TFs"),3,line=1,font=2, cex = 1)
#   multiROC(yeast_results.l, shrinkage, Mpd)
#   modelNames = c("MILLIPEDE M5", "CENTIPEDE w/ shrinkage")
#   legend("topleft", legend = modelNames, lty = 1, lwd=2, col= c("red", "blue"), bty="n", cex = 0.8)
#   #   dev.off()
#   
#   #   fileName = paste(dir_results, "humanTFs_", "fold", k, "_weights", weights, "_shrinkage", shrinkage, "_", Mpd, "_PRCs", sep = "")
#   #   pdf(file= paste(fileName, ".pdf", sep = "")) 
#   #   par(xaxs = "i")
#   #   par(yaxs = "i")
#   plot(x = c(0,1), y = c(0,1), type = "n", ylab = "Precision", xlab = "Recall")
#   mtext(paste("PR curves for yeast TFs"),3,line=1,font=2, cex = 1)
#   multiPRC(yeast_results.l, shrinkage, Mpd)
#   modelNames = c("MILLIPEDE M5", "CENTIPEDE w/ shrinkage")
#   legend("bottomleft", legend = modelNames, lty = 1, lwd=2, col= c("red", "blue"), bty="n", cex = 0.8)
  
  dev.off()
  
  cat("Yeast ROCs finished ... \n")
  
  
}



multiROC = function(results.l, shrinkage, Mpd){
  
  for (i in 1:length(results.l)){
    if (shrinkage == T){
      perf_ROC <- performance(prediction(results.l[[i]]$PostPr_Shrinkage, results.l[[i]]$label), measure='sens',x.measure='fpr')
    }else{
      perf_ROC <- performance(prediction(results.l[[i]]$PostPr_NoShrinkage, results.l[[i]]$label), measure='sens',x.measure='fpr')
    }
    plot(perf_ROC, col = "blue", lwd=2, add = T)
    
    if (Mpd == "M3"){
      perf_ROC <- performance(prediction(results.l[[i]]$M3, results.l[[i]]$label), measure='sens',x.measure='fpr')
    }else if (Mpd == "M5"){
      perf_ROC <- performance(prediction(results.l[[i]]$M5, results.l[[i]]$label), measure='sens',x.measure='fpr')
    }else if (Mpd == "M11"){
      perf_ROC <- performance(prediction(results.l[[i]]$M11, results.l[[i]]$label), measure='sens',x.measure='fpr')
    }
    plot(perf_ROC, col = "red", lwd=2, add = T)
    
  }
  
  box()
  
  
}

multiPRC = function(results.l, shrinkage, Mpd){
  
  for (i in 1:length(results.l)){
    if (shrinkage == T){
      perf_ROC <- performance(prediction(results.l[[i]]$PostPr_Shrinkage, results.l[[i]]$label), measure='prec',x.measure='rec')
    }else{
      perf_ROC <- performance(prediction(results.l[[i]]$PostPr_NoShrinkage, results.l[[i]]$label), measure='prec',x.measure='rec')
    }
    plot(perf_ROC, col = "blue", lwd=2, add = T)
    
    if (Mpd == "M3"){
      perf_ROC <- performance(prediction(results.l[[i]]$M3, results.l[[i]]$label), measure='prec',x.measure='rec')
    }else if (Mpd == "M5"){
      perf_ROC <- performance(prediction(results.l[[i]]$M5, results.l[[i]]$label), measure='prec',x.measure='rec')
    }else if (Mpd == "M11"){
      perf_ROC <- performance(prediction(results.l[[i]]$M11, results.l[[i]]$label), measure='prec',x.measure='rec')
    }
    plot(perf_ROC, col = "red", lwd=2, add = T)
    
  }
  
  box()
  
}

yeast_curves(shrinkage = T, Mpd = "M5")
# yeast_curves(shrinkage = F, Mpd = "M5")

# human_curves("site", shrinkage = T, Mpd = "M5")
human_curves("peak", shrinkage = T, Mpd = "M5")



