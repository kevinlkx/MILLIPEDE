library(GenomicRanges)
library(gplots)
library(ROCR)

setwd("~/Dropbox/projects_dropbox/MILLIPEDE/workflow/code/millipede_analysis")

source("lr_functions.R")
source("fitCentipede_footprint.R")

# Set parameters
args <- commandArgs(trailingOnly=T);
tf <- args[1]
# tf <- "reb1"

tfName <- paste(tf,"_macisaac_score4",sep = "")
flank_data <- 100
MaxSweeps <- 400
fileType <- "pdf"

chipexoList <- c("reb1", "rap1", "gal4", "phd1")

dir_dnase_data <- paste0("../../data/yeast_data/dnase_data/", tfName, "/")
dir_results <- paste0('../../output/results/', tfName, "/")
dir.create(dir_results, showWarnings = FALSE, recursive = TRUE)

dir_motifs <- paste0('../../data/yeast_data/macisaac_score4/')
dir_binding_sites <- paste0('../../data/yeast_data/binding_sites/')

myColors <- c("pink", "orange", "red", "purple", "magenta", "darkred", "darkgray", "black", "blue", "cyan")

modelNames <- c("M12", "M11", "M5", "M3", "M2", "M1", "M1 w/o PWM", "PWM only",
               "CENTIPEDE w/ shrinkage", "CENTIPEDE w/o shrinkage")

ThreshAvgMappability <- 0.8

tFPR <- 1.0/100

millipede_analysis <- function(k, weights){

  if(sum(tf == chipexoList)){
    chipexo <- T
  }else{
    chipexo <- F
  }

  if(k > 1) {
    cv <- T
  }else{
    cv <- F
  }

  cat("Analyzing:", tfName, "fold", k, "weights", weights, "\n")

  dnase_data.l <- get(load(paste(dir_dnase_data, "dnase_cuts_", flank_data, "bp.l", sep = "")))

  tf_overlap.gr <- dnase_data.l$tf.gr

  # filter sites by mappability (remove the sites with more than > 20% unmappable bases) and start position > 0
  unmappablebases.gr <- get(load("../../data/yeast_data/mappablity/unmappablebases_sacCer2.gr"))

  cat(length(tf_overlap.gr), " -> [start position > 0 and mappable] -> ")
  select_idx <- mappability(tf_overlap.gr, unmappablebases.gr, ThreshAvgMappability)

  tf_overlap.gr <- tf_overlap.gr[select_idx]
  cat(length(tf_overlap.gr),"\n")

  if(width(tf_overlap.gr[1]) > flank_data){
    tf_sites.gr <- tf_overlap.gr
    start(tf_sites.gr) <- start(tf_overlap.gr) + flank_data
    end(tf_sites.gr) <- end(tf_overlap.gr) - flank_data
  }

  rank <- order(as.integer(seqnames(tf_sites.gr)), start(tf_sites.gr))

  DNaseData <- dnase_data.l$counts[select_idx,][rank,]
  pwm <- elementMetadata(tf_sites.gr)[,"score"][rank]
  tf_sites.gr <- tf_sites.gr[rank]

  tf_macIsaac.gr <- get(load(paste0(dir_binding_sites, "/", tf, ".gr")))

  if (chipexo == TRUE){
    tf_chipexo.gr <- get(load(paste0(dir_binding_sites, "/", tf, "_chip_exo.gr")))
    tf_sites.gr <- label_macisaac(tf_sites.gr, tf_macIsaac.gr)
    tf_sites.gr <- label_chipexo(tf_sites.gr, tf_chipexo.gr)
    elementMetadata(tf_sites.gr)[,"label_joint"] <- as.integer(elementMetadata(tf_sites.gr)[,"label"] | elementMetadata(tf_sites.gr)[,"chip_exo"])
    label <- elementMetadata(tf_sites.gr)[,"label"][rank]
    label_chip_exo <- elementMetadata(tf_sites.gr)[,"chip_exo"][rank]
    label_joint <- elementMetadata(tf_sites.gr)[,"label_joint"][rank]
    cat("macIsaac sites number:", length(tf_macIsaac.gr),
        "chip-exo sites number:", length(tf_chipexo.gr), "\n")
    cat("matched macIsaac sites number:", sum(elementMetadata(tf_sites.gr)[,"label"] == 1),
        "matched chip-exo sites number:", sum(elementMetadata(tf_sites.gr)[,"chip_exo"] == 1),
        "Joint positive sites number:", sum(elementMetadata(tf_sites.gr)[,"label_joint"] == 1),"\n")

  }else if(chipexo == F){
    tf_sites.gr <- label_macisaac(tf_sites.gr, tf_macIsaac.gr)
    elementMetadata(tf_sites.gr)[,"label_joint"] <- as.integer(elementMetadata(tf_sites.gr)[,"label"])
    label_joint <- label <- elementMetadata(tf_sites.gr)[,"label"][rank]
    cat("macIsaac sites number:", length(tf_macIsaac.gr), "\n")
    cat("matched macIsaac sites number:", sum(elementMetadata(tf_sites.gr)[,"label"] == 1),"\n")
  }

  cuts_sum.l <- sum_cuts(DNaseData)
  D_s <- as.numeric(cuts_sum.l$M1)

  ## PWM only
  pred_pwm <- prediction(pwm, label_joint)

  ## DNase score only (total DNase cuts)
  pred_D_s <- prediction(D_s, label_joint)

  ## fit centipede model
  cent.df <- centipede_dnase(DNaseData, pwm)
  PostPr_Shrinkage <- cent.df$PostPr_Shrinkage
  PostPr_NoShrinkage <- cent.df$PostPr_NoShrinkage

  pred_centipede_Shrinkage <- prediction(PostPr_Shrinkage, label_joint)
  pred_centipede_NoShrinkage <- prediction(PostPr_NoShrinkage, label_joint)

  ## logistic regression with cross-validation
  n1 <- length(rank[label_joint == 1])
  n0 <- length(rank[label_joint == 0])
  fold1 <- rep(1:k,length.out = n1)[sample(n1,n1)]
  fold0 <- rep(1:k,length.out = n0)[sample(n0,n0)]
  fold <- c(fold1,fold0)

  pred_models.l <- label.l <- list()
  M12.l <- M11.l <- M5.l <- M3.l <- M2.l <- M1.l <- list()

  for ( i in 1:k ){
    testIdx  <- which(fold == i)
    if(cv == T){
      trainIdx <- which(fold != i)
    }else{
      trainIdx <- testIdx
    }

    cat(i,"cv \n")
    pred_models.l <- lr_dnase(cuts_sum.l, pwm, label_joint, trainIdx, testIdx, weights)
    label.l[[i]] <- label_joint[testIdx]

    M12.l[[i]] <- pred_models.l$scores[[1]]
    M11.l[[i]] <- pred_models.l$scores[[2]]
    M5.l[[i]] <- pred_models.l$scores[[3]]
    M3.l[[i]] <- pred_models.l$scores[[4]]
    M2.l[[i]] <- pred_models.l$scores[[5]]
    M1.l[[i]] <- pred_models.l$scores[[6]]

  }

  if ((cv == T) | (weights != "LR")){

    cat("output scores ... \n")

    rank_fold  <- order(fold)
    results.df <- data.frame(M12 = unlist(M12.l),
                             M11 = unlist(M11.l),
                             M5 = unlist(M5.l),
                             M3 = unlist(M3.l),
                             M2 = unlist(M2.l),
                             M1 = unlist(M1.l),
                             D_s = D_s[rank_fold],
                             pwm = pwm[rank_fold],
                             PostPr_Shrinkage = PostPr_Shrinkage[rank_fold],
                             PostPr_NoShrinkage = PostPr_NoShrinkage[rank_fold])

    results.l <- list(results = results.df, label = unlist(label.l), names = modelNames)

    save(results.l, file = paste(dir_results, "weights", weights, "_results.l", sep = ''))

    fileName <- paste(dir_results, simpleCap(tf), "_fold", k, "_weights", weights, "_ROC_PRs_all", sep = "")
    result_table <- result_curves(results.l$results, results.l$label, colors = myColors, results.l$names, simpleCap(tf), fileName)

    # write result tables
    fileName <- paste(dir_results, simpleCap(tf),"_fold", k, "_weights", weights, "_scores.txt", sep = "")
    write.table(result_table, file = fileName, sep = "\t", col.names = T, row.names = F, quote = F)

    fileName <- paste(dir_results, simpleCap(tf), "_fold", k, "_weights", weights, "_ROC_PRs_bins", sep = "")
    result_curves(results.l$results[,1:6], results.l$label, colors = myColors[1:6], results.l$names[1:6], simpleCap(tf), fileName)

    fileName <- paste(dir_results, simpleCap(tf), "_fold", k, "_weights", weights, "_ROC_PRs_cpd", sep = "")
    results <- cbind(results.df$PostPr_NoShrinkage, results.df$M5)
    model_name <- c("CENTIPEDE", "M5")
    result_curves(results, results.l$label, colors = c("blue","red"), model_name, simpleCap(tf), fileName)

    if ( weights != "Equal"){

      cat("printing heatmaps ...\n")

      if (chipexo == T){
        labels.df <- data.frame(joint = label_joint)
        label_name <- c("ChIP")
      }else{
        labels.df <- data.frame(macIsaac = label)
        label_name <- c("ChIP")
      }

      n0 <- which(results.l$label == 0)
      n1 <- which(results.l$label == 1)

      rank_all <- c(sample(n0,length(n0)), sample(n1,length(n1)))
      fileName <- paste(dir_results, simpleCap(tf), "_fold", k, "_weights", weights, "_heatmap_data", sep = "")
      image_data(DNaseData[rank_fold,], pwm[rank_fold], cbind(labels.df[rank_fold,]), rank_all, label_name, zMax_DNase = 20, paste("Yeast", simpleCap(tf)), fileName, pngWidth = 3)

      #       rank_subset <- c(sample(n0,length(n1)), sample(n1,length(n1)))
      #       rank_subset <- c(sample(n0,min(1000, length(n0))), sample(n1,length(n1)))

      fileName <- paste(dir_results, simpleCap(tf), "_fold", k, "_weights", weights, "_heatmap_bins", sep = "")
      image_bins(DNaseData[rank_fold,], results.df[,1:6], cbind(labels.df[rank_fold,]), rank_all, modelNames[1:6], label_name, zMax_DNase = 20, paste("Yeast", simpleCap(tf)), fileName, pngWidth = 4)

      results <- cbind(results.df$PostPr_Shrinkage, results.df$M5)
      model_name <- c("CENTIPEDE\nw/ shrinkage", "MILLIPEDE\nM5")
      fileName <- paste(dir_results, simpleCap(tf), "_fold", k, "_weights", weights, "_heatmap_cpd", sep = "")
      image_compare(DNaseData[rank_fold,], results, cbind(labels.df[rank_fold,]), rank_all, model_name, label_name, zMax_DNase = 20, paste("Yeast", simpleCap(tf)), fileName, pngWidth = 4)

    }

  }else if((cv == F) & weights == "LR"){
      # save trained model
      models.l <- pred_models.l$models
      save(models.l, file = paste(dir_results, "LRmodels.l", sep =""))

      # save model fitting summaries
      fileName <- paste(dir_results, simpleCap(tf), "_weights", weights, "_models.txt", sep ="")
      sink(fileName, append = F)
      for ( i in 1:length(models.l)){
        cat(paste("Model", modelNames[i], "\n"))
        print(summary(models.l[[i]]))
        cat("\n")
        print(anova(models.l[[i]]))
        cat("################################################################\n")
      }
      sink()

      # plot model parameters
      fileName <- paste(dir_results, simpleCap(tf), "_weights", weights, "_coefficients", sep = "")
      pdf(file= paste(fileName, ".pdf", sep = ""))
      par(mfrow=c(2,1))
      for ( i in 1:length(models.l)){
      barplot(coef(summary(models.l[[i]]))[-c(1:2),1], cex.names = 0.5, las = 2, ylab ="Regression Coefficients")
      mtext(paste("Model", modelNames[i]), 3)
      barplot(-log10(coef(summary(models.l[[i]]))[-c(1:2),4]+exp(-100)), cex.names = 0.5, las = 2, ylab ="- P value (log10)")
      abline(h = -log10(0.05+exp(-100)), col = "orange")
      abline(h = -log10(0.01+exp(-100)), col = "red")
      }
      dev.off()

  }

  cat("Completed. \n")


}

## fit CENTIPEDE model
centipede_dnase <- function(DNaseData, pwm){
  cat("fit CENTIPEDE \n")
  Xlist <- list(DNase = as.matrix(DNaseData))
  Y <- cbind(rep(1, length(pwm)), pwm)

  centFit_Shrinkage <- fitCentipede(Xlist=Xlist[1],Y=Y,sweeps=MaxSweeps, DampLambda= 0.5,DampNegBin=0.01)
  PostPr_Shrinkage <- centFit_Shrinkage$PostPr

  centFit_NoShrinkage <- fitCentipede(Xlist=Xlist[1],Y=Y,sweeps=MaxSweeps, DampLambda=0.0,DampNegBin=0.0)
  PostPr_NoShrinkage <- centFit_NoShrinkage$PostPr

  cent.df <- data.frame(PostPr_Shrinkage = PostPr_Shrinkage, PostPr_NoShrinkage = PostPr_NoShrinkage)

  return(cent.df)

}

lr_dnase = function(cuts_sum.l, pwm, label, trainIdx, testIdx, weights){
  scores.l = models.l = list()
  weights.l = list(
    M12 = c(-1,rep(1,4), -1, -1, rep(1,4),-1),
    M11 = c(-1,rep(1,4), -1, rep(1,4),-1),
    M5 = c(1, 1, -1, 1, 1),
    M3 = c(1, -1, 1),
    M2 = c(1, 1),
    M1 = c(1),
    M0 = c(1,1))

  for (i in 1: (length(cuts_sum.l))){

    data.df = data.frame(label, pwm, log2(cuts_sum.l[[i]]+0.1))
    models.l[[i]] <- glm(formula = label ~ ., data = data.df[trainIdx,], family = "binomial")
    predictor.m <- as.matrix(cbind(1, data.df[testIdx,-label]))
    if(weights == "LR"){
      predict.logit <- predictor.m %*% models.l[[i]]$coefficients
      scores.l[[i]] <- 1 / (1 + exp(-predict.logit))

    }else if(weights == "Trained"){
      ref_models.l <- get(load(paste0(dir_results, "/", "LRmodels.l")))
      predict.logit <- predictor.m %*% ref_models.l[[i]]$coefficients
      scores.l[[i]] <- 1 / (1 + exp(-predict.logit))

    }else{
      predict.logit <- predictor.m %*% c(0, 1, weights.l[[i]])
      scores.l[[i]] <- predict.logit/max(predict.logit)

    }
  }

  pred.l <- list(scores = scores.l, models = models.l)

  return(pred.l)

}

## 5-fold cross-validation, logistic regression
millipede_analysis(k = 5, weights = "LR")
#
# # get model coefficients
# millipede_analysis(k = 1, weights = "LR")
#
# # equal weights
# millipede_analysis(k = 1, weights = "Equal")
#
# # weights trained from other TFs
# millipede_analysis(k = 1, weights = "Trained")
