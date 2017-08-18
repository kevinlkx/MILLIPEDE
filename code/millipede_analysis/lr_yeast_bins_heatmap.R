library(GenomicRanges)
library(gplots)
library(ROCR)

source("~/DNaseData/script/lr_functions_heatmap.R")

# Set parameters
args<-commandArgs(trailingOnly=T);
tf = args[1]
# tf = "reb1"

tfName = paste(tf,"_macisaac_score4",sep = "")
flank_data = 100
# flank = 30
# buffer = 5
MaxSweeps = 400
fileType = "pdf"

chipexoList = c("reb1", "rap1", "gal4", "phd1")

dir_CentipedeData_DNase = paste("~/MNaseData/CentipedeData/DNase/", tfName, "/", sep ="")
dir_results = paste('~/DNaseData/Results_forPSB/', tfName, "/", sep ='')
dir.create(dir_results, showWarnings = FALSE, recursive = TRUE)
myColors <- c("pink", "orange", "red", "purple", "magenta", "darkred", "darkgray", "black", "blue", "cyan")

modelNames = c("M12", "M11", "M5", "M3", "M2", "M1", "M1 w/o PWM", "PWM only",
               "CENTIPEDE w/ shrinkage", "CENTIPEDE w/o shrinkage")

ThreshAvgMappability <- 0.8;

tFPR <- 1.0/100

analysis = function(k, weights){
#   browser()
  if(sum(tf == chipexoList)){
    chipexo = T
  }else{
    chipexo = F
  }

  if(k > 1) {
    cv = T
  }else{
    cv = F
  }

  cat("Analyzing:",tfName, "fold", k, "weights", weights, "\n")
#   browser()
  CentipedeData_DNase.l = get(load(paste(dir_CentipedeData_DNase, "dnase_cuts_", flank_data, "bp.l", sep = "")))

  tf_overlap.gr = CentipedeData_DNase.l$tf.gr

  # filter sites by mappability (remove the sites with more than > 20% unmappable bases) and start position > 0
  cat(length(tf_overlap.gr), " - start position > 0 and mappable -> ")
  select_idx = mappability(tf_overlap.gr)
  tf_overlap.gr = tf_overlap.gr[select_idx]
  cat(length(tf_overlap.gr),"\n")

  tf_macIsaac.gr = get(load(paste("~/MNaseData/tf_coord_updated/", tf,".gr", sep = "")))

  tf.gr = tf_overlap.gr
  start(tf.gr) = start(tf_overlap.gr) + flank_data
  end(tf.gr) = end(tf_overlap.gr) - flank_data

  rank = order(as.integer(seqnames(tf.gr)), start(tf.gr))

  DNaseData = CentipedeData_DNase.l$counts[select_idx,][rank,]
  pwm = elementMetadata(tf.gr)[,"score"][rank]
  tf.gr = tf.gr[rank]

  if (chipexo == T){
    tf_chipexo.gr = get(load(paste("~/MNaseData/tf_coord_updated/", tf,"_chip_exo.gr", sep = "")))

    tf.gr = label_macisaac(tf.gr, tf_macIsaac.gr)
    tf.gr = label_chipexo(tf.gr, tf_chipexo.gr)
    elementMetadata(tf.gr)[,"label_joint"] = as.integer(elementMetadata(tf.gr)[,"label"] | elementMetadata(tf.gr)[,"chip_exo"])
    label = elementMetadata(tf.gr)[,"label"][rank]
    label_chip_exo = elementMetadata(tf.gr)[,"chip_exo"][rank]
    label_joint = elementMetadata(tf.gr)[,"label_joint"][rank]
    cat("macIsaac sites number:", length(tf_macIsaac.gr),
        "chip-exo sites number:", length(tf_chipexo.gr), "\n")
    cat("matched macIsaac sites number:", sum(elementMetadata(tf.gr)[,"label"] == 1),
        "matched chip-exo sites number:", sum(elementMetadata(tf.gr)[,"chip_exo"] == 1),
        "Joint positive sites number:", sum(elementMetadata(tf.gr)[,"label_joint"] == 1),"\n")

  }else if(chipexo == F){
    tf.gr = label_macisaac(tf.gr, tf_macIsaac.gr)
    elementMetadata(tf.gr)[,"label_joint"] = as.integer(elementMetadata(tf.gr)[,"label"])
    label_joint = label = elementMetadata(tf.gr)[,"label"][rank]
    cat("macIsaac sites number:", length(tf_macIsaac.gr), "\n")
    cat("matched macIsaac sites number:", sum(elementMetadata(tf.gr)[,"label"] == 1),"\n")

  }


  cuts_sum.l = sum_cuts(DNaseData)
  D_s = as.numeric(cuts_sum.l$M1)

  pred_pwm <- prediction(pwm, label_joint)
  pred_D_s <- prediction(D_s, label_joint)

  cent.df = centipede_dnase(DNaseData, pwm)
  PostPr_Shrinkage = cent.df$PostPr_Shrinkage
  PostPr_NoShrinkage = cent.df$PostPr_NoShrinkage

  pred_centipede_Shrinkage <- prediction(PostPr_Shrinkage, label_joint)
  pred_centipede_NoShrinkage <- prediction(PostPr_NoShrinkage, label_joint)

  ####################################################################################

  n1 = length(rank[label_joint == 1])
  n0 = length(rank[label_joint == 0])
  fold1 = rep(1:k,length.out = n1)[sample(n1,n1)]
  fold0 = rep(1:k,length.out = n0)[sample(n0,n0)]
  fold = c(fold1,fold0)

  pred_models.l = label.l = list()
  M12.l = M11.l = M5.l = M3.l = M2.l = M1.l = list()

  for ( i in 1:k ){

    testIdx  = which(fold == i)
    if(cv == T){
      trainIdx = which(fold != i)
    }else{
      trainIdx = testIdx
    }

    cat(i,"cv \n")
#     browser()
    pred_models.l = lr_dnase(cuts_sum.l, pwm, label_joint, trainIdx, testIdx, weights)
    label.l[[i]] = label_joint[testIdx]
    ################
    M12.l[[i]] = pred_models.l$scores[[1]]
    M11.l[[i]] = pred_models.l$scores[[2]]
    M5.l[[i]] = pred_models.l$scores[[3]]
    M3.l[[i]] = pred_models.l$scores[[4]]
    M2.l[[i]] = pred_models.l$scores[[5]]
    M1.l[[i]] = pred_models.l$scores[[6]]

  }
#   browser()
  if ((cv == T) | (weights != "LR")){

    cat("output scores ... \n")

    rank_fold  = order(fold)
    results.df = data.frame(
      M12 = unlist(M12.l),
      M11 = unlist(M11.l),
      M5 = unlist(M5.l),
      M3 = unlist(M3.l),
      M2 = unlist(M2.l),
      M1 = unlist(M1.l),
      D_s = D_s[rank_fold],
      pwm = pwm[rank_fold],
      PostPr_Shrinkage = PostPr_Shrinkage[rank_fold],
      PostPr_NoShrinkage = PostPr_NoShrinkage[rank_fold])

    results.l = list(results = results.df, label = unlist(label.l), names = modelNames)

#     save(results.l, file = paste(dir_results, "weights", weights, "_results.l", sep = ''))

#     fileName = paste(dir_results, simpleCap(tf), "_fold", k, "_weights", weights, "_ROC_PRs_all", sep = "")
#     result_table = result_curves(results.l$results, results.l$label, colors = myColors, results.l$names, simpleCap(tf), fileName)
#
#     # write result tables
#     fileName = paste(dir_results, simpleCap(tf),"_fold", k, "_weights", weights, "_scores.txt", sep = "")
#     write.table(result_table, file = fileName, sep = "\t", col.names = T, row.names = F, quote = F)
#
#     fileName = paste(dir_results, simpleCap(tf), "_fold", k, "_weights", weights, "_ROC_PRs_bins", sep = "")
#     result_curves(results.l$results[,1:6], results.l$label, colors = myColors[1:6], results.l$names[1:6], simpleCap(tf), fileName)
#
#     fileName = paste(dir_results, simpleCap(tf), "_fold", k, "_weights", weights, "_ROC_PRs_cpd", sep = "")
#     results = cbind(results.df$PostPr_NoShrinkage, results.df$M5)
#     model_name = c("CENTIPEDE", "M5")
#     result_curves(results, results.l$label, colors = c("blue","red"), model_name, simpleCap(tf), fileName)

    if ( weights != "Equal"){
      cat("printing heatmaps ...\n")

      if (chipexo == T){
        labels.df = data.frame(joint = label_joint)
        label_name = c("ChIP\nlabel")
      }else{
        labels.df = data.frame(macIsaac = label)
        label_name = c("ChIP\nlabel")
      }

      n0 = which(results.l$label == 0)
      n1 = which(results.l$label == 1)

      rank_all = c(sample(n0,length(n0)), sample(n1,length(n1)))
      fileName = paste(dir_results, simpleCap(tf), "_fold", k, "_weights", weights, "_heatmap_data", sep = "")
      image_data(DNaseData[rank_fold,], pwm[rank_fold], cbind(labels.df[rank_fold,]), rank_all, label_name, zMax_DNase = 20, paste("Yeast", simpleCap(tf)), fileName, pngWidth = 3)

#       rank_subset = c(sample(n0,length(n1)), sample(n1,length(n1)))
#       rank_subset = c(sample(n0,min(1000, length(n0))), sample(n1,length(n1)))

      fileName = paste(dir_results, simpleCap(tf), "_fold", k, "_weights", weights, "_heatmap_bins", sep = "")
      image_bins(DNaseData[rank_fold,], results.df[,1:6], cbind(labels.df[rank_fold,]), rank_all, modelNames[1:6], label_name, zMax_DNase = 20, paste("Yeast", simpleCap(tf)), fileName, pngWidth = 4)

      results = cbind(results.df$PostPr_Shrinkage, results.df$M5)
      model_name = c("CENTIPEDE\nw/ shrinkage", "MILLIPEDE\nM5")
      fileName = paste(dir_results, simpleCap(tf), "_fold", k, "_weights", weights, "_heatmap_cpd", sep = "")
      image_compare(DNaseData[rank_fold,], results, cbind(labels.df[rank_fold,]), rank_all, model_name, label_name, zMax_DNase = 20, paste("Yeast", simpleCap(tf)), fileName, pngWidth = 4)
    }

  }else if((cv == F) & weights == "LR"){

      models.l = pred_models.l$models
#       save(models.l, file = paste(dir_results, "LRmodels.l", sep =""))

      fileName = paste(dir_results, simpleCap(tf), "_weights", weights, "_models.txt", sep ="")
      sink(fileName, append = F)
      for ( i in 1:length(models.l)){
        cat(paste("Model", modelNames[i], "\n"))
        print(summary(models.l[[i]]))
        cat("\n")
        print(anova(models.l[[i]]))
        cat("################################################################\n")
      }
      sink()
#     browser()
      fileName = paste(dir_results, simpleCap(tf), "_weights", weights, "_coefficients", sep = "")
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

centipede_dnase = function(DNaseData, pwm){

#   window = c((flank_data - flank +1 ) : (ncol(DNaseData) - flank_data + flank))
#   DNaseData = DNaseData[,window]

  Xlist = list(DNase = as.matrix(DNaseData))
  Y=cbind(rep(1, length(pwm)), pwm)

  source("~/MNaseData/R_script/fitCentipede_footprint.R")
  centFit_Shrinkage <- fitCentipede(Xlist=Xlist[1],Y=Y,sweeps=MaxSweeps, DampLambda= 0.5,DampNegBin=0.01)
  PostPr_Shrinkage = centFit_Shrinkage$PostPr

  centFit_NoShrinkage <- fitCentipede(Xlist=Xlist[1],Y=Y,sweeps=MaxSweeps, DampLambda=0.0,DampNegBin=0.0)

  PostPr_NoShrinkage = centFit_NoShrinkage$PostPr

  cent.df = data.frame(PostPr_Shrinkage = PostPr_Shrinkage, PostPr_NoShrinkage = PostPr_NoShrinkage)

  return(cent.df)

}

sum_cuts = function(DNaseData){

  flank_L1 = 81:100
  flank_L2 = flank_L1 - 20
  flank_L3 = flank_L2 - 20
  flank_L4 = flank_L3 - 20
  flank_L5 = flank_L4 - 20
  flank_R1 = (ncol(DNaseData)-99): (ncol(DNaseData)-80)
  flank_R2 = flank_R1 + 20
  flank_R3 = flank_R2 + 20
  flank_R4 = flank_R3 + 20
  flank_R5 = flank_R4 + 20

  motif_L = 101: round(ncol(DNaseData)/2)
  motif_R = (round(ncol(DNaseData)/2)+1) : (ncol(DNaseData)-100)

  M12 = data.frame(
    flank_L5 = rowSums(DNaseData[,flank_L5]),
    flank_L4 = rowSums(DNaseData[,flank_L4]),

    flank_L3 = rowSums(DNaseData[,flank_L3]),
    flank_L2 = rowSums(DNaseData[,flank_L2]),
    flank_L1 = rowSums(DNaseData[,flank_L1]),

    motif_L = rowSums(DNaseData[,motif_L]),
    motif_R = rowSums(DNaseData[,motif_R]),

    flank_R1 = rowSums(DNaseData[,flank_R1]),
    flank_R2 = rowSums(DNaseData[,flank_R2]),
    flank_R3 = rowSums(DNaseData[,flank_R3]),

    flank_R4 = rowSums(DNaseData[,flank_R4]),
    flank_R5 = rowSums(DNaseData[,flank_R5]))

  left2_sum = rowSums(M12[,c(1,2)])
  left1_sum = rowSums(M12[,c(3,4,5)])
  motif_sum = rowSums(M12[,c(6, 7)])
  right1_sum = rowSums(M12[,c(8,9, 10)])
  right2_sum = rowSums(M12[,c(11,12)])


  M11 = data.frame(M12[,c(1:5)], motif_sum, M12[,c(8:12)])

  M5 = data.frame(left2_sum, left1_sum, motif_sum, right1_sum, right2_sum)

  M3 = data.frame(left1_sum, motif_sum, right1_sum)

  M2 = data.frame(left1_sum, right1_sum)

  M1 = left1_sum + right1_sum

#   left_diff = left1_sum - motif_sum/motif_width*20
#   right_diff = right1_sum - motif_sum/motif_width*20
#   left_diff[left_diff<0] = 0
#   right_diff[right_diff<0] = 0
#
#   M0 = data.frame(left_diff, right_diff)


  cuts_sum.l = list(M12 = M12,
                M11 = M11,
                M5 = M5,
                M3 = M3,
                M2 = M2,
                M1 = M1)
  return(cuts_sum.l)

}

lr_dnase = function(cuts_sum.l, pwm, label, trainIdx, testIdx, weights){
#   browser()
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
    predictor.m = as.matrix(cbind(1, data.df[testIdx,-label]))
    if(weights == "LR"){
      predict.logit = predictor.m %*% models.l[[i]]$coefficients
      scores.l[[i]] = 1 / (1 + exp(-predict.logit))

    }else if(weights == "Trained"){
      ref_models.l = get(load("~/DNaseData/Results_forPSB/reb1_macisaac_score4/LRmodels.l"))
      predict.logit = predictor.m %*% ref_models.l[[i]]$coefficients
      scores.l[[i]] = 1 / (1 + exp(-predict.logit))

    }else{
      predict.logit = predictor.m %*% c(0, 1, weights.l[[i]])
      scores.l[[i]] = predict.logit/max(predict.logit)

    }
  }

  pred.l = list(scores = scores.l, models = models.l)

  return(pred.l)

}


analysis(k = 5, weights = "LR")

# get model coefficients
# analysis(k = 1, weights = "LR")
#
# # equal weights
# analysis(k = 1, weights = "Equal")

# weights trained from other TFs
# analysis(k = 1, weights = "Trained")
