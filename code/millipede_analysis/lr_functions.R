
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

  cuts_sum.l = list(M12 = M12,
                    M11 = M11,
                    M5 = M5,
                    M3 = M3,
                    M2 = M2,
                    M1 = M1)
  return(cuts_sum.l)

}

result_curves = function (results.df, label, colors, modelNames, tfName, fileName){
  # browser()
  auc = auprc = sens = prec = array()
  pdf(file= paste(fileName, ".pdf", sep = ""))

  par(xaxs = "i")
  par(yaxs = "i")

  plot(x = c(0,1), y = c(0,1), type = "n", ylab = "Sensitivity", xlab = "False positive rate")
  mtext(paste("ROC curves for", tfName),3,line=1,font=2, cex = 1)
  for (i in 1: ncol(results.df)){
    pred <- prediction(results.df[,i], label)
    perf_ROC <- performance(pred, measure='sens',x.measure='fpr')
    plot(perf_ROC, col = colors[i], lwd=1, add = T)
    auc[i] <- round(performance(pred,"auc")@y.values[[1]],3)
    # auprc[i] <- round(performance(pred,"auprc")@y.values[[1]],3)
    sens[i] <- round(myQ(pred),3)
    prec[i] <- round(myPrec(pred),3)
  }
  legend_auc = paste(modelNames, "auROC =", auc*100, "%")
  legend("bottomright", legend = legend_auc, lty = 1, lwd=2, col= colors, bty="o", cex = 0.6)

  plot(x = c(0,0.01), y = c(0,1), type = "n", ylab = "Sensitivity", xlab = "False positive rate")
  mtext(paste("ROC curves for", tfName),3,line=1,font=2, cex = 1)
  for (i in 1: ncol(results.df)){
    pred <- prediction(results.df[,i], label)
    perf_ROC <- performance(pred,measure='sens',x.measure='fpr')
    plot(perf_ROC, col = colors[i], lwd=1, add = T)
  }
  legend_auc = paste(modelNames, "Sens. at", tFPR*100, "% FPR =", sens*100, "%")
  legend("bottomright", legend = legend_auc, lty = 1, lwd=2, col= colors, bty="o", cex = 0.6)

  # plot(x = c(0,1), y = c(0,1), type = "n", ylab = "Precision", xlab = "Recall")
  # mtext(paste("Precision Recall curves for", tfName), 3,line=1,font=2, cex = 1)
  # for (i in 1: ncol(results.df)){
  #   pred <- prediction(results.df[,i], label)
  #   perf_PR <- performance(pred,measure='prec',x.measure='rec')
  #   plot(perf_PR, col = colors[i], lwd=1, add = T)
  # }
  #
  # legend_auprc = paste(modelNames, "auPRC =", auprc*100, "%")
  # legend("bottomleft", legend = legend_auprc, lty = 1, lwd=2, col= colors, bty="o", cex = 0.6)

  dev.off()

  # result_table = cbind(modelNames, auc, auprc, sens, prec)
  # colnames(result_table) = c("Models","AUROC", "AUPRC", "Sens", "Prec")

  result_table = cbind(modelNames, auc, sens, prec)
  colnames(result_table) = c("Models","AUROC", "Sens", "Prec")
  return(result_table)

}


myPrec <- function(pred){
  P1 <- performance(pred,measure='prec',x.measure='fpr')
  FPR <- attr(P1,"x.values")[[1]]
  Prec <- attr(P1,"y.values")[[1]]
  ii <- min(which(FPR>tFPR))
  #   ii <- min(which(Sens>99))
  max(ii,1)
  Prec[ii]
  ##print(paste(ii,[ii]))
  ##FDR[100]
  ##  browser()
}


myQ <- function(pred){
  P1 <- performance(pred,measure='sens',x.measure='fpr')
  FPR <- attr(P1,"x.values")[[1]]
  Sens <- attr(P1,"y.values")[[1]]
  right <- min(which(FPR>tFPR))
  left <- max(which(FPR<tFPR))
  #   ii <- min(which(Sens>99))
  #   max(ii,1)
  Sens = (Sens[right] + Sens[left])/2
  #   print(paste(ii,[ii]))
  #   FDR[100]
  ##  browser()
}


plotROCs <- function(pred){

  perf <- performance(pred,'sens','fpr')
  par(xaxs = "i")
  par(yaxs = "i")

  plot(perf, colorize=T, lwd=1)

}

plotPRs <- function(pred.l){

  perf <- performance(pred,'prec','rec')
  par(xaxs = "i")
  par(yaxs = "i")

  plot(perf, colorize=T, lwd=1)

}

get_scores = function(pred){
  scores = c(round(performance(pred,"auc")@y.values[[1]],3),
             round(performance(pred,"auprc")@y.values[[1]],3),
             round(myQ(pred),3),
             round(myPrec(pred),3))
}


# filter out motif instances that had more than 20% un-mappable bases in the 200bp window centered
mappability = function(tf_overlap.gr, unmappablebases.gr, ThreshAvgMappability){
  unmappableCounts = countOverlaps(tf_overlap.gr, unmappablebases.gr)
  mappableIdx = which((start(tf_overlap.gr) >0) & (unmappableCounts < width(tf_overlap.gr[1])* (1-ThreshAvgMappability)))
  return(mappableIdx)
}


label_chipexo = function(tf.gr, chipexo.gr){
  overlapIndex.df = as.data.frame(as.matrix(findOverlaps(tf.gr, chipexo.gr, type = "within")))
  elementMetadata(tf.gr)[,"chip_exo"] = 0
  elementMetadata(tf.gr)[,"chip_exo"][overlapIndex.df$query] = 1
  return(tf.gr)
}

label_macisaac = function(tf.gr, MacIsaac.gr){
  strand(MacIsaac.gr) = "*"
  overlapIndex.df = as.data.frame(as.matrix(findOverlaps(tf.gr, MacIsaac.gr, type ="within")))
  elementMetadata(tf.gr)[,"label"] = 0
  elementMetadata(tf.gr)[,"label"][overlapIndex.df$query] = 1
  return(tf.gr)
}


image_data <- function(DNaseData, pwm, labels.df, rank, label_name, zMax_DNase, tfName, fileName, pngWidth){
  #   browser()
  DNaseData[DNaseData > zMax_DNase] = zMax_DNase

  x = c(1:dim(DNaseData[rank,])[2])
  y = c(1:dim(DNaseData[rank,])[1])

  #   myColors <- colorRampPalette(c(rgb(1,1,1),rgb(0,0,1),rgb(0,0,.2)))(16)
  myColors <- colorRampPalette(c(rgb(1,1,1),rgb(0,0,1),rgb(0,0,0.5)))(1000);

  png(file = paste(fileName, ".png", sep = ""), width= pngWidth, height = 3, units='in',res=400)

  colNumber = 1+1+ncol(labels.df)
  #   def.par <- par(no.readonly = TRUE)
  #   pwidth = c(4,rep(0.8,colNumber-1))
  # #   pwidth <- c(0.35,1,1,1,0.5,rep(1,2),3,1.5,0.8,0.8,3);
  #   pwidth <- pwidth/sum(pwidth);
  #   nf <- layout(t(as.matrix(1:length(pwidth))),widths=pwidth)
  #   layout.show(nf)
  layout(matrix(c(1:colNumber), nrow = 1, ncol = colNumber, byrow = F), widths=c(3,1,1), heights = rep(1,colNumber))
  par(oma = c(2,2,2,2))
  cex.text = 0.5

  par(mai=c(0.1,0.1,0.5,0.1))
  image(x, y, z = t(DNaseData[rank,]), col = myColors, axes=FALSE)
  box()
  #   abline(v = flank_lines, lty=2,lwd=1, col = "green")
  #   abline(v = motif_lines,lty=2,lwd=1, col = "red")
  #   axis(1,at=c(0,1),labels=c(-100, 100),padj=-2.1,cex.axis=0.7,tck=-0.05, tick = F)
  #   axis(1, at = c(1, 101, dim(DNaseData)[2]-100, dim(DNaseData)[2]), labels=c(-100, "", "", 100), cex.axis = 0.3)
  #   axis(1,at=c(1, 101, ncol(DNaseData)-100, ncol(DNaseData)),labels=c(-100, "","" ,100),padj=-2.1,cex.axis=0.7,tck=-0.05, tick = T)
  axis(1,at=c(1, 101, ncol(DNaseData)-100, ncol(DNaseData)),labels=c(-100, "","" ,100),padj=-2.1,cex.axis= cex.text,tck=-0.05, tick = T, cex = cex.text)

  mtext('DNase cuts\nin the region',3,line=0.5, cex =cex.text)

  mtext('Dist. to motif (bp)',1, line=1, cex =cex.text)
  #   mtext("Motif sites ordered by ChIP labels",2,line=0.05,padj=-0.2,las=0, cex = 0.6)
  #   text(0.1, par("usr")[4]+0.21, tfName, xpd = NA, pos = 1,cex=1.0, font=2)


  par(mai=c(0.1,0.1,0.5,0.1))
  #     par(mar = c(1,1.5,2,1))
  image(t(as.matrix(pwm[rank])), col = myColors, axes=FALSE)
  box()
  mtext("PWM\nscore",3,line=0.5,cex =cex.text)
  #   text(0, par("usr")[4]+0.01, "PWM score", srt = 0, xpd = NA, pos = 4, cex = 0.6)


  for ( i in 1: ncol(labels.df)){
    par(mai=c(0.1,0.1,0.5,0.1))
    #   par(mar = c(1,1.5,2,1))
    image(t(as.matrix(labels.df[rank,i])), col = c("white","blue"), axes=FALSE)
    box()
    mtext(label_name[i],3,line=0.5, cex =cex.text)
    #     text(0, par("usr")[4]+0.01, label_name[i], srt = 0, xpd = NA, pos = 4, cex = 0.6)
  }

  mtext(tfName, side=3, outer=TRUE, line= -0.5, cex = 0.6)


  dev.off()

  cat("Image printed in ", fileName, "\n")
}

image_data_human <- function(DNaseData, pwm, ConsScore, labels.df, rank, label_name, zMax_DNase, tfName, fileName, pngWidth){
  #   browser()
  DNaseData[DNaseData > zMax_DNase] = zMax_DNase

  x = c(1:dim(DNaseData[rank,])[2])
  y = c(1:dim(DNaseData[rank,])[1])

  #   myColors <- colorRampPalette(c(rgb(1,1,1),rgb(0,0,1),rgb(0,0,.2)))(16)
  myColors <- colorRampPalette(c(rgb(1,1,1),rgb(0,0,1),rgb(0,0,0.5)))(1000);

  png(file = paste(fileName, ".png", sep = ""), width= pngWidth, height = 3, units='in',res=400)

  colNumber = 1+1+1+ncol(labels.df)

  #     def.par <- par(no.readonly = TRUE)
  #     pwidth = c(3,1,1,1)
  #   # #   pwidth <- c(0.35,1,1,1,0.5,rep(1,2),3,1.5,0.8,0.8,3);
  #     pwidth <- pwidth/sum(pwidth);
  #     nf <- layout(t(as.matrix(1:length(pwidth))),widths=pwidth)
  #     layout.show(nf)

  ##
  #   par(mai=c(2,2,2,2))
  #   plot.new();

  layout(matrix(c(1:colNumber), nrow = 1, ncol = colNumber, byrow = F), widths=c(3,1,1,1), heights = rep(1,colNumber))
  par(oma = c(2,2,2,2))
  #   par(cex=1.0)
  cex.text = 0.5

  par(mai=c(0.1,0.1,0.5,0.1))
  image(x, y, z = t(DNaseData[rank,]), col = myColors, axes=FALSE)
  box()
  #   abline(v = flank_lines, lty=2,lwd=1, col = "green")
  #   abline(v = motif_lines,lty=2,lwd=1, col = "red")
  #   axis(1,at=c(0,1),labels=c(-100, 100),padj=-2.1,cex.axis=0.7,tck=-0.05, tick = F)
  #   axis(1, at = c(1, 101, dim(DNaseData)[2]-100, dim(DNaseData)[2]), labels=c(-100, "", "", 100), cex.axis = 0.3)
  axis(1,at=c(1, 101, ncol(DNaseData)-100, ncol(DNaseData)),labels=c(-100, "","" ,100),padj=-2.1,cex.axis=0.7,tck=-0.05, tick = T)

  mtext('DNase cuts',3,line=0.5, cex =cex.text)

  mtext('Dist. to motif (bp)',1, line=1, cex =cex.text)

  #   mtext("Motif sites ordered by ChIP labels",2,line=0.05,padj=-0.2,las=0, cex = 0.6)

  #   text(0.1, par("usr")[4]+0.21, tfName, xpd = NA, pos = 1,cex=1.0, font=2)

  par(mai=c(0.1,0.1,0.5,0.1))
  #     par(mar = c(1,1.5,2,1))
  image(t(as.matrix(pwm[rank])), col = myColors, axes=FALSE)
  box()
  mtext("PWM\nscore",3,line=0.5,cex =cex.text)
  #   text(0, par("usr")[4]+0.01, "PWM score", srt = 0, xpd = NA, pos = 4, cex = 0.6)

  par(mai=c(0.1,0.1,0.5,0.1))
  #     par(mar = c(1,1.5,2,1))
  image(t(as.matrix(ConsScore[rank])), col = myColors, axes=FALSE)
  box()
  mtext("Conservation\nscore",3,line=0.5,cex =cex.text)
  #   text(0, par("usr")[4]+0.01, "PWM score", srt = 0, xpd = NA, pos = 4, cex = 0.6)


  for ( i in 1: ncol(labels.df)){
    par(mai=c(0.1,0.1,0.5,0.1))
    #   par(mar = c(1,1.5,2,1))
    image(t(as.matrix(labels.df[rank,i])), col = c("white","blue"), axes=FALSE)
    box()
    mtext(label_name[i],3,line=0.5, cex =cex.text)
    #     text(0, par("usr")[4]+0.01, label_name[i], srt = 0, xpd = NA, pos = 4, cex = 0.6)
  }

  mtext(tfName, side=3, outer=TRUE, line= -0.5, cex = 0.6)

  dev.off()

  cat("Image printed in ", fileName, "\n")
}

image_bins <- function(DNaseData, results.df, labels.df, rank, model_name, label_name, zMax_DNase, tfName, fileName, pngWidth){

  #   window = c((flank_data - flank_window +1 ) : (ncol(DNaseData) - flank_data + flank_window))
  #   DNaseData = DNaseData[,window]
  #
  #   motif_lines = c((flank+1), (ncol(DNaseData)-flank))
  #   flank1_index = c(1, (flank-buffer))
  #   flank2_index = c((ncol(DNaseData)-flank+1+buffer), ncol(DNaseData))
  #   flank_lines = c(flank1_index, flank2_index)

  DNaseData[DNaseData > zMax_DNase] = zMax_DNase

  x = c(1:dim(DNaseData[rank,])[2])
  y = c(1:dim(DNaseData[rank,])[1])

  #   myColors <- colorRampPalette(c(rgb(1,1,1),rgb(0,0,1),rgb(0,0,.2)))(16)
  myColors <- colorRampPalette(c(rgb(1,1,1),rgb(0,0,1),rgb(0,0,0.5)))(1000);

  png(file = paste(fileName, ".png", sep = ""), width=pngWidth,height = 3, units='in',res=400)

  colNumber = 1+ncol(results.df)+ncol(labels.df)

  layout(matrix(c(1:colNumber), nrow = 1, ncol = colNumber, byrow = F), widths=c(3,rep(1,colNumber-1)), heights = rep(1,colNumber))
  #   par(cex = 1)
  cex.text = 0.5
  par(oma = c(2,2,2,2))

  par(mai=c(0.1,0.05,0.5,0.05))


  image(x, y , z = t(DNaseData[rank,]), col = myColors, axes=FALSE)
  box()
  #   abline(v = flank_lines, lty=2,lwd=1, col = "green")
  #   abline(v = motif_lines,lty=2,lwd=1, col = "red")
  #   axis(1, at = c(1, 101, dim(DNaseData)[2]-100, dim(DNaseData)[2]), labels=c(-100, "", "", 100), cex.axis = 0.3)
  #   axis(1,at=c(1, 101, ncol(DNaseData)-100, ncol(DNaseData)),labels=c(-100, "","" ,100),padj=-2.1,cex.axis=0.7,tck=-0.05, tick = T, cex = 0.5)
  axis(1,at=c(1, 101, ncol(DNaseData)-100, ncol(DNaseData)),labels=c(-100, "","" ,100),padj=-2.1,cex.axis= cex.text,tck=-0.05, tick = T, cex = cex.text)

  mtext('DNase cuts\nin the region',3,line=0.5,cex =cex.text)

  mtext('Dist. to motif (bp)',1,line=1, cex = cex.text)
  #   mtext("Motif sites ordered by ChIP labels",2,line=0.05,padj=-0.2,las=0, cex = 0.6)


  for ( i in 1: ncol(results.df)){
    par(mai=c(0.1,0.05,0.5,0.05))
    #     par(mar = c(1,1.5,2,1))
    image(t(as.matrix(results.df[rank,i])), col = myColors, axes=FALSE)
    box()
    mtext(model_name[i],3,line=0.5,cex = cex.text)
    #     text(-0.9, par("usr")[4]+0.01, model_name[i], srt = 45, xpd = NA, pos = 4, cex = cex.text)

  }

  for ( i in 1: ncol(labels.df)){
    par(mai=c(0.1,0.05,0.5,0.05))
    #   par(mar = c(1,1.5,2,1))
    image(t(as.matrix(labels.df[rank,i])), col = c("white","blue"), axes=FALSE)
    box()
    mtext(label_name[i],3,line=0.5,cex = cex.text)
    #     text(-0.9, par("usr")[4]+0.01, label_name[i], srt = 45, xpd = NA, pos = 4, cex = 0.6)
  }

  mtext(tfName, side=3, outer=TRUE, line= -0.5, cex = 0.6)

  dev.off()

  cat("Image printed in ", fileName, "\n")
}

image_compare <- function(DNaseData, results.df, labels.df, rank, model_name, label_name, zMax_DNase, tfName, fileName, pngWidth){

  #   window = c((flank_data - flank_window +1 ) : (ncol(DNaseData) - flank_data + flank_window))
  #   DNaseData = DNaseData[,window]
  #
  #   motif_lines = c((flank+1), (ncol(DNaseData)-flank))
  #   flank1_index = c(1, (flank-buffer))
  #   flank2_index = c((ncol(DNaseData)-flank+1+buffer), ncol(DNaseData))
  #   flank_lines = c(flank1_index, flank2_index)

  DNaseData[DNaseData > zMax_DNase] = zMax_DNase

  x = c(1:dim(DNaseData[rank,])[2])
  y = c(1:dim(DNaseData[rank,])[1])

  #   myColors <- colorRampPalette(c(rgb(1,1,1),rgb(0,0,1),rgb(0,0,.2)))(16)
  myColors <- colorRampPalette(c(rgb(1,1,1),rgb(0,0,1),rgb(0,0,0.5)))(1000);

  png(file = paste(fileName, ".png", sep = ""), width=pngWidth,height = 4, units='in',res=400)

  colNumber = 1+ncol(results.df)+ncol(labels.df)

  layout(matrix(c(1:colNumber), nrow = 1, ncol = colNumber, byrow = F), widths=c(3,rep(1,colNumber-1)), heights = rep(1,colNumber))
  #   par(cex = 1)
  cex.text = 0.4
  par(oma = c(2,2,2,2))

  par(mai=c(0.1,0.1,0.5,0.1))


  image(x, y , z = t(DNaseData[rank,]), col = myColors, axes=FALSE)
  box()
  #   abline(v = flank_lines, lty=2,lwd=1, col = "green")
  #   abline(v = motif_lines,lty=2,lwd=1, col = "red")
  #   axis(1, at = c(1, 101, dim(DNaseData)[2]-100, dim(DNaseData)[2]), labels=c(-100, "", "", 100), cex.axis = 0.3)
  axis(1,at=c(1, 101, ncol(DNaseData)-100, ncol(DNaseData)),labels=c(-100, "","" ,100),padj=-2.1,cex.axis= cex.text,tck=-0.05, tick = T, cex = cex.text)
  mtext('DNase cuts\nin the region',3,line=0.5,cex =cex.text)

  mtext('Dist. to motif (bp)',1,line=1, cex = cex.text)
  #   mtext("Motif sites ordered by ChIP labels",2,line=0.05,padj=-0.2,las=0, cex = 0.6)


  for ( i in 1: ncol(results.df)){
    par(mai=c(0.1,0.1,0.5,0.1))
    #     par(mar = c(1,1.5,2,1))
    image(t(as.matrix(results.df[rank,i])), col = myColors, axes=FALSE)
    box()
    mtext(model_name[i],3,line=0.5,cex = cex.text)
    #         text(-0.9, par("usr")[4]+0.01, model_name[i], srt = 45, xpd = NA, pos = 4, cex = cex.text)

  }

  for ( i in 1: ncol(labels.df)){
    par(mai=c(0.1,0.1,0.5,0.1))
    #   par(mar = c(1,1.5,2,1))
    image(t(as.matrix(labels.df[rank,i])), col = c("white","blue"), axes=FALSE)
    box()
    mtext(label_name[i],3,line=0.5,cex = cex.text)
    #     text(-0.9, par("usr")[4]+0.01, label_name[i], srt = 45, xpd = NA, pos = 4, cex = 0.6)
  }

  mtext(tfName, side=3, outer=TRUE, line= -0.5, cex = 0.6)

  dev.off()

  cat("Image printed in ", fileName, "\n")
}


# Capitalizing - toupper every first letter of the TF names
simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2), sep="", collapse=" ")
}
