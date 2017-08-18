rm(list = ls())
require(CENTIPEDE)
# require(GenomicRanges)
require(ROCR)
source("~/DNaseData/script/lr_functions_heatmap.R")

##
args<-commandArgs(trailingOnly=T);
tfName = args[1]
goldstandard = args[2]

# tfName = "REST"
# goldstandard = "peak"
# goldstandard = "peak"
# k = 5
# weights = T

if (tfName == "CTCF"){
  PwmName <- "M01200"   ## CTCF 
  ChIPseqName <- 'BroadGm12878Ctcf'
}else if(tfName == "NRSF" | tfName == "REST"){
  PwmName <- "M00256"  ## NRSF
  ChIPseqName <- 'HudsonalphaGm12878Nrsf'
}else if(tfName == "SRF"){
  PwmName <- "MA0083"   ## Srf
  ChIPseqName <- 'HudsonalphaGm12878Srf'
}else if(tfName == "GABP" | tfName == "GABPA"){
  PwmName <- "M00108"   ## GABP
  ChIPseqName <- 'HudsonalphaGm12878Gabp'
}else if(tfName == "JUND"){
  PwmName <- "M00199"   ## Jund
  ChIPseqName <- 'YaleGm12878Jund'
}else if(tfName == "MAX"){
  PwmName <- "M00118"   ## Max
  ChIPseqName <- 'YaleGm12878Max'
}

CellName <- "BcellAll"
Treatment="DNase2"
dir_results = paste('~/DNaseData/Results_forPSB/', tfName, "/", goldstandard, "/", sep ='')
dir.create(dir_results, showWarnings = FALSE, recursive = TRUE)
# myColors <- c("yellow", "pink", "orange", "palevioletred", "magenta", "red", "darkred", "darkgray", "black", "blue", "cyan")
myColors <- c("yellow", "pink", "orange", "red", "purple", "magenta", "darkred", "darkgray", "black", "blue", "cyan")

modelNames = c("M24", "M12", "M11", "M5", "M3", "M2","M1", "M1 w/o PWM", "PWM only", 
               "CENTIPEDE w/ shrinkage", "CENTIPEDE w/o shrinkage")
ThreshAvgMappability <- 0.8;
Window <- 100;
NewWindow <- 100;
MaxSweeps <- 400;

tFPR <- 1.0/100


analysis = function(k, weights){
  
  if(k > 1) {
    cv = T
  }else{
    cv = F
  }
  
  cat("Analyzing:",tfName, PwmName,"on",CellName,"using",Treatment,"Window",NewWindow, "fold", k, "weights", weights, "\n");
  
  ## Find name of the PWM
  aa<-read.table('~/DNaseData/CentipedeFiguresV2/data/PwmIDandNames.txt.gz',header=F,stringsAsFactors=F)
  Tname<-aa[match(PwmName,aa$V1),2];
  if (is.na(Tname)){
    Tname <- PwmName
  }
  print(Tname);
  rm(aa);
  ##
  
  CutSiteFolder <- paste('~/DNaseData/CentipedeFiguresV2/data/',CellName,'/',Treatment,'/',sep="")
  print(CutSiteFolder)
  
  ## File definition
  AnnotationFile <- paste('~/DNaseData/CentipedeFiguresV2/data/Annotated/',PwmName,'.bed.gz',sep="")
  CutSitesFile<- paste(CutSiteFolder,PwmName,'.bed.CutSites.',Window,'.txt.gz',sep="")
  MappabilityFile <- paste('~/DNaseData/CentipedeFiguresV2/data/mappability/',PwmName,'.bed.CutSites.',100,'.txt.gz',sep="")
  SeqFile <- paste('../data/Sequences/',PwmName,'.seq.gz',sep="")
  
  ## Opening files 
  AnnoTable <- read.table(gzfile(AnnotationFile),header=F, stringsAsFactors=F)[,1:5];
  names(AnnoTable)<-c("Chrom","ChromStart","ChromEnd","Strand","PwmScore")
  ##names(AnnoTable) <- scan('/data/share/AllPwmSites/Annotated/ColNames.txt',character());
  Mlen <- AnnoTable$ChromEnd[1]-AnnoTable$ChromStart[1]+1
  
  ## Mappability 
  Mappability <- read.table(gzfile(MappabilityFile),header=F, sep="\t",stringsAsFactors=F)
  S2<-dim(Mappability)[2]
  Mappability<-Mappability[,2:S2];
  ##Mappability[Mappability>1]<-0
  Mappability <- apply(Mappability, 2, function(i){i[i > 1] <- 0; return(i)})
  
  AvgMappability <- rowMeans(Mappability);
  IndCenter <- ((Window-10):(Window+Mlen+10))
  IndCenter <- c(IndCenter,IndCenter+S2/2+1)
  AvgMappAtCenter <- rowMeans(Mappability[,IndCenter]);
  IndMappable <- ((AvgMappability>0.8));## & (AvgMappAtCenter>0.8))
  
  Mappability <- Mappability[IndMappable,]
  AnnoTable <- AnnoTable[IndMappable,]
  
  ## ## PWM
  pwm = AnnoTable$PwmScore
  
  ## ## TSS
  TssAnno <- read.table('~/DNaseData/CentipedeFiguresV2/data/Ensembl.txt',header=F,as.is=T)
  names(TssAnno) <- c("Chr","Pos","Strand","Id")
  TssAnno$Chr <- factor(TssAnno$Chr)
  TssAnno <- split(TssAnno,TssAnno$Chr)
  
  findClosestTss <- function(Chr,Center,Tss){
    min(abs(Tss[[Chr]]$Pos-Center))
  }
  TssDistance <- mapply(findClosestTss,AnnoTable$Chrom,(AnnoTable$ChromStart+AnnoTable$ChromEnd)/2,MoreArgs=list(Tss=TssAnno))
  TssScore <- as.numeric(1/(1+abs(TssDistance/1000)))
  
  ## Conservation
  ConsFile <- paste('~/DNaseData/CentipedeFiguresV2/data/PhastCons/',PwmName,'.bed.phastCons44way.placental.motif_only.txt.gz',sep="")
  ConsScore <- read.table(gzfile(ConsFile),header=F, stringsAsFactors=F);
  ConsScore <- ConsScore$V4*ConsScore$V5
  ConsScore <- ConsScore[IndMappable]
  
  ## Read DNase cut-sites
  CutSitesDnase <- read.table(gzfile(CutSitesFile),header=F, sep="\t",stringsAsFactors=F)
  S<-dim(CutSitesDnase)[2]
  CutSitesDnase<-CutSitesDnase[,2:S];
  S<-dim(CutSitesDnase)[2]
  Window <- (S/2-Mlen)/2;
  
  ## Readjust Window length for Dnase CutSites ## Implement as a function!!
  NewWindowInd <- rep(0,S)
  NewWindowInd[(Window-NewWindow+1):(S/2-(Window-NewWindow))]<-1
  NewWindowInd[(S/2+1+Window-NewWindow):(S-(Window-NewWindow))] <- 1
  NewWindowInd <- NewWindowInd==1
  Window <- NewWindow;
  
  CutSitesDnase <- CutSitesDnase[IndMappable,NewWindowInd]
  ## Clip extreme values
  ## CutSitesDnase[CutSitesDnase>20] <- 20;
#   CutSitesDnase <- apply(CutSitesDnase, 2, function(i){i[i > 20] <- 20; return(i)})
  
  S <- dim(CutSitesDnase)[2]
  
  ## load Centipede defined ChIP labels
  chip_seq = read.table(
    paste('~/DNaseData/CentipedeGoldStandard/Results_',ChIPseqName, "_", PwmName, "_", CellName, '.Centipede.bed',sep=""),
                        header = T, na.strings = "NA")
  
  IdxPos = which(chip_seq$PeakCall == 1)
  if (goldstandard == "site") {
    IdxNeg = which(chip_seq$PeakCall == 0)
  }else{
    IdxNeg = which(chip_seq$PeakCall == 0 | is.na(chip_seq$PeakCall))
  }
  
  IdxAll <- c(IdxPos,IdxNeg)
  label = c(rep(1,length(IdxPos)),rep(0,length(IdxNeg)))
  # label_chip = chip_seq$PeakCall[IdxAll]
  
  cat("negative sites:", length(IdxNeg), "positive sites:", length(IdxPos), "\n");
  
  
  DNaseData = CutSitesDnase[IdxAll,]
  pwm = pwm[IdxAll]
  ConsScore = ConsScore[IdxAll]
  TssScore = TssScore[IdxAll]
    
  cuts_sum.l = sum_cuts(DNaseData)
  D_s = as.numeric(cuts_sum.l$M1)
  
#   window_cuts.df = sum_cuts(CutSitesDnase, buffer)$window_cuts.df
  
  Xlist <- c(list(Dnase=DNaseData));
  Xlist <- lapply(Xlist,as.matrix);
  
  Y <- data.frame(IntCept=1, Pwm = pwm,Cons=ConsScore);
#   Y <- data.frame(IntCept=1, Pwm = pwm,Cons=ConsScore,Tss=TssScore);
  
  Y <- as.matrix(Y)
  
  source("~/MNaseData/R_script/fitCentipede_footprint.R")
  centFit_Shrinkage <- fitCentipede(Xlist=Xlist[1],Y=Y,sweeps=MaxSweeps, DampLambda=0.5,DampNegBin=0.01)  
  cat("correlation Shrinkage vs. No Shrinkage: ", cor(centFit_Shrinkage$PostPr, chip_seq$PostPr[IdxAll]), "\n")
  
  centFit_NoShrinkage <- fitCentipede(Xlist=Xlist[1],Y=Y,sweeps=MaxSweeps, DampLambda=0.0,DampNegBin=0.0)  
  cat("correlation Shrinkage vs. No Shrinkage: ", cor(centFit_Shrinkage$PostPr, centFit_NoShrinkage$PostPr), "\n")
  
  PostPr_Shrinkage = centFit_Shrinkage$PostPr
  PostPr_NoShrinkage = centFit_NoShrinkage$PostPr
  
  pred_pwm <- prediction(pwm, label)
  pred_D_s <- prediction(D_s, label)
  pred_centipede_Shrinkage <- prediction(PostPr_Shrinkage, label)
  pred_centipede_NoShrinkage <- prediction(PostPr_NoShrinkage, label)
  
  
####################################################################################
  
  n1 = sum(label == 1)
  n0 = sum(label == 0)
  fold1 = rep(1:k,length.out = n1)[sample(n1,n1)] 
  fold0 = rep(1:k,length.out = n0)[sample(n0,n0)] 
  fold = c(fold1,fold0)
  
  pred_models.l = label.l = list()
  M24.l = M12.l = M11.l = M5.l = M3.l = M2.l = M1.l = list()
    
  for ( i in 1:k ){
    
    testIdx  = (fold == i)
    if(cv == T){
      trainIdx = (fold != i)
    }else{
      trainIdx = testIdx
    }
    cat(i,"\n")
      
#     browser()
            
    pred_models.l = lr_dnase(cuts_sum.l, pwm, ConsScore, TssScore, label, trainIdx, testIdx, weights)
    label.l[[i]] = label[testIdx]
    ################
    M24.l[[i]] = pred_models.l$scores[[1]]
    M12.l[[i]] = pred_models.l$scores[[2]]
    M11.l[[i]] = pred_models.l$scores[[3]]
    M5.l[[i]] = pred_models.l$scores[[4]]
    M3.l[[i]] = pred_models.l$scores[[5]]
    M2.l[[i]] = pred_models.l$scores[[6]]
    M1.l[[i]] = pred_models.l$scores[[7]]
  }
  #     browser()
  if ((cv == T) | (weights != "LR")){
#     browser()
    cat("output scores ...\n")
    
    rank_fold  = order(fold)
    results.df = data.frame(
      M24 = unlist(M24.l),
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
#     save(results.l, file = paste(dir_results, "weights", weights, "results.l", sep = ''))
#     
#     fileName = paste(dir_results, "fold", k, "_weights", weights, "_ROC_PRs_all", sep = "")
#     result_table = result_curves(results.l$results, results.l$label, colors = myColors, results.l$names, tfName, fileName)
#     
#     # write result tables
#     fileName = paste(dir_results, "fold", k, "_weights", weights, "_scores.txt", sep = "")
#     write.table(result_table, file = fileName, sep = "\t", col.names = T, row.names = F, quote = F)
#     
#     fileName = paste(dir_results, "fold", k, "_weights", weights, "_ROC_PRs_bins", sep = "")
#     result_curves(results.l$results[,1:7], results.l$label, colors = myColors[1:7], results.l$names[1:7], tfName, fileName)
#     
#     fileName = paste(dir_results, "fold", k, "_weights", weights, "_ROC_PRs_cpd", sep = "")
#     results = cbind(results.df$PostPr_Shrinkage, results.df$M5)
#     model_name = c("CENTIPEDE", "M5")
#     result_curves(results, results.l$label, colors = c("blue","red"), model_name, tfName, fileName)
    
    if ( weights != "Equal"){
      cat("printing heatmaps ... \n")
      
      DNaseData = DNaseData[1:(ncol(DNaseData)/2)] + DNaseData[ncol(DNaseData):(ncol(DNaseData)/2+1)]
      label_name = c("ChIP\nlabel")
      labels.df = data.frame(label = label)
      
      n0 = which(results.l$label == 0)
      n1 = which(results.l$label == 1)
#       browser()
      rank_all = c(sample(n0,length(n0)), sample(n1,length(n1)))
      fileName = paste(dir_results, "fold", k, "_weights", weights, "_heatmap_data", sep = "")
      image_data(DNaseData[rank_fold,], pwm[rank_fold], cbind(labels.df[rank_fold,]), rank_all, label_name, zMax_DNase = 7, paste("Human", "REST/NRSF"), fileName, pngWidth = 3)
      
#       rank_subset = c(sample(n0,min(1000, length(n0))), sample(n1,length(n1)))
      fileName = paste(dir_results, "fold", k, "_weights", weights, "_heatmap_bins", sep = "")
      image_bins(DNaseData[rank_fold,], results.df[,1:7], cbind(labels.df[rank_fold,]), rank_all, modelNames[1:7], label_name, zMax_DNase = 7, paste("Human", tfName), fileName, pngWidth = 4)
      
      results = cbind(results.df$PostPr_Shrinkage, results.df$M5)
      model_name = c("CENTIPEDE\nw/ shrinkage", "MILLIPEDE\nM5")
      fileName = paste(dir_results, "fold", k, "_weights", weights, "_heatmap_cpd", sep = "")
      image_compare(DNaseData[rank_fold,], results, cbind(labels.df[rank_fold,]), rank_all, model_name, label_name, zMax_DNase = 7, paste("Human", tfName), fileName, pngWidth = 4)
    }
    
    
  }else if((cv == F) & weights == "LR" & goldstandard == "site"){

      models.l = pred_models.l$models
#       save(models.l, file = paste(dir_results, "LRmodels.l", sep =""))
      
      fileName = paste(dir_results, "weights", weights, "_models.txt", sep ="")
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
      fileName = paste(dir_results, "weights", weights, "_coefficients", sep = "")
      pdf(file= paste(fileName, ".pdf", sep = "")) 
      par(mfrow=c(2,1))
      for ( i in 1:length(models.l)){
        barplot(coef(summary(models.l[[i]]))[-c(1:3),1], cex.names = 0.5, las = 2, ylab ="Regression Coefficients")
        mtext(paste("Model", modelNames[i]), 3)
        barplot(-log10(coef(summary(models.l[[i]]))[-c(1:3),4]+exp(-100)), cex.names = 0.5, las = 2, ylab ="- P value (log10)")
        abline(h = -log10(0.05+exp(-100)), col = "orange")
        abline(h = -log10(0.01+exp(-100)), col = "red")
        
      }
      dev.off()
    
  }

  
  cat("Completed. \n")
  
}


sum_cuts = function(DNaseData){
    
  flank_LF1 = 81:100
  flank_LF2 = flank_LF1 - 20
  flank_LF3 = flank_LF2 - 20
  flank_LF4 = flank_LF3 - 20
  flank_LF5 = flank_LF4 - 20
  flank_RF1 = (ncol(DNaseData)/2-99): (ncol(DNaseData)/2-80)
  flank_RF2 = flank_RF1 + 20
  flank_RF3 = flank_RF2 + 20
  flank_RF4 = flank_RF3 + 20
  flank_RF5 = flank_RF4 + 20
  
  flank_LR1 = ncol(DNaseData)/2 + flank_LF1
  flank_LR2 = ncol(DNaseData)/2 + flank_LF2
  flank_LR3 = ncol(DNaseData)/2 + flank_LF3
  flank_LR4 = ncol(DNaseData)/2 + flank_LF4
  flank_LR5 = ncol(DNaseData)/2 + flank_LF5
  flank_RR1 = ncol(DNaseData)/2 + flank_RF1
  flank_RR2 = ncol(DNaseData)/2 + flank_RF2
  flank_RR3 = ncol(DNaseData)/2 + flank_RF3
  flank_RR4 = ncol(DNaseData)/2 + flank_RF4
  flank_RR5 = ncol(DNaseData)/2 + flank_RF5
  
  motif_LF = 101: round(ncol(DNaseData)/4)
  motif_RF = (round(ncol(DNaseData)/4)+1) : (ncol(DNaseData)/2-100)
  motif_LR = ncol(DNaseData)/2 + motif_LF
  motif_RR = ncol(DNaseData)/2 + motif_RF
  
  M24 = data.frame(
    flank_LF5 = rowSums(DNaseData[,flank_LF5]),
    flank_LF4 = rowSums(DNaseData[,flank_LF4]),
    flank_LF3 = rowSums(DNaseData[,flank_LF3]),
    flank_LF2 = rowSums(DNaseData[,flank_LF2]),
    flank_LF1 = rowSums(DNaseData[,flank_LF1]), 
    
    motif_LF = rowSums(DNaseData[,motif_LF]),
    motif_RF = rowSums(DNaseData[,motif_RF]),
    
    flank_RF1 = rowSums(DNaseData[,flank_RF1]),
    flank_RF2 = rowSums(DNaseData[,flank_RF2]),
    flank_RF3 = rowSums(DNaseData[,flank_RF3]),
    flank_RF4 = rowSums(DNaseData[,flank_RF4]),
    flank_RF5 = rowSums(DNaseData[,flank_RF5]),

    flank_LR5 = rowSums(DNaseData[,flank_LR5]),
    flank_LR4 = rowSums(DNaseData[,flank_LR4]),
    flank_LR3 = rowSums(DNaseData[,flank_LR3]),
    flank_LR2 = rowSums(DNaseData[,flank_LR2]),
    flank_LR1 = rowSums(DNaseData[,flank_LR1]),
    
    motif_LR = rowSums(DNaseData[,motif_LR]),
    motif_RR = rowSums(DNaseData[,motif_RR]),
    
    flank_RR1 = rowSums(DNaseData[,flank_RR1]),
    flank_RR2 = rowSums(DNaseData[,flank_RR2]),
    flank_RR3 = rowSums(DNaseData[,flank_RR3]),
    flank_RR4 = rowSums(DNaseData[,flank_RR4]),
    flank_RR5 = rowSums(DNaseData[,flank_RR5]))
  
  # Roger Pigue-Regi did case 2 of counts but showed incorrect combination in figure 1
  M12 = M24[,c(1:12)] + M24[,c(24:13)]
  
  left2_sum = rowSums(M12[,c(1,2)])
  left1_sum = rowSums(M12[,c(3,4,5)])
  motif_sum = rowSums(M12[,c(6, 7)])
  right1_sum = rowSums(M12[,c(8,9,10)])
  right2_sum = rowSums(M12[,c(11,12)])
  
  M11 = data.frame(M12[,c(1:5)], motif_sum, M12[,c(8:12)])
  
  M5 = data.frame(left2_sum, left1_sum, motif_sum, right1_sum, right2_sum)
  
  M3 = data.frame(left1_sum, motif_sum, right1_sum)
  
  M2 = data.frame(left1_sum, right1_sum)
  
  M1 = left1_sum + right1_sum
  
  
  cuts_sum.l = list(M24 = M24,
                    M12 = M12, 
                    M11 = M11, 
                    M5 = M5,
                    M3 = M3,
                    M2 = M2,
                    M1 = M1)
 
  
  return(cuts_sum.l)
  
}

lr_dnase = function(cuts_sum.l, pwm, ConsScore, TssScore, label, trainIdx, testIdx, weights){

  scores.l = models.l = list()
  weights.l = list(
    M24 = c(-1,rep(1,4), -1,-1, rep(1,4),-1, 
            -1,rep(1,4), -1,-1, rep(1,4),-1),
    M12 = c(-1,rep(1,4), -1,-1, rep(1,4),-1),
    M11 = c(-1,rep(1,4), -1, rep(1,4),-1),
    M5 = c(1, 1, -1, 1, 1),
    M3 = c(1, -1, 1),
    M2 = c(1, 1),
    M1 = c(1))
  
  
  for (i in 1: (length(cuts_sum.l))){
    
    data.df = data.frame(label, pwm, ConsScore, log2(cuts_sum.l[[i]]+0.1))
    models.l[[i]] <- glm(formula = label ~ ., data = data.df[trainIdx,], family = "binomial")
    predictor.m = as.matrix(cbind(1, data.df[testIdx,-label]))
    if(weights == "LR"){
      predict.logit = predictor.m %*% models.l[[i]]$coefficients
    }else if(weights == "Trained"){
      ref_models.l = get(load("~/DNaseData/Results_forPSB/REST/LRmodels.l"))
      predict.logit = predictor.m %*% ref_models.l[[i]]$coefficients
    }else{
      predict.logit = predictor.m %*% c(-10, 1, 1, weights.l[[i]])
    }
    scores.l[[i]] = 1 / (1 + exp(-predict.logit))
  }
 
  pred.l = list(scores = scores.l, models = models.l)
  
  return(pred.l)
  
  
}



analysis(k = 5, weights = "LR")
# 
# # get model coefficients
# analysis(k = 1, weights = "LR")
# 
# # equal weights
# analysis(k = 1, weights = "Equal")

# weights trained from other TFs
# analysis(k = 1, weights = "Trained")


