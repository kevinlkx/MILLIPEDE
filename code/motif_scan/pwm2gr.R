### convert pwm scores to GenomicRanges object

library(GenomicRanges)

## parameters
args <- commandArgs(trailingOnly = TRUE)
motifID <- args[1]
tfName <- args[2]
pwm_score <- args[3]

## load pwm scores and convert to a GenomicRanges object
tf_pwm.df <- read.table(paste0("~/MNaseData/Model/PWM/motif_scan/PWMoutput/motif", motifID, "_score", pwm_score, ".txt"), header = TRUE)

tf.gr <- GRanges(
  seqnames = Rle(as.integer(as.roman(as.character(tf_pwm.df$chr)))),
  ranges = IRanges(start = tf_pwm.df$start, end = tf_pwm.df$end),
  strand = Rle(strand(tf_pwm.df$strand)), score = tf_pwm.df$score)

dir_tf <- "~/MNaseData/tf_coord_updated/"
dir.create(dir_tf, showWarnings = FALSE, recursive = TRUE)

save(tf.gr, file = paste0(dir_tf, tfName, "_score", pwm_score,".gr"))

cat("finish", tfName)


