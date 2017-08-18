library(GenomicRanges)
setwd("~/Dropbox/projects_dropbox/MILLIPEDE/workflow/code/millipede_analysis")

source("lr_functions.R")
source("fitCentipede_footprint.R")

tf <- "reb1"

tfName <- paste(tf,"_macisaac_score4",sep = "")
flank_data <- 100

chipexoList <- c("reb1", "rap1", "gal4", "phd1")

dir_dnase_data <- paste0("../../data/yeast_data/dnase_data/", tfName, "/")
dir_motifs <- paste0('../../data/yeast_data/macisaac_score4/')
dir_binding_sites <- paste0('../../data/yeast_data/binding_sites/')

ThreshAvgMappability <- 0.8

reb1_sites_test.df <- read.table("../../data/test/reb1.cutoff.4.binding.sites.sorted.txt")
colnames(reb1_sites_test.df) <- c("chr", "start", "stop", "strand", "score")
reb1_sites_test.df$chr <- as.factor(as.numeric(as.roman(gsub("chr", "", reb1_sites_test.df$chr))))
reb1_sites_test.gr <- makeGRangesFromDataFrame(reb1_sites_test.df, keep.extra.columns = TRUE)

# reb1_chipexo_test.df <- read.table("../../data/test/reb1.saccer2.chipexo.tsv", sep = '\t', header = TRUE, comment.char = "!")
tf_chipexo.gr <- get(load(paste0(dir_binding_sites, "/", tf, "_chip_exo.gr")))

tf_macIsaac.gr <- get(load(paste0(dir_binding_sites, "/", tf, ".gr")))

reb1_sites_test.gr <- label_macisaac(reb1_sites_test.gr, tf_macIsaac.gr, "within")
reb1_sites_test.gr <- label_chipexo(reb1_sites_test.gr, tf_chipexo.gr, "within")
elementMetadata(reb1_sites_test.gr)[,"label_joint"] <- as.integer(elementMetadata(reb1_sites_test.gr)[,"label"] | elementMetadata(reb1_sites_test.gr)[,"chip_exo"])

cat("matched macIsaac sites number:", sum(elementMetadata(reb1_sites_test.gr)[,"label"] == 1),
    "matched chip-exo sites number:", sum(elementMetadata(reb1_sites_test.gr)[,"chip_exo"] == 1),
    "Joint positive sites number:", sum(elementMetadata(reb1_sites_test.gr)[,"label_joint"] == 1),"\n")


## load reb1 sites (filtered by mappablity) used in PSB paper

dnase_data.l <- get(load(paste(dir_dnase_data, "dnase_cuts_", flank_data, "bp.l", sep = "")))

tf_overlap.gr <- dnase_data.l$tf.gr

tf_sites.gr <- tf_overlap.gr
start(tf_sites.gr) <- start(tf_overlap.gr) + flank_data
end(tf_sites.gr) <- end(tf_overlap.gr) - flank_data

findOverlaps(reb1_sites_test.gr, tf_sites.gr, type="equal")

## load all reb1 motif matches

tf_all_motifs.gr <- get(load(paste0(dir_motifs, "/", tfName, ".gr")))

findOverlaps(reb1_sites_test.gr, tf_all_motifs.gr, type="equal")

findOverlaps(tf_sites.gr, tf_all_motifs.gr, type="equal")
