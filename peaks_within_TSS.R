library(tidyverse)
library(bedtoolsr)
library(openxlsx)

cutRun.peaks <- readRDS("../Input_sub/toxo_cdc/rds_ME49_59/CutRunPeaks.rds")
nrow(cutRun.peaks)

atac.peaks <- readRDS("../Input_sub/toxo_cdc/rds_ME49_59/atac_peaks.rds")
nrow(atac.peaks)


# ## distance of motifs relative to TSS
gtf.file <- "../Input_sub/Genomes/ToxoDB-59_TgondiiME49.gtf"
gtf <- read.table(gtf.file, header = F, sep = '\t', quote = NULL)
gtf.filt <- gtf %>% dplyr::filter(!grepl('KE.*', V1))
gtf.filt$gene_name <- gsub("\\;.*", "", gsub("transcript_id ", "", gsub("-t.*", "", gtf.filt$V9)))
gtf.filt$gene_name <- gsub("\"", "", gtf.filt$gene_name)
gtf.filt <- gtf.filt %>% filter(V3 == "transcript")

gtf.pos <- gtf.filt %>% filter(V7 == "+")
gtf.pos <- gtf.pos %>% mutate(tss_ext1 = V4 - 2000, tss_ext2 = V4 + 2000)

gtf.neg <- gtf.filt %>% filter(V7 == "-")
gtf.neg <- gtf.neg %>%  mutate(tss_ext1 = V5 - 2000 , tss_ext2 = V5 + 2000) 

gtf.all <- rbind(gtf.pos, gtf.neg) 
gtf.TSS <- gtf.all %>% select(V1, tss_ext1, tss_ext2, gene_name, V7, V4, V5)
gtf.TSS$tss_ext1[which(gtf.TSS$tss_ext1 < 0) ]  <- 0

options(bedtools.path = "/Users/kourosh.zarringhalam/miniconda3/bin/")

peak.cutRun.TSS <- bedtoolsr::bt.intersect(a = cutRun.peaks, b = gtf.TSS, wb = T)
length(unique(peak.cutRun.TSS$V4)) / length(unique(cutRun.peaks$V4))


peak.atac.TSS <- bedtoolsr::bt.intersect(a = atac.peaks, b = gtf.TSS, wb = T)
length(unique(peak.atac.TSS$V4)) / length(unique(atac.peaks$V4))


