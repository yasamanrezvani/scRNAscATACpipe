
##Cut and Run quality check
in.dir <- "../Input_sub/toxo_cdc/cutNrun/all_macs2_old_new_batch_NO_Filter/"
files <- list.files(in.dir, pattern = ".narrowPeak")
f.names <- gsub("\\.narrowPeak", "", files)

# read all called peaks (new AP2_TY (new) vs 4 controls + MiseqA_2mm (old) vs 4 controls)
qval <- 0.05
all.peaks <- list()
for (f in files) {
  tmp <- read.table(paste(in.dir, f, sep = ''), header=F, sep="\t", quote = NULL)
  tmp <- tmp %>% filter(!grepl("KE.*", V1)) 
  tmp <- tmp %>% filter(V9 > -log10(qval))
  all.peaks <- c(all.peaks, list(tmp))
}
names(all.peaks) <- f.names
all.peaks.tab <- do.call("rbind", all.peaks[1:4])
all.peaks.tab$CasevsControl <- gsub("_peak.*", "", all.peaks.tab$V4) 
all.peaks.tab$Case <- unlist(lapply(strsplit(all.peaks.tab$CasevsControl, "_vs_"), "[[", 1))
all.peaks.tab$Control <- unlist(lapply(strsplit(all.peaks.tab$CasevsControl, "_vs_"), "[[", 2))
all.peaks.tab <- all.peaks.tab %>% `rownames<-`(NULL)
names(all.peaks.tab)[1:10] <- c('chr', "strt", "end", "peak_name", "V5", "V6", "SignalValue", "-log10(pVal)", "-log10(qVal)", "peak")

## Mieq with 6 different concentration 
in.dir.Miseq <- "../Input_copy_2/toxo_cdc/cutNrun/paired_end_230207/MACS2/"
files.Miseq <- list.files(in.dir.Miseq, pattern = ".narrowPeak")
f.names.Miseq <- gsub("\\.narrowPeak", "", files.Miseq)

# read all called peaks (new AP2_TY (new) vs 4 controls + MiseqA_2mm (old) vs 4 controls)
qval <- 0.05
all.peaks.Miseq <- list()
for (f in files.Miseq) {
  tmp <- read.table(paste(in.dir.Miseq, f, sep = ''), header=F, sep="\t", quote = NULL)
  tmp <- tmp %>% filter(!grepl("KE.*", V1)) 
  tmp <- tmp %>% filter(V9 > -log10(qval))
  all.peaks.Miseq <- c(all.peaks.Miseq, list(tmp))
}
names(all.peaks.Miseq) <- f.names.Miseq

all.peaks.tab.Miseq <- do.call("rbind", all.peaks.Miseq)
all.peaks.tab.Miseq$CasevsControl <- gsub("_peak.*", "", all.peaks.tab.Miseq$V4) 
all.peaks.tab.Miseq$Case <- unlist(lapply(strsplit(all.peaks.tab.Miseq$CasevsControl, "_vs_"), "[[", 1))
all.peaks.tab.Miseq$Control <- unlist(lapply(strsplit(all.peaks.tab.Miseq$CasevsControl, "_vs_"), "[[", 2))
all.peaks.tab.Miseq <- all.peaks.tab.Miseq %>% `rownames<-`(NULL)
names(all.peaks.tab.Miseq)[1:10] <- c('chr', "strt", "end", "peak_name", "V5", "V6", "SignalValue", "-log10(pVal)", "-log10(qVal)", "peak")

peaks.tab <- rbind(all.peaks.tab, all.peaks.tab.Miseq)
#peaks.tab <- peaks.tab %>% filter(Control == "RH_Negative_S2")

ggplot(peaks.tab, aes(x =  Case, y = SignalValue, color = Case)) +
  geom_boxplot() +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
