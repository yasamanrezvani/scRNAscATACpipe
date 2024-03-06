


getGOtab <- function(dir) {
  
  all.files <- list.files(dir)
  
  all.clust.items <- list()
  for(f in all.files){
    nn <- gsub('\\.tsv', '', f)
    cluster <- strsplit(nn, split = '_')[[1]][2]
    GF <- strsplit(nn, split = '_')[[1]][1]
    tmp <- read_tsv(paste(dir, f, sep = ''))
    tmp$GF <- GF
    tmp$cluster <- cluster
    all.clust.items <- c(all.clust.items, list(tmp))
  }
  all.clust.items <- do.call(rbind, all.clust.items)
  colnames(all.clust.items) <- gsub("P-value", "pval", colnames(all.clust.items))
  
  filtered.Go <- all.clust.items
  filtered.Go$cluster <- factor(all.clust.items$cluster, levels = sort(unique(filtered.Go$cluster)))
  filtered.Go$ID <- factor(filtered.Go$ID, level=unique(filtered.Go$ID))
  filtered.Go$Name <- factor(filtered.Go$Name, level=unique(filtered.Go$Name))
  colnames(filtered.Go) <- gsub(" ", "_", colnames(filtered.Go))
  
  
  return(all.clust.items)
}


########################
## transition color ####
########################

in.dir <- '../Output_YR/toxo_cdc/ME49_59/tables/rna_sig_markers_rna_trans_sum_GO_toxodb_OutPut/'
all.clust.items <- getGOtab(dir  = in.dir)
write.xlsx(all.clust.items, '../Output/toxo_cdc/ME49_59/tables/all_GO/rna_sig_markers_rna_trans_GO_term.xlsx')



filtered.Go <- all.clust.items %>% arrange(cluster, Benjamini) %>% distinct() %>%
  group_by(cluster) %>% mutate(rank = row_number()) %>%
  #dplyr::filter(pval < 0.05 )
  dplyr::filter(Benjamini < 0.1 & rank < 15) %>% arrange(cluster, Benjamini)

filtered.Go$cluster <- factor(filtered.Go$cluster, levels = sort(unique(filtered.Go$cluster)))
filtered.Go$ID <- factor(filtered.Go$ID, level=unique(filtered.Go$ID))
filtered.Go$Name <- factor(filtered.Go$Name, level = rev(unique(filtered.Go$Name)))
write.xlsx(filtered.Go, '../Output/toxo_cdc/ME49_59/tables/all_GO/rna_sig_markers_rna_trans_GO_term_filt.xlsx')

p <- ggplot(filtered.Go, aes(x = cluster, y = Name,fill = cluster, color = cluster)) + 
  geom_point(aes( size =  -log(Benjamini))) +
  theme_bw(base_size = 14) +
  scale_color_manual(values = c("T1" = "#ff9a00", 'T2' = '#9ca820', 'T3' = '#615FB1', 'T4' = '#8f139f')) +
  scale_fill_manual(values = c("T1" = "#ff9a00", 'T2' = '#9ca820', 'T3' = '#615FB1', 'T4' = '#8f139f')) +
  #scale_colour_gradient(limits=c(0, 0.01), low="red") +
  ylab(NULL) + xlab(NULL) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face="bold", colour = "black")) + 
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face="bold", colour = "black")) +
  #theme(legend.position="none") +
  theme(strip.background = element_rect(colour="black", fill="white", 
                                        size=0.5, linetype="solid")) + guides(color = FALSE)+
  #ggtitle(strsplit(in.dir,split = "/")[[1]][5]) +
  theme(plot.title = element_text(size = 10)) + 
  xlab("cluster") + ylab("GO term") +
  theme(axis.title = element_text(size = 14, colour = "black", face = "bold"))

plot(p)

ggsave(filename="../Output/toxo_cdc/ME49_59/figures_paper/rna_sig_markers_rna_trans_GO_term_top.pdf",
       plot=p ,
       width = 8, height = 8,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)
##############################
## rna transition markers cluster based ##
##############################
in.dir <- '../Output_YR/toxo_cdc/ME49_59/tables/rna_sig_markers_rna_trans_clust_based_sum/'
all.clust.items <- getGOtab(dir  = in.dir)
all.clust.items$group <- gsub("C[0-9]","", all.clust.items$cluster)
all.clust.items$cluster.RNA <-   gsub("T[0-9]","", all.clust.items$cluster)

all.clust.items <- all.clust.items %>% 
  mutate(new.ord.clust = 
           case_when(group == "T1" & cluster.RNA == "C1" ~ "D4",
                     group == "T1" & cluster.RNA == "C2" ~ "D2", 
                     group == "T1" & cluster.RNA == "C3" ~ "D3", 
                     group == "T1" & cluster.RNA == "C4" ~ "D1",
                     group == "T2" & cluster.RNA == "C1" ~ "D4",
                     group == "T2" & cluster.RNA == "C2" ~ "D2", 
                     group == "T2" & cluster.RNA == "C3" ~ "D1", 
                     group == "T2" & cluster.RNA == "C4" ~ "D3", 
                     group == "T3" & cluster.RNA == "C1" ~ "D4",
                     group == "T3" & cluster.RNA == "C2" ~ "D3",
                     group == "T3" & cluster.RNA == "C3" ~ "D2", 
                     group == "T3" & cluster.RNA == "C4" ~ "D1", 
                     group == "T4" & cluster.RNA == "C1" ~ "D3", 
                     group == "T4" & cluster.RNA == "C2" ~ "D4", 
                     group == "T4" & cluster.RNA == "C3" ~ "D1", 
                     group == "T4" & cluster.RNA == "C4" ~ "D2", 
                     TRUE ~ "NA"))
all.clust.items$new.cluster.RNA <- gsub("D", "C",all.clust.items$new.ord.clust)
all.clust.items$cluster2 <- paste(all.clust.items$group, all.clust.items$new.cluster.RNA, sep = "")


write.xlsx(all.clust.items, '../Output/toxo_cdc/ME49_59/tables/all_GO/rna_sig_markers_rna_trans_clust_based_GO_term_ordered_clust.xlsx')


filtered.Go <- all.clust.items 
filtered.Go <- filtered.Go %>% arrange(cluster2, Benjamini) %>% distinct() %>%
  group_by(cluster2) %>% mutate(rank = row_number()) %>%
  #dplyr::filter(pval < 0.05 )
  dplyr::filter(Benjamini < 0.1 & rank < 5) %>% arrange(cluster2, Benjamini)

filtered.Go$ID <- factor(filtered.Go$ID, level=unique(filtered.Go$ID))
filtered.Go$Name <- factor(filtered.Go$Name, level=rev(unique(filtered.Go$Name)))
filtered.Go <- filtered.Go %>%  arrange(cluster2, Benjamini)
write.xlsx(filtered.Go, '../Output_YR/toxo_cdc/ME49_59/tables/all_GO/rna_sig_markers_rna_trans_clust_based_GO_term_filt_ordered_clust.xlsx')


# transition color clust based
p <- ggplot(filtered.Go, aes(x = cluster2, y = Name,fill = cluster2, color = cluster2)) + 
  geom_point(aes( size =  -log(Benjamini))) +
  theme_bw(base_size = 14) +
  scale_color_manual(values = c("T1C1" = "#ff9a00", 'T1C2' = '#ff9a00', 'T1C3' = '#ff9a00', 'T1C4' = '#ff9a00',
                                'T2C1' = '#9ca820', 'T2C2' = '#9ca820', 'T2C3' = '#9ca820', 'T2C4' = '#9ca820', 
                                'T3C1' = '#615FB1', 'T3C2' = '#615FB1', 'T3C3' = '#615FB1', 'T3C4' = '#615FB1',
                                'T4C1' = '#8f139f', 'T4C2' = '#8f139f', 'T4C3' = '#8f139f', 'T4C4' = '#8f139f')) +
  scale_fill_manual(values = c("T1C1" = "#ff9a00", 'T1C2' = '#ff9a00', 'T1C3' = '#ff9a00', 'T1C4' = '#ff9a00',
                               'T2C1' = '#9ca820', 'T2C2' = '#9ca820', 'T2C3' = '#9ca820', 'T2C4' = '#9ca820', 
                               'T3C1' = '#615FB1', 'T3C2' = '#615FB1', 'T3C3' = '#615FB1', 'T3C4' = '#615FB1',
                               'T4C1' = '#8f139f', 'T4C2' = '#8f139f', 'T4C3' = '#8f139f', 'T4C4' = '#8f139f')) +
  #scale_colour_gradient(limits=c(0, 0.01), low="red") +
  ylab(NULL) + xlab(NULL) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face="bold", colour = "black")) + 
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 14, face="bold", colour = "black")) +
  #theme(legend.position="none") +
  theme(strip.background = element_rect(colour="black", fill="white", 
                                        size=0.5, linetype="solid")) + guides(color = FALSE)+
  #ggtitle(strsplit(in.dir,split = "/")[[1]][5]) +
  theme(plot.title = element_text(size = 10)) + 
  xlab("cluster") + ylab("GO term") +
  theme(axis.title = element_text(size = 14, colour = "black", face = "bold"))

plot(p)
ggsave(filename="../Output_YR//toxo_cdc/ME49_59/figures_paper/rna_sig_markers_rna_trans_clust_based_GO_term_ordered_clust.pdf",
       plot=p ,
       width = 8, height = 10,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)



#################################################
### cut&run and global/phase based DEGs #########
#################################################
in.dir <- '../Input_YR//toxo_cdc/cutNrun/all_macs2_old_new_batch_NO_Filter/BAMM_analysis/HighConfPeaks/GO/' ## phase based
all.clust.items <- getGOtab(dir  = in.dir)
write.xlsx(all.clust.items, '../Output/toxo_cdc/ME49_59/tables/all_GO/HighConfPeaks_phase_based_DEG.xlsx')

filtered.Go <- all.clust.items %>% arrange(cluster, Benjamini) %>% distinct() %>%
  group_by(cluster) %>% mutate(rank = row_number()) %>%
  #dplyr::filter(pval < 0.05 )
  dplyr::filter(Benjamini < 0.1 & rank < 30) %>% arrange(cluster, Benjamini)

filtered.Go$cluster <- factor(filtered.Go$cluster, levels = sort(unique(filtered.Go$cluster)))
filtered.Go$ID <- factor(filtered.Go$ID, level=unique(filtered.Go$ID))
filtered.Go$Name <- factor(filtered.Go$Name, level = rev(unique(filtered.Go$Name)))
write.xlsx(filtered.Go, '../Output_YR/toxo_cdc/ME49_59/tables/all_GO/HighConfPeaks_phase_based_DEG_filt.xlsx')

p <- ggplot(filtered.Go, aes(x = cluster, y = Name)) + 
  geom_point(aes(colour = cluster, size =  -log(Benjamini))) +
  theme_bw(base_size = 14) +
  scale_color_manual(values = c("downReg" = "#6565bf","upReg" ='#ee5d6c')) +
  scale_fill_manual(values = c("downReg" = "#6565bf","upReg" ='#ee5d6c'))+
  #scale_colour_gradient(limits=c(0, 0.01), low="red") +
  ylab(NULL) + xlab(NULL) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, face="bold")) + 
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  #theme(legend.position="none") +
  theme(strip.background = element_rect(colour="black", fill="white", 
                                        size=0.5, linetype="solid")) + guides(color = FALSE)+
  #ggtitle(strsplit(in.dir,split = "/")[[1]][5]) +
  ggtitle("High Conf CutRun peaks Global DEGs in KD_vs_WT")+
  theme(plot.title = element_text(size = 9))

plot(p)

ggsave(filename="../Output_YR/toxo_cdc/ME49_59/figures_paper/HighConfPeaks_phase_based_DEG_GO_term.pdf",
       plot=p ,
       width = 7, height = 10,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)



## global
in.dir <- '../Input/toxo_cdc/cutNrun/all_macs2_old_new_batch_NO_Filter/BAMM_analysis/HighConfPeaks_globalDEG/GO/'
all.clust.items <- getGOtab(dir  = in.dir)
filtered.Go <- all.clust.items %>% arrange(cluster, Benjamini) %>% distinct() %>%
  group_by(cluster) %>% mutate(rank = row_number()) %>%
  #dplyr::filter(pval < 0.05 )
  dplyr::filter(Benjamini < 0.1 & rank < 30) %>% arrange(cluster, Benjamini)

write.xlsx(all.clust.items, '../Output/toxo_cdc/ME49_59/tables/all_GO/HighConfPeaks_globalDEG.xlsx')


filtered.Go$cluster <- factor(filtered.Go$cluster, levels = sort(unique(filtered.Go$cluster)))
filtered.Go$ID <- factor(filtered.Go$ID, level=unique(filtered.Go$ID))
filtered.Go$Name <- factor(filtered.Go$Name, level = rev(unique(filtered.Go$Name)))
write.xlsx(filtered.Go, '../Output/toxo_cdc/ME49_59/tables/all_GO/HighConfPeaks_globalDEG_filt.xlsx')

p <- ggplot(filtered.Go, aes(x = cluster, y = Name)) + 
  geom_point(aes(colour = cluster, size =  -log(Benjamini))) +
  theme_bw(base_size = 14) +
  scale_color_manual(values = c("downReg" = "#6565bf","upReg" ='#ee5d6c')) +
  scale_fill_manual(values = c("downReg" = "#6565bf","upReg" ='#ee5d6c'))+
  #scale_colour_gradient(limits=c(0, 0.01), low="red") +
  ylab(NULL) + xlab(NULL) + 
  theme(axis.text.x = element_text(angle = 0,colour = "black",hjust = 0.5, size = 14, face="bold", )) + 
  theme(axis.text.y = element_text(angle = 0,colour = "black" ,size = 12, hjust = 1, face="bold")) +
  #theme(legend.position="none") +
  theme(strip.background = element_rect(colour="black", fill="white", 
                                        size=0.5, linetype="solid")) + guides(color = FALSE)+
  #ggtitle(strsplit(in.dir,split = "/")[[1]][5]) +
  #ggtitle("High Conf CutRun peaks Global DEGs in KD_vs_WT")+
  theme(plot.title = element_text(size = 9))

plot(p)

ggsave(filename="../Output/toxo_cdc/ME49_59/tables/all_GO/HighConfPeaks_globalDEG_GO_term.pdf",
       plot=p ,
       width = 9, height = 9,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)


###########################################################
############## atac sub clusters ##########################
###########################################################
in.dir <- '../Output/toxo_cdc/ME49_59/tables/atac_clusters_within_rna_tran_dtw_clusters_manual_GO_result/'
all.clust.items <- getGOtab(dir  = in.dir)

all.clust.items$cluster[all.clust.items$cluster == 'T3C2C2'] <- 'T3C2C1'
all.clust.items$cluster[all.clust.items$cluster == 'T4C2C2'] <- 'T4C2C1'
all.clust.items$cluster[all.clust.items$cluster == 'T4C3C3'] <- 'T4C3C1'



all.clust.items.ord <- all.clust.items %>% 
  mutate(atac.clust.ord = 
           case_when(cluster == "T1C1C1" ~ "D1",
                     cluster == "T1C2C1" ~ "D2", 
                     cluster == "T1C2C2" ~ "D1", 
                     cluster == "T1C2C3" ~ "D3", 
                     cluster == "T1C3C1" ~ "D1", 
                     cluster == "T1C3C2" ~ "D2",
                     cluster == "T1C3C3" ~ "D3",
                     cluster == "T1C4C1" ~ "D2",
                     cluster == "T1C4C2" ~ "D1",
                     cluster == "T2C1C1" ~ "D1",
                     cluster == "T2C2C1" ~ "D3",
                     cluster == "T2C2C2" ~ "D2",
                     cluster == "T2C2C3" ~ "D1",
                     cluster == "T2C3C1" ~ "D2",
                     cluster == "T2C3C2" ~ "D3",
                     cluster == "T2C3C3" ~ "D1",
                     cluster == "T2C4C1" ~ "D2",
                     cluster == "T2C4C2" ~ "D1",
                     cluster == "T3C1C1" ~ "D1",
                     cluster == "T3C2C1" ~ "D1",
                     cluster == "T3C3C1" ~ "D1",
                     cluster == "T3C3C2" ~ "D3",
                     cluster == "T3C3C3" ~ "D2",
                     cluster == "T3C4C1" ~ "D2",
                     cluster == "T3C4C2" ~ "D1",
                     cluster == "T4C1C1" ~ "D2",
                     cluster == "T4C1C2" ~ "D1",
                     cluster == "T4C2C1" ~ "D1",
                     cluster == "T4C3C1" ~ "D1",
                     cluster == "T4C4C1" ~ "D2",
                     cluster == "T4C4C2" ~ "D1",
                     TRUE ~ "NA"))
all.clust.items.ord$atac.clust.ord <- gsub("\\D", "C", all.clust.items.ord$atac.clust.ord) 
all.clust.items.ord$cluster.ord <- paste(substr(all.clust.items.ord$cluster, 1, 4),all.clust.items.ord$atac.clust.ord , sep = "")
write.xlsx(all.clust.items, '../Output/toxo_cdc/ME49_59/tables/all_GO/atac_sub_clusters_GO_term_ordered.xlsx')


filtered.Go <- all.clust.items.ord %>% arrange(cluster.ord, Benjamini) %>% distinct() %>%
  group_by(cluster.ord) %>% mutate(rank = row_number()) %>%
  #dplyr::filter(pval < 0.05 )
  dplyr::filter(Benjamini < 0.1 & rank < 5) %>% arrange(cluster.ord, Benjamini)

filtered.Go$cluster.ord <- factor(filtered.Go$cluster.ord, levels = sort(unique(filtered.Go$cluster.ord)))
filtered.Go$ID <- factor(filtered.Go$ID, level=unique(filtered.Go$ID))
filtered.Go$Name <- factor(filtered.Go$Name, level = rev(unique(filtered.Go$Name)))
write.xlsx(filtered.Go, '../Output/toxo_cdc/ME49_59/tables/all_GO/atac_sub_clusters_filt_GO_term_ordered.xlsx')

 #scale_colour_gradient(limits=c(0, 0.01)
mycolors <- c(rep("#ff9a00",8),rep('#9ca820',6), rep('#615FB1',5), rep('#8f139f',4))
names(mycolors) <- factor(unique(filtered.Go$cluster.ord), levels = unique(filtered.Go$cluster.ord))

p <- ggplot(filtered.Go, aes(x = cluster.ord, y = Name, color = cluster.ord)) + 
  geom_point(aes( size =  -log(Benjamini))) +
  theme_bw(base_size = 14) +
  scale_colour_manual(name = filtered.Go$cluster.ord, values = mycolors)+
  #scale_colour_gradient(limits=c(0, 0.01), low="red") +
  ylab(NULL) + xlab(NULL) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16, face="bold")) + 
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(legend.position="none") +
  theme(strip.background = element_rect(colour="black", fill="white", 
                                        size=0.5, linetype="solid")) + guides(color = FALSE)+
  #ggtitle(strsplit(in.dir,split = "/")[[1]][5]) +
  theme(plot.title = element_text(size = 10))

plot(p)

ggsave(filename ="../Output/toxo_cdc/ME49_59/tables/all_GO/atac_sub_clusters_GO_term_ordered.pdf",
       plot=p ,
       width = 14, height = 12,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)

