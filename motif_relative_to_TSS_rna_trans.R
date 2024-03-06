


in.dir <- "../Input_sub/toxo_cdc/BAMM/rna_transitions_v4/BAMM/all/"
occ.files <- list.files(in.dir)[grep('occurrence', list.files(in.dir))]
in.file <- paste0(in.dir, occ.files)
all.motifs.tab <- lapply(1:length(in.file), function(i){
  mofit <- gsub('cut_run_union_ribo_targs_', '', gsub('\\..*', '', occ.files[i]))
  
  occ.tab <- read.table(in.file[i], header = T, sep = '\t')
  locs <- strsplit(gsub('>', '', occ.tab$seq), split = ':')
  chr <- unlist(lapply(locs, `[[`, 1))
  pos_strd <- unlist(lapply(locs, `[[`, 2))
  strd <- gsub('\\)', '', gsub('.*\\(', '', pos_strd))
  tmp <- strsplit(gsub('\\(.*', '', pos_strd), split = '-')
  strt <- as.numeric(lapply(tmp, `[[`, 1))
  stp <- as.numeric(lapply(tmp, `[[`, 2))
  motif.strd <- occ.tab$strand
  tmp <- strsplit(occ.tab$start..end, split = '\\.')
  motif.strt <- as.numeric(lapply(tmp, `[[`, 1))
  motif.stp <- as.numeric(lapply(tmp, `[[`, 3))
  df <- data.frame(id = gsub('>', '', occ.tab$seq), seq.len = occ.tab$length,
                   chr = chr, start_peak = strt, end_peak = stp, strd = strd, 
                   motif.strt = motif.strt, motif.stp = motif.stp, motif.strd = motif.strd, 
                   pattern = occ.tab$pattern, motif = mofit)
  
  df <- df %>% mutate(motif.strt.rel = ifelse(motif.strd == '-', 2 * seq.len - motif.stp, motif.strt),
                      motif.stp.rel = ifelse(motif.strd == '-', 2 * seq.len - motif.strt, motif.stp))
  
  ## new 
  df <- df %>% mutate(motif.strt.abs = start_peak + motif.strt.rel,
                      motif.stp.abs = start_peak + motif.stp.rel)
  
  return(df)
})

all.motifs.tab <- bind_rows(all.motifs.tab)


peak.genes <- read.xlsx("../Input_sub/toxo_scATAC_MJ_ME49_59/peak_gene_assigned_final.xlsx")
colnames(peak.genes) <- gsub("gene_name", "TGME49", colnames(peak.genes))
peak.genes <- peak.genes %>%
  mutate(id = paste0(V1.x, ":", start_peak, "-", end_peak, "(", V11, ")"))

all.motifs.tab.genes <- full_join(all.motifs.tab, peak.genes, by = 'id')
all.motifs.tab.genes <- all.motifs.tab.genes %>% arrange(TGME49, motif.strt.rel, motif.stp.rel)
all.motifs.tab.genes <- all.motifs.tab.genes[!is.na(all.motifs.tab.genes$motif),]


length(unique(all.motifs.tab.genes$TGME49[!is.na(all.motifs.tab.genes$motif)])) 

## NEW
## TSS dist 
all.motifs.tab.genes.motif <- all.motifs.tab.genes %>% 
  dplyr::select(id, TGME49, motif, chr, pattern, motif.strt.abs, motif.stp.abs)

# ## distance of motifs relative to TSS
gtf.file <- "../Input_sub/Genomes/ToxoDB-59_TgondiiME49.gtf"
gtf <- read.table(gtf.file, header = F, sep = '\t', quote = NULL)
gtf.filt <- gtf %>% dplyr::filter(!grepl('KE.*', V1))
gtf.filt$gene_name <- gsub("\\;.*", "", gsub("transcript_id ", "", gsub("-t.*", "", gtf.filt$V9)))
gtf.filt$gene_name <- gsub("\"", "", gtf.filt$gene_name)
gtf.filt <- gtf.filt %>% filter(V3 == "transcript")

all.motifs.tab.genes.motif <- left_join(all.motifs.tab.genes.motif, gtf.filt, by = c("TGME49" = "gene_name") )

all.motifs.tab.genes.motif <-  all.motifs.tab.genes.motif %>% mutate(tss = ifelse((V7 == "+"), V4, V5))
all.motifs.tab.genes.motif <- all.motifs.tab.genes.motif %>% mutate(tss.dist = tss - motif.strt.abs)

all.motifs.tab.genes.motif <- all.motifs.tab.genes.motif %>% filter(motif != "T1_motif_3")

saveRDS(all.motifs.tab.genes.motif, "../Input_sub/toxo_cdc/rds_ME49_59/motif_dist_from_TSS_plt.rds")
p <- ggplot(all.motifs.tab.genes.motif, aes(x = tss.dist)) + 
  geom_histogram(aes(y = ..density..),
                 colour = 1, fill = "white") +
  geom_density(alpha=.2, fill="#FF6666",lwd = 0.7,linetype = 1,colour = 2) +
  theme_bw() + 
  facet_wrap(motif ~ ., ncol = 3) +
  theme(strip.background=element_rect(fill='white', color = 'black'),
        panel.spacing = unit(0.5, "lines"), 
        strip.text.x=element_text(angle=0, hjust=0.5,vjust=0.5, size = 14,face = 'bold'),
        strip.text = element_text(face = "bold", size = 24,  angle = 0), 
        strip.placement = "outside",
        plot.title = element_text(size=16, face = "bold.italic", color = 'black'),
        axis.title = element_blank(),
        axis.title.x = element_text(size=20, face="bold"),
        axis.title.y = element_text(size=20, face="bold"),
        axis.text.y =  element_text(size = 16, face = "bold", colour = "black"),
        axis.text.x =  element_text(size = 16, face = "bold", colour = "black", angle = 35))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(0.5, "lines"),
        ) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,size = 16, face = 'bold', color = 'black'))+
  theme(axis.text.y = element_text(size = 16, face = 'bold', color = 'black'))+
  theme(axis.title = element_text(size = 20, face = "bold", color = "black")) +
  theme(legend.text = element_text(face = "bold", size = 10),
        legend.title = element_text(face = "bold", size = 14))

p

ggsave("../OutPut/toxo_cdc/ME49_59/figures_paper/ribo_motifs/rna_transition_motif_TSS_dist_faceted.pdf", 
       plot = p,height = 8, width = 8, dpi = 300)

all.motifs.tab.genes.motif.list <- split(all.motifs.tab.genes.motif, f = all.motifs.tab.genes.motif$motif)
tss.plt <- lapply(1:length(all.motifs.tab.genes.motif.list), function(i){
  
  p <- ggplot(all.motifs.tab.genes.motif.list[[i]], aes(x = tss.dist)) + 
    geom_histogram(aes(y = ..density..),
                   colour = 1, fill = "white") +
    geom_density(lwd = 0.7,linetype = 1,colour = 2) +
    theme_bw() + 
    facet_wrap(~ motif, scales='fixed') +
    theme(strip.background=element_rect(fill='white', color = 'black'),
          panel.spacing = unit(1.5, "lines"), 
          strip.text.x=element_text(angle=0, hjust=0.5,vjust=0.5, size = 14,face = 'bold'),
          plot.title = element_text(size=16, face = "bold.italic", color = 'black'),
          axis.title = element_blank(),
          #axis.title.x = element_text(size=20, face="bold"),
          #axis.title.y = element_text(size=20, face="bold"),
          axis.text.y =  element_text(size = 10, face = "bold", colour = "black"),
          axis.text.x =  element_text(size = 10, face = "bold", colour = "black", angle = 35))+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          #panel.border = element_blank(),
          panel.spacing = unit(0.5, "lines"),
          strip.text = element_text(face = "bold", size = 24,  angle = 0), 
          strip.placement = "outside") +
    theme(legend.text = element_text(face = "bold", size = 10),
          legend.title = element_text(face = "bold", size = 14))+
    theme(panel.spacing = unit(1.5, "lines"))
  p  
  
  return(p)
})

pp <- grid.arrange(grobs = tss.plt, ncol = 3)
ggsave( "../Output/toxo_cdc/ME49_59/figures_paper/motifs_dist_tss_V2.pdf", 
        plot = pp,
        height = 8,width = 8, dpi = 300)
