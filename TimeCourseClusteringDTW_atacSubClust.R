
library(Seurat)
library(ggplot2)
library(ggfortify)
library(jcolors)
require(gridExtra)
library(grid)
library(gam)
library(princurve)
library(parallel)
library(tidyverse)
library(MyEllipsefit)
library(sctransform)
library(openxlsx)
library(doParallel)
library(tidytext)
library(ggrepel)
library(dtwclust)
library(geomtextpath)
library(bigmemory)


source('./util_funcs.R')


## IDs
prod.desc  <- read.xlsx('../Input_sub/toxo_genomics/genes/ProductDescription_GT1.xlsx')
TGGT1_ME49 <- read.xlsx('../Input_sub/toxo_genomics/Orthologs/TGGT1_ME49 Orthologs.xlsx')


rna_sub <- readRDS('../Input_sub/toxo_cdc/rds_ME49_59/S.O_intra_lables_pt.rds')
atac_sub <- readRDS('../Input_sub/toxo_cdc/rds_ME49_59/S.O_intra_atac_lables_pt.rds')


## Splines
sc.rna.spline.fits <- readRDS('../Input_sub/toxo_cdc/rds_ME49_59/sc_rna_spline_fits_all_genes.rds')
sc.atac.spline.fits <- readRDS('../Input_sub/toxo_cdc/rds_ME49_59/sc_atac_spline_fits_all_genes.rds')

## Turn the data into wide format (time by gene) and center & scale each gene
sc.rna.dtw.wide <- sc.rna.spline.fits %>% 
  pivot_wider(names_from = 'GeneID', values_from = 'y') %>% 
  mutate_at(vars(matches("TGME")), ~scale(., center = T, scale = T)) %>%
  as.data.frame()

na.ind <- which(apply(sc.rna.dtw.wide, 2, function(x) any(is.na(x))))
if(length(na.ind)){
  sc.rna.dtw.wide <- sc.rna.dtw.wide[,-na.ind]
}

sc.atac.dtw.wide <- sc.atac.spline.fits %>% 
  pivot_wider(names_from = 'GeneID', values_from = 'y') %>% 
  mutate_at(vars(matches("TGME")), ~scale(., center = T, scale = T)) %>%
  as.data.frame()

na.ind <- which(apply(sc.atac.dtw.wide, 2, function(x) any(is.na(x))))
if(length(na.ind)){
  sc.atac.dtw.wide <- sc.atac.dtw.wide[,-na.ind]
}

sc.rna.mu.scale <- sc.rna.dtw.wide %>% 
  pivot_longer(-x, names_to = 'GeneID', values_to = 'expr')

sc.atac.mu.scale <- sc.atac.dtw.wide %>% 
  pivot_longer(-x, names_to = 'GeneID', values_to = 'expr')



trans.list <- readRDS("../Input_sub/toxo_cdc/rds_ME49_59/rna_markers_rna_transitions_dtw_clust_list_ordered.rds")
trans.all <- do.call("rbind", trans.list)
trans.all <- trans.all %>% dplyr::select(GeneID, cluster.RNA.ordered, group,Name) %>% distinct()
trans.all$trans.cluster.rna <- paste(trans.all$group, trans.all$cluster.RNA.ordered, sep = "_")
colnames(trans.all) <- gsub("GeneID", "gene_name", colnames(trans.all))
colnames(trans.all) <- gsub("Name", "ProductDescription", colnames(trans.all))
trans.all.sum <- trans.all %>% 
  group_by(trans.cluster.rna) %>% summarise(genes = list(gsub("-", "_", gene_name)), total = n())


trans.clust.rna.list  <- split(trans.all, f = trans.all$trans.cluster.rna)
atac.list <- c()

out.dir <- "../Output/toxo_cdc/ME49_59/tables/atac_clusters_within_rna_tran_dtw_clusters_manual/"

k <- c(2,3,3,2,2,3,3,2,2,2,3,2,2,2,2,2)
atac.list <- lapply(1:length(trans.clust.rna.list), function(i) {
  my.df <- trans.clust.rna.list[[i]]
  df <- clust.atac.df(my.df, num.clust = k[i])
  df$cluster.ATAC <- gsub(" ","", df$cluster.ATAC)
  df$group <- names(trans.clust.rna.list)[i]
  df$trans.rna.atac.clust <- paste(df$group, df$cluster.ATAC, sep = "_")
  df <- left_join(df, prod.desc, by = "GeneID")
  p1 <- plot_atac_trand(df) +
    ggtitle(names(trans.clust.rna.list[i]))
  p1
  
  df.Go.tab <- df %>% dplyr::select(GeneID, trans.cluster.rna, trans.rna.atac.clust) %>% distinct()
  df.Go.tab
})

## T1, C1
## KZ manual
atac_sub_clust <- bind_rows(atac.list)
atac_sub_clust$trans.rna.atac.clust[atac_sub_clust$trans.rna.atac.clust == 'T1_C1_C2'] <- 'T1_C1_C1'
atac_sub_clust$trans.rna.atac.clust[atac_sub_clust$trans.rna.atac.clust == 'T2_C1_C2'] <- 'T2_C1_C1'
atac_sub_clust$trans.rna.atac.clust[atac_sub_clust$trans.rna.atac.clust == 'T3_C1_C2'] <- 'T3_C1_C1'
atac_sub_clust$trans.rna.atac.clust[atac_sub_clust$trans.rna.atac.clust == 'T3_C2_C2'] <- 'T3_C2_C1'
atac_sub_clust$trans.rna.atac.clust[atac_sub_clust$trans.rna.atac.clust == 'T4_C2_C2'] <- 'T4_C2_C1'
atac_sub_clust$trans.rna.atac.clust[atac_sub_clust$trans.rna.atac.clust == 'T4_C3_C2'] <- 'T4_C3_C1'

out.dir <- '../Output/toxo_cdc/ME49_59/tables/atac_clusters_within_rna_tran_dtw_clusters_manual/'
write.xlsx(atac_sub_clust, paste(out.dir, 'rna_atac_sub_clusters', ".xlsx", sep = ""))



## manual plotting of atac subclusters profiles
i <- 7
my.df <- trans.clust.rna.list[[i]]

expr.atac.tab <- get_rna_atac_profile(rna.splines = sc.rna.spline.fits, 
                                      atac.splines = sc.atac.spline.fits, 
                                      genes.tab = my.df, scale = T)

## if you wan the atac profiles as well, then do not filter data == "scRNA
expr.tab <- expr.atac.tab %>% filter(data == "scATAC")
p1 <- plot_rna_atac(expr.tab)
p1 
ggsave("../OutPut/toxo_cdc/ME49_59/figures_paper/T4C3C1.pdf", 
       plot = p1, height = 3, width = 5, dpi = 300)

################################################################

## manual reordering the sub atac clusters 
out.dir <- '../Output/toxo_cdc/ME49_59/tables/atac_clusters_within_rna_tran_dtw_clusters_manual/'
atac_sub_clust <- read.xlsx(paste(out.dir, 'rna_atac_sub_clusters', ".xlsx", sep = ""))
atac_sub_clust$atac.clust <- unlist(lapply(strsplit(atac_sub_clust$trans.rna.atac.clust, "_"), "[[", 3))

## order the rna clusters manually according to their peaks

atac_sub_clust <- atac_sub_clust %>% 
  mutate(atac.clust.ord = 
           case_when(trans.rna.atac.clust == "T1_C1_C1" ~ "D1",
                     trans.rna.atac.clust == "T1_C2_C1" ~ "D2", 
                     trans.rna.atac.clust == "T1_C2_C2" ~ "D1", 
                     trans.rna.atac.clust == "T1_C2_C3" ~ "D3", 
                     trans.rna.atac.clust == "T1_C3_C1" ~ "D1", 
                     trans.rna.atac.clust == "T1_C3_C2" ~ "D2",
                     trans.rna.atac.clust == "T1_C3_C3" ~ "D3",
                     trans.rna.atac.clust == "T1_C4_C1" ~ "D2",
                     trans.rna.atac.clust == "T1_C4_C2" ~ "D1",
                     trans.rna.atac.clust == "T2_C1_C1" ~ "D1",
                     trans.rna.atac.clust == "T2_C2_C1" ~ "D3",
                     trans.rna.atac.clust == "T2_C2_C2" ~ "D2",
                     trans.rna.atac.clust == "T2_C2_C3" ~ "D1",
                     trans.rna.atac.clust == "T2_C3_C1" ~ "D2",
                     trans.rna.atac.clust == "T2_C3_C2" ~ "D3",
                     trans.rna.atac.clust == "T2_C3_C3" ~ "D1",
                     trans.rna.atac.clust == "T2_C4_C1" ~ "D2",
                     trans.rna.atac.clust == "T2_C4_C2" ~ "D1",
                     trans.rna.atac.clust == "T3_C1_C1" ~ "D1",
                     trans.rna.atac.clust == "T3_C2_C1" ~ "D1",
                     trans.rna.atac.clust == "T3_C3_C1" ~ "D1",
                     trans.rna.atac.clust == "T3_C3_C2" ~ "D3",
                     trans.rna.atac.clust == "T3_C3_C3" ~ "D2",
                     trans.rna.atac.clust == "T3_C4_C1" ~ "D2",
                     trans.rna.atac.clust == "T3_C4_C2" ~ "D1",
                     trans.rna.atac.clust == "T4_C1_C1" ~ "D2",
                     trans.rna.atac.clust == "T4_C1_C2" ~ "D1",
                     trans.rna.atac.clust == "T4_C2_C1" ~ "D1",
                     trans.rna.atac.clust == "T4_C3_C1" ~ "D1",
                     trans.rna.atac.clust == "T4_C4_C1" ~ "D2",
                     trans.rna.atac.clust == "T4_C4_C2" ~ "D1",
                     TRUE ~ "NA"))
atac_sub_clust$atac.clust.ord <- gsub("\\D", "C", atac_sub_clust$atac.clust.ord) 
atac_sub_clust$trans.rna.atac.clust.ord <- paste(atac_sub_clust$trans.cluster.rna, atac_sub_clust$atac.clust.ord, sep = "_")

out.dir <- '../Output/toxo_cdc/ME49_59/tables/atac_clusters_within_rna_tran_dtw_clusters_manual/'
write.xlsx(atac_sub_clust, paste(out.dir, "rna_atac_sub_clusters_ordered.xlsx", sep = ""))

### plot ordered atac sub clusters 



plot_atac_sub_clust_ordered <- function(tab){
  p  <- ggplot(tab, aes(x= time,y=normExpr)) +
    geom_path(aes(color = gene_name),alpha = 0.8, size = 0.8)+ 
    theme_bw() +
    theme(legend.position = "right") +
    ylab('normExpr') + xlab('Time') +
    # theme(axis.text.x = element_text(angle = 0, size = 20, face="bold", colour = "black")) +
    # theme(axis.text.y = element_text(angle = 0, size = 20, face="bold", colour = "black")) +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(), 
          axis.title = element_blank()) + 
    theme(strip.background = element_rect(colour="black", fill="white",
                                          size=0.5, linetype="solid")) +
    theme(strip.text = element_text(size = 20, face="bold", angle = 0)) + 
    
    coord_cartesian(xlim = c(0,6.5)) + 
    facet_grid(atac.clust.ord ~ ., scales = 'free', space = 'free') +
    theme(
      plot.title = element_text(size=14, face = "bold.italic", color = 'red'),
      # axis.title.x = element_text(size=20, face="bold", hjust = 1),
      # axis.title.y = element_text(size=20, face="bold")
    ) + 
    theme(#legend.position = c(0.15, 0.85),
      legend.position = 'none',
      legend.title = element_text(colour="black", size=12, 
                                  face="bold"),
      legend.text = element_text(colour="black", size=12, 
                                 face="bold"))
  return(p)
  
} 

out.dir <- '../Output/toxo_cdc/ME49_59/tables/atac_clusters_within_rna_tran_dtw_clusters_manual/'
atac_sub_clust <- read.xlsx(paste(out.dir, "rna_atac_sub_clusters_ordered.xlsx", sep = ""))
colnames(atac_sub_clust)[1] <- "TGME49"

## NEW
prod.desc <- read.xlsx("../Input_sub/Toxo_genomics/genes/MJ_annotation_07_27_2023.xlsx")
prod.desc$TGME49 <- gsub("_", "-", prod.desc$TGME49)

atac_sub_clust <- inner_join(atac_sub_clust, prod.desc,  , by = "TGME49")
atac_sub_clust.supp <- atac_sub_clust %>% dplyr::select(TGME49, trans.cluster.rna, trans.rna.atac.clust.ord, ProductDescription,new.name )

write.xlsx(atac_sub_clust.supp, "../OutPut/toxo_cdc/ME49_59/tables/rna_atac_sub_clusters_prod_desc_ordered_final.xlsx")

atac_sub_clust.list <- split(atac_sub_clust, f= atac_sub_clust$trans.cluster.rna)
names(atac_sub_clust.list)

atac.sub.clust.plt <- lapply(1:length(atac_sub_clust.list), function(i){
  
  my.df <- atac_sub_clust.list[[i]] 
  colnames(my.df) <- gsub("TGME49", "gene_name", colnames(my.df))
  clust.name <- unique(my.df$trans.rna.atac.clust.ord)
  
  expr.atac.tab <- get_rna_atac_profile(rna.splines = sc.rna.spline.fits, 
                                        atac.splines = sc.atac.spline.fits, 
                                        genes.tab = my.df, scale = T)
  expr.atac.tab <- inner_join(my.df, expr.atac.tab, by = c("gene_name" = "GeneID"))
  expr.tab <- expr.atac.tab %>% filter(data == "scATAC")
  
  p <- plot_atac_sub_clust_ordered(expr.tab)
  
  return(p)
})
names(atac.sub.clust.plt) <- names(atac_sub_clust.list)

out.dir <- "../OutPut/toxo_cdc/ME49_59/figures_paper/atac_sub_cluster_profiles_ordered/"

lapply(1:length(atac.sub.clust.plt), function(i){
  
  p <- atac.sub.clust.plt[[i]]
  clust.name <- names(atac_sub_clust.list[i])
  ggsave(paste(out.dir, paste(clust.name, ".pdf", sep = ""), sep = ""), 
         plot = p, height = 3, width = 5, dpi = 300)
  
})

i <- 15
p <- atac.sub.clust.plt[[i]]
ggsave(paste(out.dir, paste("T4_C3_C1", ".pdf", sep = ""), sep = ""), 
       plot = p, height = 1.5, width = 5, dpi = 300)


write.xlsx(atac_sub_clust.supp, "../OutPut/toxo_cdc/ME49_59/tables/rna_atac_sub_clusters_prod_desc_ordered_final.xlsx")

atac_sub_clust.list <- split(atac_sub_clust, f= atac_sub_clust$trans.cluster.rna)
names(atac_sub_clust.list)

atac.sub.clust.tab <- lapply(1:length(atac_sub_clust.list), function(i){
  
  my.df <- atac_sub_clust.list[[i]] 
  colnames(my.df) <- gsub("TGME49", "gene_name", colnames(my.df))
  clust.name <- unique(my.df$trans.rna.atac.clust.ord)
  
  expr.atac.tab <- get_rna_atac_profile(rna.splines = sc.rna.spline.fits, 
                                        atac.splines = sc.atac.spline.fits, 
                                        genes.tab = my.df, scale = T)
  expr.atac.tab <- inner_join(my.df, expr.atac.tab, by = c("gene_name" = "GeneID"))
  expr.tab <- expr.atac.tab %>% filter(data == "scATAC")
  
  #p <- plot_atac_sub_clust_ordered(expr.tab)
  
  return(expr.tab)
})
names(atac.sub.clust.tab) <- names(atac_sub_clust.list)

saveRDS(atac.sub.clust.tab, "../Input_sub/toxo_cdc/rds_ME49_59/atac_sub_clusters_access_profiles_list_ordered.rds")



atac.sub.clust.tab <- readRDS( "../Input_sub/toxo_cdc/rds_ME49_59/atac_sub_clusters_access_profiles_list_ordered.rds")

atac.sub.clust.plt <- c()
atac.sub.clust.plt <- lapply(1:length(atac.sub.clust.tab), function(i){
  
  p <- plot_atac_sub_clust_ordered(atac.sub.clust.tab[[i]])
  p
})
names(atac.sub.clust.plt) 
pp <- grid.arrange(grobs = atac.sub.clust.plt, ncol = 4)
ggsave("../OutPut/toxo_cdc/ME49_59/figures_paper/test.pdf", height = 18, width = 16, plot = pp)
