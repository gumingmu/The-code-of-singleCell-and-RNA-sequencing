# title: bulk.R
# author: gumingmu
# date: 2024-04-07



# set parameters--------------------------------
rm(list = ls())
gc()

library(tidyverse)
source("D:/work/SLE+circRNA/bioanalysis/Rdata/.R/function/function.R")
load("D:/source/color/color_single.Rdata")

wd <- "D:/work/singlecell/3.rnaseq"
matrix_file <- "D:/work/singlecell/3.rnaseq/rawData/matrix/count.txt"
outdir <- "D:/work/singlecell/3.rnaseq/"

setwd(wd)



# de -------------------------------------
library(edgeR)

exp.d <- read.table(matrix_file, sep = "\t", header = T, check.names = F)
exp.df <- exp.d %>% 
  column_to_rownames("Geneid") %>% 
  dplyr::select(starts_with("./")) %>% 
  rename_all(.funs = funs(paste0(rep(c("PTT+EXOSOME INHIBITOR", "PTT", "PBS"), each = 1), "_", rep(1:3, times = 4))))

anno.df <- data.frame(sample = colnames(exp.df), group = my.strsplit(colnames(exp.df), split = "_", ori = 1), order = my.strsplit(colnames(exp.df), split = "_", ori = 2))

outdir1 <- paste0(outdir, "rawData/annotation/")
dir.create(outdir1)
write.table(anno.df, paste0(outdir1, "anno.txt"), sep = "\t", quote = F)

exp.df <- exp.df[, order(anno.df$group, anno.df$order)]
anno.df <- anno.df %>% 
  arrange(match(sample, colnames(exp.df)))

exp.dge <- DGEList(counts = exp.df, group = anno.df$group)
filter1 <- filterByExpr(exp.dge)
filter2 <- grepl("protein_coding", rownames(exp.df))
filter <- filter1&filter2
exp.dge <- exp.dge[filter, , keep.lib.size = F]
exp.dge <- calcNormFactors(exp.dge)

outdir1 <- paste0(outdir, "pipeline_data/")
dir.create(outdir1)
saveRDS(exp.dge, paste0(outdir1, "edgeR.rds"))

cpm.df <- cpm(exp.dge)

outdir1 <- paste0(outdir, "rawData/matrix/")
write.table(cpm.df, paste0(outdir1, "cpm.txt"), quote = F, sep = "\t")

length <- exp.d %>% 
  filter(Geneid %in% rownames(exp.dge$counts)) %>% 
  dplyr::select(Length) %>% 
  unlist()
tpm.df <- exp.dge$counts %>% 
  sweep(., 1, length, FUN = "/")
colSum <- colSums(tpm.df)/10^6
tpm.df <- sweep(tpm.df, 2, colSum, FUN = "/")

outdir1 <- paste0(outdir, "rawData/matrix/")
write.table(tpm.df, paste0(outdir1, "tpm.txt"), quote = F, sep = "\t")

exp.dge <- readRDS("pipeline_data/edgeR.rds")
groupPair <- combn(unique(exp.dge$samples$group) %>% as.character, 2) %>% as.data.frame() %>% as.list()

outdir1 <- paste0(outdir, "1.de/edgeR_pad/sample3/FDR0.05/")
dir.create(outdir1)
speices_num <- 2
goidchange <- F
for(i in 1:length(groupPair)){
  outdir2 <- paste0(outdir1, groupPair[[i]][2], "-", groupPair[[i]][1],"/")
  dir.create(outdir2)
  
  temp.dge <- exp.dge[, exp.dge$samples$group %in% groupPair[[i]]]
  design <- model.matrix(~temp.dge$samples$group)
  
  temp.dge <- estimateDisp(temp.dge, design, robust = T)
  temp.fit <- glmFit(temp.dge)
  temp.lrt <- glmLRT(temp.fit)
  temp.de <- topTags(temp.lrt,n=nrow(temp.lrt$table))$table
  temp.de <- temp.de %>% 
    mutate(sign = ifelse(FDR < 0.05, "Y", "N"), change = ifelse(logFC > 0, ifelse(logFC > 1, "2x Up", "Up"), ifelse(logFC < -1, "2x Down", "Down")))
  temp.summ <- table(temp.de[, 6:7])
  
  write.table(temp.de, paste0(outdir2, groupPair[[i]][2], "-", groupPair[[i]][1],".txt"), sep = "\t", quote = F)
  write.table(temp.summ, paste0(outdir2, groupPair[[i]][2], "-", groupPair[[i]][1],"_summary.txt"), sep = "\t", quote = F)
  
  ## plot
  plot.df <- temp.de %>% 
    mutate(change1 = ifelse(logFC > 0, "Up", "Down"))
  ggplot(plot.df[plot.df$sign == "Y", ], aes(x = change1))+
    geom_bar(aes(fill = change1), stat = "count", width = 0.6, position = "stack")+
    scale_fill_manual(values = color_pair_fill[c(3,1)])+
    geom_text(aes(label = ..count..), stat = "count", nudge_y = 0.05*max(table(plot.df$change1[plot.df$sign == "Y"])), size = 5/14*6)+
    labs(x = "", y = "Counts")+
    theme_minimal()+
    theme( text = element_text(size = 6),panel.grid = element_line(size = 0.1), plot.margin = unit(c(0.1,0.1,0.1,0.1),      "cm"),legend.position = "")
  ggsave(paste0(outdir2, groupPair[[i]][2], "-", groupPair[[i]][1],"_bar.pdf"), width = 5, height = 5, units = "cm")
  
  # enrich ----------------------------------------------------
  for(j in 1:2){
    if(j == 1){
      outdir3 <- paste0(outdir2, "Up/")
      input <- rownames(temp.de)[temp.de$FDR<0.05&temp.de$logFC>1] 
    }else{
      outdir3 <- paste0(outdir2, "Down/")
      input <- rownames(temp.de)[temp.de$FDR<0.05&temp.de$logFC<-1] 
    }
    
    dir.create(outdir3, recursive = T)
    input <- my.strsplit(input)
    backgene <- rownames(temp.de) %>% my.strsplit()
    
    
    if(speices_num == 1){
      suppressPackageStartupMessages(library(org.Hs.eg.db))
      lib_org <- "org.Hs.eg.db"
    }else{
      suppressPackageStartupMessages(library(org.Mm.eg.db))
      lib_org <- "org.Mm.eg.db"
    }
    
    suppressPackageStartupMessages(library(ggplot2))
    
    DEG.gene_symbol <- input
    DEG.entrez_id = mapIds(get(lib_org), DEG.gene_symbol, 'ENTREZID', 'SYMBOL') %>% suppressMessages()
    DEG.entrez_id = na.omit(DEG.entrez_id)
    
    black.gene_symbol <- backgene
    black.entrez_id = mapIds(get(lib_org), black.gene_symbol, 'ENTREZID', 'SYMBOL') %>% suppressMessages()
    black.entrez_id = na.omit(black.entrez_id)
    
    
    go_method <- c("BP")
    for(Method in go_method){
      outplotfile <- paste0(outdir3, "GO_", Method, ".pdf")
      outgofile <- paste0(outdir3, "GO_", Method, ".txt")
      
      # erich.go.BP = enrichGO(gene = DEG.entrez_id,
      #                        OrgDb = get(lib_org),
      #                        keyType = "ENTREZID",
      #                        ont = Method ,
      #                        pvalueCutoff = 0.5,
      #                        qvalueCutoff = 0.5
      #                        ,universe = black.entrez_id
      # )
      
      erich.go.BP = enrichKEGG(gene = DEG.entrez_id,
                               keyType = "kegg",
                               organism = "mmu",
                               pvalueCutoff = 1,
                               qvalueCutoff = 1,
                               universe = black.entrez_id
      )
      #进行entreid转换
      test <- erich.go.BP@result
      
      if(goidchange){
        test[,10] <- ""
        colnames(test)[10] <- "symbol"
        for(i in 1:nrow(test)){
          entre <- test[i, 8] %>% strsplit(., split = "/") %>% unlist
          symbol <- mapIds(get(lib_org), entre,'SYMBOL', 'ENTREZID') %>% na.omit() %>% suppressMessages()
          test[i, 10] <- paste0(symbol, collapse = "/")
        }
      }
      
      
      write.table(test, file = outgofile, quote = F, col.names = F, sep = "\t")

      test=as.data.frame(erich.go.BP@result)
      test=test[1:30,]

      library(stringr)
      gr1 <- as.numeric(str_split(test$GeneRatio,"/",simplify = T)[,1])
      gr2 <- as.numeric(str_split(test$GeneRatio,"/",simplify = T)[,2])
      bg1 <- as.numeric(str_split(test$BgRatio,"/",simplify = T)[,1])
      bg2 <- as.numeric(str_split(test$BgRatio,"/",simplify = T)[,2])
      test$fold <- (gr1/gr2)/(bg1/bg2)
      test <- arrange(test,-log10(p.adjust))
      test$Description = factor(test$Description,levels = test$Description,ordered = T)
      p <- ggplot(test,aes(x = -log10(p.adjust),y = Description))+
        geom_point(aes(size = fold,color = Count))+
        scale_color_gradient(high = "#fb9182", low = "#40d4dc") +
        xlab("-log10(p.adjust)")+
        labs(size = "Fold Enrichment", color = "Counts")+
        guides(color = guide_colorbar(reverse = F))+
        ylab("")+
        geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "#ff5000", size = 0.8)+
        theme(panel.background = element_blank(), panel.grid = element_blank(),panel.grid.major.y = element_line(linetype = 2, colour = "grey80"),
              legend.background = element_rect(fill = "white"),text = element_text(size = 12),
              legend.key = element_blank(), axis.text = element_text( color = "black", size = 12),
              axis.line.x = element_line(),
              axis.line.y = element_line(),
        )
      print(p)
      
      ggsave(outplotfile, p, width = 16, height = 9)
    }
    
    
  }
  
}



# enrich plot
# PTT-PBS---------------------------
test <- fread("D:/work/singlecell/3.rnaseq/1.de/edgeR_pad/sample3/PTT-PBS/Up/GO_BP.txt")[, -1]
colnames(test)[c(1:6, 9)] <- c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "Count")

gr1 <- as.numeric(str_split(test$GeneRatio,"/",simplify = T)[,1])
gr2 <- as.numeric(str_split(test$GeneRatio,"/",simplify = T)[,2])
bg1 <- as.numeric(str_split(test$BgRatio,"/",simplify = T)[,1])
bg2 <- as.numeric(str_split(test$BgRatio,"/",simplify = T)[,2])
test$fold <- (gr1/gr2)/(bg1/bg2)
test <- arrange(test,p.adjust)
test$Description = factor(test$Description,levels = test$Description %>% rev(),ordered = T)

ggplot(test[3:12, ], aes(x = -log10(p.adjust), y = Description))+
  geom_point(aes(color = fold, size = Count, fill = fold), shape = 21)+
  scale_color_gradient(low = "#3900f7", high = "#ee0040")+
  scale_fill_gradient(low = "#3900f7", high = "#ee0040")+
  scale_x_continuous(expand = c(0.2, 0, 0.2, 0))+
  scale_size_continuous(range = c(0.5,1)*4)+
  labs(x = "-log10(FDR)", y = "", color = "GeneRatio", fill = "GeneRatio")+
  theme_bw()+
  theme(text = element_text(size = 8), axis.text.y = element_text(color = "black"), legend.key.size = unit(0.25, "cm"), 
        legend.title = element_text(size = 7))
ggsave(paste0(outdir_f, "fig_enrich_1.pdf"), width = 12, height = 5, units = "cm")

#PTT+EXOSOME INHIBITOR-PBS--------------
test <- fread("D:/work/singlecell/3.rnaseq/1.de/edgeR_pad/sample3/PTT+EXOSOME INHIBITOR-PBS/Up/GO_BP.txt")[, -1]
colnames(test)[c(1:6, 9)] <- c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "Count")

gr1 <- as.numeric(str_split(test$GeneRatio,"/",simplify = T)[,1])
gr2 <- as.numeric(str_split(test$GeneRatio,"/",simplify = T)[,2])
bg1 <- as.numeric(str_split(test$BgRatio,"/",simplify = T)[,1])
bg2 <- as.numeric(str_split(test$BgRatio,"/",simplify = T)[,2])
test$fold <- (gr1/gr2)/(bg1/bg2)
test <- arrange(test,p.adjust)
test$Description = factor(test$Description,levels = test$Description %>% rev(),ordered = T)

ggplot(test[c(9, 19, 24, 30, 36, 38, 49, 50, 56, 75), ], aes(x = -log10(p.adjust), y = Description))+
  geom_point(aes(color = fold, size = Count, fill = fold), shape = 21)+
  scale_color_gradient(low = "#3900f7", high = "#ee0040")+
  scale_fill_gradient(low = "#3900f7", high = "#ee0040")+
  scale_x_continuous(expand = c(0.2, 0, 0.2, 0))+
  scale_size_continuous(range = c(0.5,1)*4)+
  labs(x = "-log10(FDR)", y = "", color = "GeneRatio", fill = "GeneRatio")+
  theme_bw()+
  theme(text = element_text(size = 8), axis.text.y = element_text(color = "black"), legend.key.size = unit(0.2, "cm"), 
        legend.title = element_text(size = 7))
ggsave(paste0(outdir_f, "fig_enrich_2.pdf"), width = 12, height = 5, units = "cm")





# immune cell analysis-----------------
library(estimate)
outdir1 <- paste0(outdir, "immune/")
dir.create(outdir1)

cpm.df <- read.table("rawData/matrix/cpm.txt", header = T, sep = "\t", check.names = F)
cpm.df <- cpm.df %>% 
  rownames_to_column("id") %>% 
  mutate(id = my.strsplit(id)) %>% 
  group_by(id) %>% 
  summarise_all(median) %>% 
  column_to_rownames("id")
write.table(cpm.df, paste0(outdir1, "estimate_input.txt"), sep = "\t", quote = F)

### immuneCellAi-mouse
result.im <- read.table("immune/ImmuCellAI_mouse_abundance_result1.txt", sep = "\t", header = T,row.names = 1)
rownames(result.im) <- colnames(cpm.df)
result.im <- t(result.im) %>% as.data.frame()
library(pheatmap)
color_pheatmap_zscore2 <- colorRampPalette(c(rep(c("#153e6a"), 1),"#32aecc","gray95","#f45b54",rep(c("#89111c"), 1)))(50)
p <- pheatmap(result.im, scale = "row", cluster_cols = F, color = color_pheatmap_zscore2, angle_col = 45, gaps_col = c(3))
ggsave(paste0(outdir1, "cellAi_heat_sample.pdf"),p,  width = 6.5*1.2, height = 8*1.2)

