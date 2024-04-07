# title: singleCell.R
# author: gumingmu
# date: 2024-04-07

# set parameters-----------------------
rm(list = ls())
gc()

wd <- "D:/work/singlecell/"

wd_c <- paste0(wd, "seurat")
if(!dir.exists(wd_c)){
  dir.create(wd_c)
}
setwd(wd)



# loading package--------------------
library(Seurat)
library(ggplot2)
library(SingleR)
library(dplyr)
library(clustree)
library(celldex)
library(RColorBrewer)
library(SingleCellExperiment)
library(EnsDb.Mmusculus.v79)
library(showtext)
library(reshape)
library(monocle)
library(ggsci)
font.add("kaishu", "simkai.ttf")
showtext_auto()



# loading data----------------------------
dir <- "./matrix/A-1_matrix/"
outdir <- paste0("./seurat/", strsplit(dir, "/")[[1]][3])
dir.create(outdir)
adj.matrix <- Read10X(dir)
srat <- CreateSeuratObject(adj.matrix,project = "pbmc10k") 
dim(srat)
gene_anno.d <- read.table("matrix/A-1_matrix/genes.tsv")


adj.matrix <- NULL
meta <- srat@meta.data
dim(meta)
summary(meta$nCount_RNA)
summary(meta$nFeature_RNA)



dir <- list.dirs("./matrix/", recursive = F)[1:2]
outdir <- paste0("./seurat/", "2_double/")
dir.create(outdir)
adj.matrix <- Read10X(dir)
srat <- CreateSeuratObject(adj.matrix) 
dim(srat)
table(srat$orig.ident)
# 16591 10072 
gene_anno.d <- read.table("matrix/A-1_matrix/genes.tsv")




# qc--------------------------
srat[["percent.mt"]] <- PercentageFeatureSet(srat, pattern = "^mt-")
srat[["percent.rb"]] <- PercentageFeatureSet(srat, pattern = "^Rp[sl]")

VlnPlot(srat, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rb"),ncol = 4,pt.size = 0.1) &
  scale_x_discrete(labels = rep("", 4)) &
  theme(plot.title = element_text(size=10))
ggsave(paste0(outdir, "/vln.pdf"), width = 8, height = 6)
ggplot(srat@meta.data[, 2:5] %>% melt, aes(x = variable, y = value, fill = variable))+
  geom_boxplot(outlier.color = NA)+
  scale_fill_d3()+
  xlab("")+
  ylab("")+
  facet_wrap( ~variable, scale = "free", nrow = 1)+
  theme_bw()+
  theme(legend.position = "none")
ggsave(paste0(outdir, "/box.pdf"), width = 8, height = 6)


FeatureScatter(srat, feature1 = "nCount_RNA", feature2 = "percent.mt")
ggsave(paste0(outdir, "/cor.pdf"), width = 6, height = 6)

cut_gene <- 500
cut_mt <- 8

srat[['QC']] <- "Pass"
# srat[['QC']] <- ifelse(srat@meta.data$Is_doublet == 'True','Doublet','Pass')
srat[['QC']] <- ifelse(srat@meta.data$nFeature_RNA < cut_gene & srat@meta.data$QC == 'Pass','Low_nFeature',srat@meta.data$QC)
srat[['QC']] <- ifelse(srat@meta.data$nFeature_RNA < cut_gene & srat@meta.data$QC != 'Pass' & srat@meta.data$QC != 'Low_nFeature',paste('Low_nFeature',srat@meta.data$QC,sep = ','),srat@meta.data$QC)
srat[['QC']] <- ifelse(srat@meta.data$percent.mt > cut_mt & srat@meta.data$QC == 'Pass','High_MT',srat@meta.data$QC)
srat[['QC']] <- ifelse(srat@meta.data$nFeature_RNA < cut_gene & srat@meta.data$QC != 'Pass' & srat@meta.data$QC != 'High_MT',paste('High_MT',srat@meta.data$QC,sep = ','),srat@meta.data$QC)
table(srat@meta.data[,c(1,6)])

ensdb_genes <- genes(EnsDb.Mmusculus.v79)
pc_name <-  names(ensdb_genes)[ensdb_genes$gene_biotype == "protein_coding"]
index1 <- which(gene_anno.d$V1 %in% pc_name)
re <- c()
for(i in 1:50){
  temp <- which(rowSums(srat@assays$RNA@counts >= 1) >=i) %>% length
  re <- c(re, temp)
}
ggplot(data.frame(x = 1:50, y = re), aes(x, y))+
  geom_point()+
  theme_bw()
ggsave(paste0(outdir, "/geneInCell.pdf"), width = 8, height = 6)

index2 <- which(rowSums(srat@assays$RNA@counts >= 1) >=10)
index <- intersect(index1, index2)

srat <- subset(srat, subset = QC == "Pass")
srat <- srat[index, ]
dim(srat)

table(srat@meta.data$orig.ident)
saveRDS(srat, paste0(outdir, "/data_qc.rdf"))




if(!"srat" %in% ls()){
  srat <- readRDS(paste0(outdir, "/data_qc.rdf"))
}
dim(srat)
table(srat$orig.ident)



# resample----------
set.seed(1998)
index1 <- lapply(split(1:ncol(srat), srat$orig.ident),function(x){
  sample(x, 9000)
}) %>% unlist()
srat <- srat[, index1]
table(srat$orig.ident)



# normalize----------------------
srat <- NormalizeData(srat)

srat <- FindVariableFeatures(srat, selection.method = "vst", nfeatures = 3000)
top10 <- head(VariableFeatures(srat), 10)
top10 
plot1 <- VariableFeaturePlot(srat)
LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
ggsave(paste0(outdir, "/featureGene.pdf"), width = 6, height = 6)

#scale
all.genes <- rownames(srat)
srat <- ScaleData(srat, features = all.genes)

# saveRDS(srat, paste0(outdir, "/data_norm.rdf"))
gc()






# reduction-------------------
if(!"srat" %in% ls()){
  srat <- readRDS(paste0(outdir, "/data_norm.rdf"))
}

srat <- RunPCA(srat, features = VariableFeatures(object = srat))
DimPlot(srat, reduction = "pca")
ElbowPlot(srat)

dim <- 15
srat <- FindNeighbors(srat, dims = 1:dim)

#resolution
res.used <- seq(0.8,2.7,by=0.3)

for(i in res.used){
  srat <- FindClusters(object = srat, verbose = T, resolution = res.used)
}

clus.tree.out <- clustree(srat) +
  theme(legend.position = "bottom") +
  scale_color_brewer(palette = "Set1") +
  scale_edge_color_continuous(low = "grey80", high = "red")
ggsave(paste0(outdir, "/cluster.pdf"), width = 10, height = 10)


resu <- 1.7
srat <- FindClusters(srat, resolution = resu)

srat <- RunUMAP(srat, dims = 1:dim, verbose = F)
plot1 <- DimPlot(srat,label.size = 4,repel = T,label = T)
srat <- RunTSNE(srat, dims = 1:dim)
plot2 <- DimPlot(srat, reduction = "tsne",label.size = 4,repel = T,label = T)
plot <- plot1 + plot2
ggsave(paste0(outdir, "/reduct.pdf"), width = 12, height = 5)

saveRDS(srat, paste0(outdir, "/umap.rds"))



# annotation------------------------------------
if(!"srat" %in% ls()){
  srat <- readRDS(paste0(outdir, "/umap.rds"))
}
if(!"anno_cl" %in% ls() & file.exists(paste0(outdir, "/gapdh.txt"))){
  temp <- read.table(paste0(outdir, "/anno_immu.txt"), header = F, sep = "\t")
}

# plot1 <- FeaturePlot(srat,"Gapdh", reduction = "tsne") + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
# plot2 <- FeaturePlot(srat,"Eno1", reduction = "tsne") + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
# plot3 <- plot1+plot2
# ggsave(paste0(outdir, "/gapdh.pdf"), width = 5, height = 10)
# index1 <- grep("\\bGapdh\\b", srat@assays$RNA@counts@Dimnames[[1]])
# index2 <- grep("\\bEno1\\b", srat@assays$RNA@counts@Dimnames[[1]])
# mean1 <- data.frame(cluster = srat@meta.data$seurat_clusters, exp = srat@assays$RNA@data[index1, ] %>% as.numeric) %>% group_by(cluster) %>% summarise(gap_mean = round(median(exp), 3))
# mean2 <- data.frame(cluster = srat@meta.data$seurat_clusters, exp = srat@assays$RNA@data[index2, ] %>% as.numeric) %>% group_by(cluster) %>% summarise(eno_mean = round(median(exp), 3))
# mean <- mean1 %>% 
#   mutate(eno_mean = mean2$eno_mean)
# write.table(mean, paste0(outdir, "gapdh.txt"), row.names = F, quote = F, sep = "\t")


# marker
Tcell <- c("Cd3e")
DC <- c("Flt3")
B <- c("Cd79b", "Ly6d")
mp <- c("Adgre1", "Itgam", "Csf1r")
Np <- c("Ly6g", "Arg2", "S100a9")

marker <- list("Tcell" = Tcell,  "DC" = DC, "B"= B, "NK" = NK,"mp" = mp, "Np" = Np)

DotPlot(srat, features = unique(marker %>% unlist), group.by = "seurat_clusters")+
  RotatedAxis()+
  scale_x_discrete("")+
  scale_y_discrete("")
ggsave(paste0(outdir, "/marker.pdf"), width = 8, height = 5)



#singleR -------------------------
mouseImmu <- get(load("SingleR_ref/ref_Mouse_imm.RData"))

ref <- mouseImmu
label <- ref$label.fine
singleR.re <- SingleR(test = srat@assays$RNA@data, ref = ref, labels = label, clusters = srat@meta.data$seurat_clusters)

temp <- singleR.re$pruned.labels %>% as.data.frame()
temp$cluster <- rep(0:(length(unique(srat$seurat_clusters))-1))

write.table(temp, paste0(outdir, "/anno_immu.txt"), row.names = F, quote = F, sep = "\t", col.names = F)

DotPlot(srat, features = unique(marker %>% unlist), group.by = "seurat_clusters")+
  RotatedAxis()+
  scale_x_discrete("")+
  scale_y_discrete("", label = temp[,1] %>% unlist)
ggsave(paste0(outdir, "/marker_singleR.pdf"), width = 8*1.3, height = 5*1.3)

srat@meta.data$monaco.fine <- factor(srat@meta.data$seurat_clusters, labels = singleR.re$pruned.labels)

srat <- SetIdent(srat, value = "monaco.fine")
plot1 <- DimPlot(srat, label = T , repel = T, label.size = 3)+NoLegend()
plot2 <- DimPlot(srat, label = T , repel = T, label.size = 3, reduction = "tsne")+NoLegend()
plot <- plot1 + plot2
ggsave(paste0(outdir, "/anno.pdf"), plot, width = 13, height = 6)




# filter unknow cluster------------------
srat <- srat[, !srat$seurat_clusters %in% c(26, 27, 28)]

cluster_m <- c("Neutrophils", "Macrophages", rep("Neutrophils", 2), "Macrophages", "Neutrophils", "Macrophages", rep("Neutrophils", 2), "Macrophages",  "T cells", rep("Macrophages", 6), "T cells", "Macrophages", "Neutrophils", "Macrophages", "DC", "Macrophages", "B cell", "T cells", "Macrophages")
cluster_m <- factor(srat$seurat_clusters, labels = cluster_m)
srat$type <- cluster_m



# T cell subtype-------------------------
outdir_sub = paste0(outdir, "T cells/")
dir.create(outdir_sub)
subsrat <- srat[, srat$type == "T cells"]
dim(subsrat)

subsrat <- FindVariableFeatures(subsrat, selection.method = "vst", nfeatures = 3000)
top10 <- head(VariableFeatures(subsrat), 10)
top10 
plot1 <- VariableFeaturePlot(subsrat)
LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
ggsave(paste0(outdir_sub, "/featureGene.pdf"), width = 6, height = 6)
#scale
all.genes <- rownames(subsrat)
subsrat <- ScaleData(subsrat, features = all.genes)

subsrat <- RunPCA(subsrat, features = VariableFeatures(object = subsrat))
DimPlot(subsrat, reduction = "pca")
ElbowPlot(subsrat)

dim <- 18
subsrat <- FindNeighbors(subsrat, dims = 1:dim)

#resolution
res.used <- seq(0.1, 1,by=0.1)

for(i in res.used){
  subsrat <- FindClusters(object = subsrat, verbose = T, resolution = res.used)
}

clus.tree.out <- clustree(subsrat) +
  theme(legend.position = "bottom") +
  scale_color_brewer(palette = "Set1") +
  scale_edge_color_continuous(low = "grey80", high = "red")
ggsave(paste0(outdir_sub, "/cluster.pdf"), width = 15, height = 15)


resu <- 0.25
resu <- 5

subsrat <- FindClusters(subsrat, resolution = resu)
table(subsrat$seurat_clusters)


subsrat <- RunUMAP(subsrat, dims = 1:dim, verbose = F)
plot1 <- DimPlot(subsrat,label.size = 4,repel = T,label = T)
subsrat <- RunTSNE(subsrat, dims = 1:dim, verbose = F)
plot2 <- DimPlot(subsrat, reduction = "tsne",label.size = 4,repel = T,label = T)
plot <- plot1 + plot2
ggsave(paste0(outdir_sub, "/reduct.pdf"), width = 12, height = 5)

saveRDS(subsrat, paste0(outdir_sub, "/umap.rds"))
# subsrat <- readRDS(paste0(outdir_sub, "/umap.rds"))

subsrat$type <- factor(subsrat$seurat_clusters, labels = paste0("T cells", 1:6))

# plot
color_cd8 <- c("#b912d8", "#488de0", "#cb021b", "#417207", "#f3a71f", "#03bbf3")
marker_T <- c( "Cd3e","Cd8a", "Cd4", "Pdcd1", "Lag3","Gzmb","Ifng", "Ncr1", "Klrk1", "Mki67", "Cdk1", "Cd44", "Il7r", "Cxcr6", "Cd163l1" )

VlnPlot(subsrat, marker_T, group.by = "type", pt.size = 0, cols = color_cd8, split.by = "type", stack  =T, flip  = T )+
  theme(legend.position = "none")
ggsave(paste0(outdir_f, "fig4.pdf"), width = 6*1.5, height = 9*1.5, units = "cm")


color_cd8 <- c("#b912d8", "#488de0", "#cb021b", "#417207", "#f3a71f", "#03bbf3")
input.df <- data.frame(cluster = subsrat$type, group = subsrat$orig.ident)

ggplot(input.df, aes(x = group, fill = cluster))+
  geom_bar(stat = "count", width = 0.7, position = "fill")+
  scale_fill_manual(name = "", values = color_cd8)+
  scale_y_continuous(expand = c(0, 0, 0, 0), labels = scales::percent)+
  scale_x_discrete(labels = group)+
  labs(y = "Percentage of T cells\nin each cluster", x = "")+
  theme_classic()+
  theme(text = element_text(size = 7), legend.key.size = unit(0.3, "cm"))
ggsave(paste0(outdir_f, "fig5.pdf"), width = 6, height = 6, units = "cm")



# T cell DE enrich-------------------------
de.df <- FindMarkers(subsrat, ident.1 = "2", ident.2 = "1", group.by = "orig.ident", logfc.threshold = 0,min.pct = 0)
fwrite(de.df, paste0(outdir_f, "fig6_de.txt"), quote = F, sep = "\t", row.names = T)

de.df <- read.table(paste0(outdir_f, "fig6_de.txt"))
go.bp <- read.gmt("D:/source/GO/gsea/mouse/m2.cp.v2022.1.Mm.symbols.gmt")
fclist <- arrange(de.df, -avg_log2FC)[, 2]
names(fclist) <- arrange(de.df, -avg_log2FC) %>% rownames()

gsea.r <- GSEA(fclist,TERM2GENE =  go.bp, minGSSize = 0, pvalueCutoff = 1, eps = 0)
fwrite(gsea.r@result, paste0(outdir_f, "fig6_gsea_h2.txt"), sep = "\t")
#plot
library(gggsea)
a1 <- read.table("3.sum/gsea/1/gsea_go.bp.txt", header = T, sep = "\t", row.names = 1)[, c(1,5,6,8)]
term <- a1$ID[6]

go.bp <- read.gmt("D:/source/GO/gsea/mouse/m2.cp.v2022.1.Mm.symbols.gmt")
gene.list <- go.bp$gene[go.bp$term == "BIOCARTA_TCYTOTOXIC_PATHWAY"] %>% list
names(gene.list) <- "TCYTOTOXIC PATHWAY"

plot.df <- fread(paste0(outdir_f, "fig6_gsea_h2.txt"))[, c(1, 5:7)] %>% as.data.frame()
colnames(plot.df) <- c("pathway", "NES", "pval", "FDR")
plot.df <- plot.df[plot.df$pathway == "BIOCARTA_TCYTOTOXIC_PATHWAY", ]

plot1.df <- gseaCurve(fclist, gene.list)
plot1.df$color =  colorRampPalette(c(rep(c("#153e6a"), 1),"#32aecc","white","#f45b54",rep(c("#89111c"), 1)))(nrow(plot1.df)) %>% rev()

ggplot()+
  geom_gsea(plot1.df, linecolor = color_mutli[3], prettyGSEA = T, linesize = 0.5, ticksize = 0.1,)+
  geom_text(aes(x =10500, y = 0.75, label = paste0("NES = " ,plot.df$NES %>% round(., 2)
                                                   ,"\nFDR = ", plot.df$FDR %>% round(., 3) ))
            , size = 7*5/14)+
  labs(y = "Enrichment score (ES)", x = "Rank in Ordered Dataset")+
  theme_bw()+
  theme(text = element_text(size = 7))
ggsave(paste0(outdir_f, "fig6.pdf"), width = 6, height = 5, units = "cm")



#DC-----------------------
outdir <- paste0("./seurat/", "2_double/")
outdir_sub = paste0(outdir, "DC/")
dir.create(outdir_sub)
subsrat <- srat[, srat$type == "DC"]
dim(subsrat)

subsrat <- FindVariableFeatures(subsrat, selection.method = "vst", nfeatures = 3000)

top10 <- head(VariableFeatures(subsrat), 10)
top10 
plot1 <- VariableFeaturePlot(subsrat)
LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
ggsave(paste0(outdir_sub, "/featureGene.pdf"), width = 6, height = 6)
#scale
all.genes <- rownames(subsrat)
subsrat <- ScaleData(subsrat, features = all.genes)

subsrat <- RunPCA(subsrat, features = VariableFeatures(object = subsrat))
DimPlot(subsrat, reduction = "pca")
ElbowPlot(subsrat)

dim <- 18
subsrat <- FindNeighbors(subsrat, dims = 1:dim)

#resolution
res.used <- seq(0.1, 0.7,by=0.05)
# 
for(i in res.used){
  subsrat <- FindClusters(object = subsrat, verbose = T, resolution = res.used)
}
# 
clus.tree.out <- clustree(subsrat) +
  theme(legend.position = "bottom") +
  scale_color_brewer(palette = "Set1") +
  scale_edge_color_continuous(low = "grey80", high = "red")
ggsave(paste0(outdir_sub, "/cluster.pdf"), width = 15, height = 15)

resu <- 0.25
resu <- 0.35

subsrat <- FindClusters(subsrat, resolution = resu)
table(subsrat$seurat_clusters)

subsrat <- RunUMAP(subsrat, dims = 1:dim, verbose = F)
plot1 <- DimPlot(subsrat,label.size = 4,repel = T,label = T)
subsrat <- RunTSNE(subsrat, dims = 1:dim, verbose = F)
plot2 <- DimPlot(subsrat, reduction = "tsne",label.size = 4,repel = T,label = T)
plot <- plot1 + plot2
ggsave(paste0(outdir_sub, "/reduct.pdf"), width = 12, height = 5)

saveRDS(subsrat, paste0(outdir_sub, "/umap.rds"))
# subsrat <- readRDS(paste0(outdir_sub, "/umap.rds"))

subsrat$type <- factor(subsrat$seurat_clusters, labels = paste0("DC_", 1:4))

subsrat <- readRDS(paste0("seurat/2_double/DC", "/umap.rds"))
subsrat$type <- factor(subsrat$seurat_clusters, labels = paste0("DC_", 1:4))

color_cd8 <- c("#b912d8", "#488de0", "#cb021b", "#417207")
DimPlot(subsrat, reduction = "tsne",pt.size = 0.4,repel = F,label = T, group.by = "type", label.size = 8*5/14, cols = color_cd8,label.color = "black")+
  scale_x_continuous(expand = c(0.1, 0.1))+
  scale_y_continuous(expand = c(0.1, 0.1))+
  labs(x = "t-SNE_1", y = "t-SNE_2",title = "") +
  theme_classic()+
  theme(legend.position = "none")
ggsave(paste0(outdir_f, "fig7.pdf"), width = 6*1.5, height = 4.5*1.5, units = "cm")

# marker
marker_T <- c( "Xcr1", "Cadm1", "H2-DMb2", "Itgax", "Ccr7", "Fscn1", "Cd80", "Cd86", "Itgae")

VlnPlot(subsrat, marker_T, group.by = "type", pt.size = 0, cols = color_cd8, split.by = "type", stack  =T, flip  = T )+
  theme(legend.position = "none")
ggsave(paste0(outdir_f, "fig8.pdf"), width = 5*1.5, height = 6*1.5, units = "cm")


input.df <- data.frame(cluster = subsrat$type, group = subsrat$orig.ident)

ggplot(input.df, aes(x = group, fill = cluster))+
  geom_bar(stat = "count", width = 0.7, position = "fill")+
  scale_fill_manual(name = "", values = color_cd8)+
  scale_y_continuous(expand = c(0, 0, 0, 0), labels = scales::percent)+
  scale_x_discrete(labels = group)+
  labs(y = "Percentage of DCs\nin each cluster", x = "")+
  theme_classic()+
  theme(text = element_text(size = 7), legend.key.size = unit(0.3, "cm"))
ggsave(paste0(outdir_f, "fig9.pdf"), width = 6, height = 6, units = "cm")


FeaturePlot(subsrat, c("Itgae", "Cd80"), reduction = "tsne", cols = c("lightgrey", "darkred"), pt.size = 0.5, label.size = 0.5)
ggsave(paste0(outdir_f, "fig10.pdf"), width = 12*1.5, height = 6*1.5, units = "cm")

subsrat$group <- factor(subsrat$orig.ident, labels = group[1:2])
DimPlot(subsrat, group.by = "group", reduction = "tsne", cols = color_group, pt.size = 0.5)+
  labs(title = "")+
  theme(legend.position = "top")
ggsave(paste0(outdir_f, "fig10-1.pdf"), width = 6*1.5, height = 6.5*1.5, units = "cm")

# gsea
de.df <- FindMarkers(subsrat, ident.1 = "2", ident.2 = "1", group.by = "orig.ident", logfc.threshold = 0,min.pct = 0)
fwrite(de.df, paste0(outdir_f, "fig11_de.txt"), quote = F, sep = "\t", row.names = T)

de.df <- read.table(paste0(outdir_f, "fig11_de.txt"))
go.bp <- read.gmt("D:/source/GO/gsea/mouse/m5.go.bp.v2022.1.Mm.symbols.gmt")
fclist <- arrange(de.df, -avg_log2FC)[, 2]
names(fclist) <- arrange(de.df, -avg_log2FC) %>% rownames()
# fclist <- fclist[abs(fclist) < 0.1]
gsea.r <- GSEA(fclist,TERM2GENE =  go.bp, pvalueCutoff = 1, eps = 0)
fwrite(gsea.r@result, paste0(outdir_f, "fig11_gsea_h2.txt"), sep = "\t")
#plot
library(gggsea)

go.bp <- read.gmt("D:/source/GO/gsea/mouse/m5.go.bp.v2022.1.Mm.symbols.gmt")
gene.list <- go.bp$gene[go.bp$term == "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION"] %>% list
names(gene.list) <- "ANTIGEN PROCESSING AND PRESENTATION"

plot.df <- fread(paste0(outdir_f, "fig11_gsea_h2.txt"))[, c(1, 5:7)] %>% as.data.frame()
colnames(plot.df) <- c("pathway", "NES", "pval", "FDR")
plot.df <- plot.df[plot.df$pathway == "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION", ]

plot1.df <- gseaCurve(fclist, gene.list)
plot1.df$color =  colorRampPalette(c(rep(c("#153e6a"), 1),"#32aecc","white","#f45b54",rep(c("#89111c"), 1)))(nrow(plot1.df)) %>% rev()

ggplot()+
  geom_gsea(plot1.df, linecolor = color_mutli[3], prettyGSEA = T, linesize = 0.5, ticksize = 0.1,)+
  geom_text(aes(x =11500, y = 0.45, label = paste0("NES = " ,plot.df$NES %>% round(., 2)
                                                   ,"\nFDR = ", plot.df$FDR %>% round(., 3) ))
            , size = 7*5/14)+
  labs(y = "Enrichment score (ES)", x = "Rank in Ordered Dataset")+
  theme_bw()+
  theme(text = element_text(size = 7))
ggsave(paste0(outdir_f, "fig11.pdf"), width = 6, height = 5, units = "cm")



