#Load Required Libraries
library(Nebulosa)
library(Seurat)
library(SingleCellExperiment)
library(scater)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(hdf5r)
library(dsb)
library(data.table)
library(cowplot)
library(ggplot2)
library(gridExtra)
library(SeuratWrappers)
library(Azimuth)
library(ggrepel)
library(patchwork)
library(scCustomize)
library(reticulate)
library(circlize)
library(ComplexHeatmap)
library(readxl)
##set path to load data


setwd("C:/Users/axi313/Documents/ECCITE_Seq_Katherine/WNN")
load.path <- "C:/Users/axi313/Documents/ECCITE_Seq_Katherine/saved_R_data/"

load(paste0(load.path,'Seuratv5_isotype_Assay3.RData'))

################################################# Integrate RNA #############################################
DefaultAssay(seurat_isotype) <- 'RNA'
DefaultAssay(seurat_isotype)


rna.list <- SplitObject(seurat_isotype, split.by = "orig.ident")

for (i in 1:length(rna.list)) {
  rna.list[[i]] <- NormalizeData(rna.list[[i]],assay = 'RNA')
  rna.list[[i]] <- FindVariableFeatures(rna.list[[i]], selection.method = "vst", 
                                        nfeatures = 2000,assay = 'RNA')
}

rna.anchors <- FindIntegrationAnchors(object.list = rna.list, dims = 1:30)

rna.integrated <- IntegrateData(anchorset = rna.anchors, dims = 1:30)

rna.integrated <- RenameAssays(rna.integrated, integrated = 'rna.integrated')

# Run the standard workflow for visualization and clustering
rna.integrated <- ScaleData(rna.integrated,assay = 'rna.integrated')
rna.integrated <- RunPCA(rna.integrated, npcs = 30, reduction.name = "pca.rna.integrated")
rna.integrated <- RunUMAP(rna.integrated, reduction = "pca.rna.integrated", dims = 1:30, reduction.name = "umap.rna.integrated")

save(rna.integrated, file=paste0(load.path,"Seuratv5_WNN_rna_integrated.RData"))
#load(paste0(load.path,'Seuratv5_WNN_rna_integrated.RData'))

################################################# Integrate ADT ###############################################
DefaultAssay(rna.integrated) <- 'ADT'
DefaultAssay(rna.integrated)
adt.list <- SplitObject(rna.integrated, split.by = "orig.ident")

# define proteins to use in clustering (non-isptype controls)
prots <- rownames(seurat_isotype@assays$ADT@data)
isotype_genes <- c('Mouse-IgG1', 'Mouse-IgG2a', 'Mouse-IgG2b', 'Rat-IgG2b', 'Rat-IgG1', 'Rat-IgG2a', 'Hamster-IgG')
prots <- setdiff(prots, isotype_genes)
for (i in 1:length(adt.list)) {
  VariableFeatures(adt.list[[i]]) = prots
}

features <- SelectIntegrationFeatures(object.list = adt.list)

adt.list <- lapply(X = adt.list, FUN = function(x) {
  x <- ScaleData(x, features = features)
  x <- RunPCA(x, features = features)
})

anchors <- FindIntegrationAnchors(object.list = adt.list, reduction = "rpca", 
                                  dims = 1:30)

rna_adt_integrated <- IntegrateData(anchorset = anchors, dims = 1:30,new.assay.name = "integrated")
rna_adt_integrated <- RenameAssays(rna_adt_integrated, integrated = 'adt.integrated')

rna_adt_integrated <- ScaleData(rna_adt_integrated,assay = 'adt.integrated')
rna_adt_integrated <- RunPCA(rna_adt_integrated, npcs = 30, reduction.name = "pca.adt.integrated")
rna_adt_integrated <- RunUMAP(rna_adt_integrated, reduction = "pca.adt.integrated", dims = 1:30, reduction.name = "umap.adt.integrated")
save(rna_adt_integrated, file=paste0(load.path,"Seuratv5_WNN_rna_adt_integrated.RData"))

####################################### Seurat WNN default with PCA on dsb normalized protein ###############

seurat_isotype <- rna_adt_integrated


########### Check PCA variance
# SET rna features 
rna.features <- rownames(seurat_isotype@assays$rna.integrated)
seurat_isotype <- ScaleData(seurat_isotype,assay = 'rna.integrated')
seurat_isotype <- RunPCA(seurat_isotype, npcs = 30, reduction.name = "pca.rna.integrated",features=rna.features, assay = 'rna.integrated')
seurat_isotype <- RunUMAP(seurat_isotype, reduction = "pca.rna.integrated", dims = 1:30, reduction.name = "umap.rna.integrated")

### RNA

ElbowPlot(object = seurat_isotype, 
          ndims = 50, reduction = 'pca.rna.integrated')


# Determine percent of variation associated with each PC
pct <- seurat_isotype[["pca.rna.integrated"]]@stdev / sum(seurat_isotype[["pca.rna.integrated"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
co2

# Minimum of the two calculation
pcs <- min(co1, co2)

pcs
# Create a dataframe with values
plot_df <- data.frame(pct = pct, 
                      cumu = cumu, 
                      rank = 1:length(pct))

# Elbow plot to visualize 

ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()

### PCs=16

### ADT
ElbowPlot(object = seurat_isotype, 
          ndims = 50, reduction = 'pca.adt.integrated')

# Determine percent of variation associated with each PC
pct <- seurat_isotype[["pca.adt.integrated"]]@stdev / sum(seurat_isotype[["pca.adt.integrated"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
co2

# Minimum of the two calculation
pcs <- min(co1, co2)

pcs
# Create a dataframe with values
plot_df <- data.frame(pct = pct, 
                      cumu = cumu, 
                      rank = 1:length(pct))

# Elbow plot to visualize 

ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()


### apca=11

# run WNN 
seurat_isotype = FindMultiModalNeighbors(
  seurat_isotype, reduction.list = list("pca.rna.integrated", "pca.adt.integrated"), 
  dims.list = list(1:16, 1:11)
)

# cluster 
seurat_isotype <-RunUMAP(seurat_isotype, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
seurat_isotype <- FindClusters(seurat_isotype, graph.name = "wsnn", 
                               algorithm = 3, 
                               resolution = 0.8,
                               random.seed = 1990)

save(seurat_isotype, file=paste0(load.path,"Seuratv5_WNN_Complete.RData"))


DimPlot(seurat_isotype, reduction='wnn.umap',label=T,label.size = 5,repel=T, group.by = 'predicted.celltype.l2', split.by = 'stim')
DefaultAssay(seurat_isotype) <- 'ADT'
FeaturePlot(seurat_isotype, reduction = 'wnn.umap', features = "CD40LG")
VlnPlot(seurat_isotype, features = "TNFRSF4",  group.by = 'predicted.celltype.l2', split.by = 'stim')

########## Annotation #############

### Cluster Plots

seurat_isotype@meta.data$wsnn_res.0.8
setwd('C:/Users/axi313/Documents/ECCITE_Seq_Katherine/WNN/Cluster_Plots')

DimPlot(seurat_isotype, reduction='wnn.umap',label=T,label.size = 5,repel=T, group.by = 'predicted.celltype.l1')
ggsave('WNN_ClusterPlot_Azimuth1.png',dpi=500, width = 13)

DimPlot(seurat_isotype, reduction='wnn.umap',label=T,label.size = 5,repel=T, group.by = 'predicted.celltype.l2')
ggsave('WNN_ClusterPlot_Azimuth2.png',dpi=500, width = 13)

DimPlot(seurat_isotype, reduction='wnn.umap',label=T,label.size = 5,repel=T, group.by = 'predicted.celltype.l3')
ggsave('WNN_ClusterPlot_Azimuth3.png',dpi=500, width = 13)

DimPlot(seurat_isotype, reduction='wnn.umap',label=T,label.size = 5,repel=T)
ggsave('WNN_ClusterPlot_res0.8.png',dpi=500, width = 13)

DimPlot(seurat_isotype, reduction='wnn.umap',label=T,label.size = 5,repel=T, split.by = 'stim')
ggsave('WNN_ClusterPlot_res0.8_bystim.png',dpi=500, width = 13)

DimPlot(seurat_isotype, reduction='wnn.umap',label=T,label.size = 5,repel=T, split.by = 'PID')
ggsave('WNN_ClusterPlot_res0.8_byPID.png',dpi=500, width = 13)

DimPlot(seurat_isotype, reduction='wnn.umap',label=T,label.size = 5,repel=T, split.by = 'condition')
ggsave('WNN_ClusterPlot_res0.8_bycondition.png',dpi=500, width = 13)




#################################### Heatmaps ##############################################
# create multimodal heatmap 


# find marker genes for the joint clusters 
Idents(seurat_isotype) = "wsnn_res.0.8"
DefaultAssay(seurat_isotype)  = "RNA"
rnade = FindAllMarkers(seurat_isotype, features = rna.features, only.pos = TRUE)
gene_plot = rnade %>% 
  dplyr::filter(avg_log2FC > 1 ) %>%  
  dplyr::group_by(cluster) %>% 
  dplyr::top_n(3) %$% gene %>% unique 


cite_data = GetAssayData(seurat_isotype,slot = 'data',assay = 'ADT') %>% t()
rna_subset = GetAssayData(seurat_isotype,assay = 'RNA',slot = 'data')[gene_plot, ] %>%
  as.data.frame() %>% 
  t() %>% 
  as.matrix()

# combine into dataframe 
d_r = cbind(seurat_isotype@meta.data, rna_subset) 
d_p = cbind(seurat_isotype@meta.data, cite_data) 
# calculate the median protein expression per cluster
dat_plot_r = d_r %>% 
  dplyr::group_by(wsnn_res.0.8) %>% 
  dplyr::summarize_at(.vars = c(gene_plot), .funs = median) %>% 
  tibble::remove_rownames() %>% 
  tibble::column_to_rownames("wsnn_res.0.8") 

dat_plot_p = d_p %>% 
  dplyr::group_by(wsnn_res.0.8) %>% 
  dplyr::summarize_at(.vars = c(prots), .funs = median) %>% 
  tibble::remove_rownames() %>% 
  tibble::column_to_rownames("wsnn_res.0.8") 

# protein heatmap 
# protein heatmap 
prot_col = circlize::colorRamp2(breaks = seq(-1,25, by = 1), 
                                colors = viridis::viridis(n = 27, option = "B"))
p1 = Heatmap(t(dat_plot_p)[prots, ], 
             name = "protein", 
             col = prot_col, 
             use_raster = T,
             row_names_gp = gpar(color = "black", fontsize = 5)
)
p1


# mRNA heatmap 
mrna = t(dat_plot_r)[gene_plot, ]
rna_col = circlize::colorRamp2(breaks = c(-2,-1,0,1,2), 
                               colors = colorspace::diverge_hsv(n = 5))

zero_variance_rows <- which(apply(mrna, 1, var) == 0)
mrna_cleaned <- mrna[-zero_variance_rows, ]


p2 = Heatmap(t(scale(t(mrna_cleaned))), 
             name = "mRNA", 
             col = rna_col,
             use_raster = T, 
             clustering_method_columns = 'average',
             column_names_gp = gpar(color = "black", fontsize = 7), 
             row_names_gp = gpar(color = "black", fontsize = 5))


p2

# heatmaps
setwd('C:/Users/axi313/Documents/ECCITE_Seq_Katherine/WNN/Heatmap')
tiff("ADT_Heatmap_ALLPROT.tiff", width = 2500, height = 4000, res = 300) # Adjust width, height, and resolution as needed
draw(p1)
dev.off()
tiff("RNA_Heatmap_TopGene.tiff", width = 2000, height = 2000, res = 300) # Adjust width, height, and resolution as needed
draw(p2)
dev.off()

Idents(seurat_isotype)

########################################## ADT Feature Plots and VLN Plots ###############################################
DefaultAssay(seurat_isotype) <- 'ADT'

#Subset HIV+ Twins

seurat_isotype_h <- subset(seurat_isotype, subset = condition == 'HIV+')

features <- seurat_isotype@assays$ADT@counts@Dimnames[[1]]

# Violin Plot

setwd('C:/Users/axi313/Documents/ECCITE_Seq_Katherine/WNN/Violin_Plot/ADT')
for (i in features) {
  vln.pl <- VlnPlot(seurat_isotype_h, features = i, assay = 'ADT', split.by = 'stim')
  ggsave(paste0(i,'_VLNplot_HIVonly.png'),dpi=500, width = 13, vln.pl)
  vln.pl.2 <- VlnPlot(seurat_isotype_h, features = i, assay = 'ADT', split.by = 'stim', group.by = 'predicted.celltype.l2')
  ggsave(paste0(i,'_VLNplot_HIVonly_Azimuth.png'),dpi=500, width = 13, vln.pl.2)
}

# Feature Plots

setwd('C:/Users/axi313/Documents/ECCITE_Seq_Katherine/WNN/Feature_Plot/ADT')
for (i in features) {
  fea.pl <- FeaturePlot(seurat_isotype_h, reduction = 'wnn.umap', features = i)
  ggsave(paste0(i,'_Featureplot_HIVonly.png'),dpi=500, width = 10, fea.pl)
}



###########################
# load(paste0(load.path,"Seuratv5_WNN_Complete.RData"))

################################## Dimplots ############################################
setwd('C:/Users/axi313/Desktop/Analysis/Pediatric/Output/CITEseqV3/WNN/Cluster_Plots/By_Cluster')
DimPlot(subset(seurat_isotype, idents = 0), reduction='wnn.umap',label=T,label.size = 5,repel=T)
for (i in 0:35) {
  x <- DimPlot(subset(seurat_isotype, idents = i), reduction='wnn.umap',label=T,label.size = 5,repel=T)+xlim(-20, 12)+ ylim(-20, 12)
  ggsave(paste0('Cluster_',i,'_Plot.png'),dpi=500, width = 13, x)
}

################################ Old RNA VLN PLOT and Nebulosa ############################################
fe.2 <- c('CD14','FCGR2B','SERPING1','CCR7','CD27','TCF7','CCL5','FCGR3A','PRF1','CD40LG','IRF8','TNFRSF4',
          'CD8A','TNFRSF9','XCL2','CD7','CD8B','NELL2','C1QBP','CD3E','ICOS','IGFBP2','IGFBP4','LDHA',
          'CCND3','MIR155HG','NR4A1','CTLA4','FOXP3','IL2RA','CD19','CD79A','IGHM','EBI3','HLA-DPA1',
          'HLA-DRB1','CTSW','KLRC1','TNFRSF18','CCR4','IRF4','MALAT1','IKZF2','TRDV1','TRGC2',
          'CD3D','CXCR3','GZMK','CCL2','HLA-DRA','SERPINA1','GNLY','NKG7','TIGIT','LTB','MAL','SELL',
          'CCL4L2','CD70','IFNG','IL2RB','KLRD1','TRBC1','HAVCR2','LGALS1','NCAM1','CD36','CD4','IFI30',
          'CXCL8','ITGAX','IL18BP','TNF','TRDV2','TRGV9','FABP5','MT-ND1','MT-ND5','CCL3','IL1B','TNFAIP2',
          'CD40','MS4A1','XCL1','HIST1H4C','LTA','MKI67')

setwd('C:/Users/axi313/Desktop/Analysis/Pediatric/Output/CITEseqV3/WNN/Violin_Plot/RNA_2')

for (i in fe.2) {
  vln.pl <- VlnPlot(seurat_isotype, features = i, assay = 'RNA')
  ggsave(paste0(i,'_VLNplot.png'),dpi=500, width = 13, vln.pl)
}

# Feature Plots

setwd('C:/Users/axi313/Desktop/Analysis/Pediatric/Output/CITEseqV3/WNN/Feature_Plot/RNA_2')
DefaultAssay(seurat_isotype) <- 'RNA'

for (i in fe.2) {
  fea.pl <- FeaturePlot(seurat_isotype, reduction = 'wnn.umap', features = i)
  ggsave(paste0(i,'_Featureplot.png'),dpi=500, width = 10, fea.pl)
}

# Nebulosa
setwd('C:/Users/axi313/Desktop/Analysis/Pediatric/Output/CITEseqV3/WNN/Nebulosa_Plot/RNA_2')
DefaultAssay(seurat_isotype) <- 'RNA'
for (i in fe.2) {
  neb.pl <- plot_density(seurat_isotype, reduction = 'wnn.umap', features = i)
  ggsave(paste0(i,'_Featureplot_Nebulosa.png'),dpi=500, width = 10, neb.pl)
}
