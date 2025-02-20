#Load Required Libraries
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


setwd("C:/Users/ammas/Documents/ECCITE_Seq_Katherine/WNN")
load.path <- "C:/Users/ammas/Documents/ECCITE_Seq_Katherine/saved_R_data/"

load(paste0(load.path,'seurat_isotype_WNN.Rdata'))


### Cluster Plots

setwd('C:/Users/ammas/Documents/ECCITE_Seq_Katherine/WNN/Cluster_Plots')

DimPlot_scCustom(
  seurat_isotype,
  reduction = "wnn.umap",
  group.by = "snn.louvianmlr_1",
  split.by = 'stim',
  label = TRUE,
  repel = TRUE,
  label.size = 5
)

ggsave('WNN_ClusterPlot_res0.8_bystim.png',dpi=500, width = 13)


DimPlot_scCustom(
  seurat_isotype,
  reduction = "wnn.umap",
  group.by = "snn.louvianmlr_1",
  split.by = 'PID',
  label = TRUE,
  repel = TRUE,
  label.size = 5
)
ggsave('WNN_ClusterPlot_res0.8_byPID.png',dpi=500, width = 13)

DimPlot_scCustom(
  seurat_isotype,
  reduction = "wnn.umap",
  group.by = "snn.louvianmlr_1",
  split.by = 'condition',
  label = TRUE,
  repel = TRUE,
  label.size = 5
)
ggsave('WNN_ClusterPlot_res0.8_bycondition.png',dpi=500, width = 13)

Idents(seurat_isotype) <- "snn.louvianmlr_1"
Idents(seurat_isotype)
########################################## Feature Plots and VLN Plots ###############################################
### ADT

DefaultAssay(seurat_isotype) <- 'ADT'

#Subset HIV+ Twins

features <- rownames(seurat_isotype[["ADT"]]$data)

# VLN Plots 

# Clusters

setwd('C:/Users/ammas/Documents/ECCITE_Seq_Katherine/WNN/Violin_Plot/ADT')

for (i in features) {
  vln.pl <- VlnPlot_scCustom(seurat_isotype, features = i, assay = 'ADT',  plot_median = T) & NoLegend()
  ggsave(paste0(i,'_VLNplot.png'),dpi=500, width = 15, vln.pl)
  vln.pl.2 <- VlnPlot_scCustom(seurat_isotype, features = i, assay = 'ADT',  split.by = 'stim') 
  ggsave(paste0(i,'_VLNplot_by_stim.png'),dpi=500, width = 15, vln.pl.2)
  vln.pl.3 <- VlnPlot_scCustom(seurat_isotype, features = i, assay = 'ADT',  split.by = 'condition') 
  ggsave(paste0(i,'_VLNplot_by_condition.png'),dpi=500, width = 15, vln.pl.3)
  vln.pl.4 <- VlnPlot_scCustom(seurat_isotype, features = i, assay = 'ADT',  split.by = 'PID')
  ggsave(paste0(i,'_VLNplot_by_PID.png'),dpi=500, width = 15, vln.pl.4)
}

# Azimuth
setwd('C:/Users/ammas/Documents/ECCITE_Seq_Katherine/WNN/Violin_Plot/ADT_Azimuth')

for (i in features) {
  vln.pl <- VlnPlot_scCustom(seurat_isotype, features = i, assay = 'ADT', group.by = 'predicted.celltype.l2', plot_median = T) & NoLegend()
  ggsave(paste0(i,'_VLNplot.png'),dpi=500, width = 15, vln.pl)
  vln.pl.2 <- VlnPlot_scCustom(seurat_isotype, features = i, assay = 'ADT', group.by = 'predicted.celltype.l2', split.by = 'stim') 
  ggsave(paste0(i,'_VLNplot_by_stim.png'),dpi=500, width = 15, vln.pl.2)
  vln.pl.3 <- VlnPlot_scCustom(seurat_isotype, features = i, assay = 'ADT', group.by = 'predicted.celltype.l2', split.by = 'condition') 
  ggsave(paste0(i,'_VLNplot_by_condition.png'),dpi=500, width = 15, vln.pl.3)
  vln.pl.4 <- VlnPlot_scCustom(seurat_isotype, features = i, assay = 'ADT', group.by = 'predicted.celltype.l2', split.by = 'PID')
  ggsave(paste0(i,'_VLNplot_by_PID.png'),dpi=500, width = 15, vln.pl.4)
}


# Feature Plots

setwd('C:/Users/ammas/Documents/ECCITE_Seq_Katherine/WNN/Feature_Plot/ADT')
for (i in features) {
  fea.pl <- FeaturePlot_scCustom(seurat_isotype, reduction = 'wnn.umap', features = i)
  ggsave(paste0(i,'_Featureplot.png'),dpi=500, width = 6, fea.pl)
  fea.pl.2 <- FeaturePlot_scCustom(seurat_isotype, reduction = 'wnn.umap', features = i, split.by = "orig.ident",
                                   num_columns = 4)
  ggsave(paste0(i,'_Featureplot_by_PID.png'),dpi=500, height=14, width = 21, fea.pl.2)
  fea.pl.3 <- FeaturePlot_scCustom(seurat_isotype, reduction = 'wnn.umap', features = i, split.by = "stim",
                                   num_columns = 3)
  ggsave(paste0(i,'_Featureplot_by_Stim.png'),dpi=500, width = 16, fea.pl.3)
  fea.pl.4 <- FeaturePlot_scCustom(seurat_isotype, reduction = 'wnn.umap', features = i, split.by = "condition",
                                   num_columns = 2)
  ggsave(paste0(i,'_Featureplot_by_Condition.png'),dpi=500, width = 11, fea.pl.4)
}

### Average Expression, Heatmap and Dotplot
setwd('C:/Users/ammas/Documents/ECCITE_Seq_Katherine/WNN/Average_Expression/ADT')

ADT_AE_Azimuth <- AverageExpression(seurat_isotype,assays = 'ADT', group.by = c("predicted.celltype.l2","orig.ident"))
write.csv(ADT_AE_Azimuth$ADT,"Average_Expression_ADT_split_by_sample_Azimuth_Cluster.csv")

ADT_AE_Clust <- AverageExpression(seurat_isotype,assays = 'ADT', group.by = c("snn.louvianmlr_1","orig.ident"))
write.csv(ADT_AE_Azimuth$ADT,"Average_Expression_ADT_split_by_sample_snn.louvianmlr_1.csv")

adt_dp <- Clustered_DotPlot(seurat_object = seurat_isotype, features = features, group.by= 'predicted.celltype.l2', 
                  plot_km_elbow=F, plot_padding = TRUE)

CairoPNG("AdT_Dotplot_Clustered.png", width = 1700, height = 4000, dpi=200)

# Draw the heatmap
draw(adt_dp)

# Close the device and save the PNG
dev.off()


### RNA

DefaultAssay(seurat_isotype) <- 'RNA'

features <- c('CD14','FCGR2B','SERPING1','CCR7','CD27','TCF7','CCL5','FCGR3A','PRF1','CD40LG','IRF8','TNFRSF4',
              'CD8A','TNFRSF9','XCL2','CD7','CD8B','NELL2','C1QBP','CD3E','ICOS','IGFBP2','IGFBP4','LDHA',
              'CCND3','MIR155HG','NR4A1','CTLA4','FOXP3','IL2RA','CD19','CD79A','IGHM','EBI3','HLA-DPA1',
              'HLA-DRB1','CTSW','KLRC1','TNFRSF18','CCR4','IRF4','MALAT1','IKZF2','TRDV1','TRGC2',
              'CD3D','CXCR3','GZMK','CCL2','HLA-DRA','SERPINA1','GNLY','NKG7','TIGIT','LTB','MAL','SELL',
              'CCL4L2','CD70','IFNG','IL2RB','KLRD1','TRBC1','HAVCR2','LGALS1','NCAM1','CD36','CD4','IFI30',
              'CXCL8','ITGAX','IL18BP','TNF','TRDV2','TRGV9','FABP5','MT-ND1','MT-ND5','CCL3','IL1B','TNFAIP2',
              'CD40','MS4A1','XCL1','HIST1H4C','LTA','MKI67')
# VLN Plots 

# Clusters

setwd('C:/Users/ammas/Documents/ECCITE_Seq_Katherine/WNN/Violin_Plot/RNA')

for (i in features) {
  vln.pl <- VlnPlot_scCustom(seurat_isotype, features = i, assay = 'RNA',  plot_median = T) & NoLegend()
  ggsave(paste0(i,'_VLNplot.png'),dpi=500, width = 15, vln.pl)
  vln.pl.2 <- VlnPlot_scCustom(seurat_isotype, features = i, assay = 'RNA',  split.by = 'stim') 
  ggsave(paste0(i,'_VLNplot_by_stim.png'),dpi=500, width = 15, vln.pl.2)
  vln.pl.3 <- VlnPlot_scCustom(seurat_isotype, features = i, assay = 'RNA',  split.by = 'condition') 
  ggsave(paste0(i,'_VLNplot_by_condition.png'),dpi=500, width = 15, vln.pl.3)
  vln.pl.4 <- VlnPlot_scCustom(seurat_isotype, features = i, assay = 'RNA',  split.by = 'PID')
  ggsave(paste0(i,'_VLNplot_by_PID.png'),dpi=500, width = 15, vln.pl.4)
}

# Azimuth
setwd('C:/Users/ammas/Documents/ECCITE_Seq_Katherine/WNN/Violin_Plot/RNA_Azimuth')

for (i in features) {
  vln.pl <- VlnPlot_scCustom(seurat_isotype, features = i, assay = 'RNA', group.by = 'predicted.celltype.l2', plot_median = T) & NoLegend()
  ggsave(paste0(i,'_VLNplot.png'),dpi=500, width = 15, vln.pl)
  vln.pl.2 <- VlnPlot_scCustom(seurat_isotype, features = i, assay = 'RNA', group.by = 'predicted.celltype.l2', split.by = 'stim') 
  ggsave(paste0(i,'_VLNplot_by_stim.png'),dpi=500, width = 15, vln.pl.2)
  vln.pl.3 <- VlnPlot_scCustom(seurat_isotype, features = i, assay = 'RNA', group.by = 'predicted.celltype.l2', split.by = 'condition') 
  ggsave(paste0(i,'_VLNplot_by_condition.png'),dpi=500, width = 15, vln.pl.3)
  vln.pl.4 <- VlnPlot_scCustom(seurat_isotype, features = i, assay = 'RNA', group.by = 'predicted.celltype.l2', split.by = 'PID')
  ggsave(paste0(i,'_VLNplot_by_PID.png'),dpi=500, width = 15, vln.pl.4)
}


# Feature Plots

setwd('C:/Users/ammas/Documents/ECCITE_Seq_Katherine/WNN/Feature_Plot/RNA')
for (i in features) {
  fea.pl <- FeaturePlot_scCustom(seurat_isotype, reduction = 'wnn.umap', features = i)
  ggsave(paste0(i,'_Featureplot.png'),dpi=500, width = 6, fea.pl)
  fea.pl.2 <- FeaturePlot_scCustom(seurat_isotype, reduction = 'wnn.umap', features = i, split.by = "orig.ident",
                                   num_columns = 4)
  ggsave(paste0(i,'_Featureplot_by_PID.png'),dpi=500, height=14, width = 21, fea.pl.2)
  fea.pl.3 <- FeaturePlot_scCustom(seurat_isotype, reduction = 'wnn.umap', features = i, split.by = "stim",
                                   num_columns = 3)
  ggsave(paste0(i,'_Featureplot_by_Stim.png'),dpi=500, width = 16, fea.pl.3)
  fea.pl.4 <- FeaturePlot_scCustom(seurat_isotype, reduction = 'wnn.umap', features = i, split.by = "condition",
                                   num_columns = 2)
  ggsave(paste0(i,'_Featureplot_by_Condition.png'),dpi=500, width = 11, fea.pl.4)
}

### Average Expression, Heatmap and Dotplot
setwd('C:/Users/ammas/Documents/ECCITE_Seq_Katherine/WNN/Average_Expression/RNA')

RNA_AE_Azimuth <- AverageExpression(seurat_isotype,assays = 'RNA', group.by = c("predicted.celltype.l2","orig.ident"))
write.csv(RNA_AE_Azimuth$RNA,"Average_Expression_RNA_split_by_sample_Azimuth_Cluster.csv")

RNA_AE_Clust <- AverageExpression(seurat_isotype,assays = 'RNA', group.by = c("snn.louvianmlr_1","orig.ident"))
write.csv(RNA_AE_Clust$RNA,"Average_Expression_RNA_split_by_sample_snn.louvianmlr_1.csv")

RNA_dp <- Clustered_DotPlot(seurat_object = seurat_isotype, features = features, group.by= 'predicted.celltype.l2', 
                            plot_km_elbow=F, plot_padding = TRUE)

CairoPNG("RNA_Dotplot_Clustered.png", width = 1700, height = 4000, dpi=200)

# Draw the heatmap
draw(RNA_dp)

# Close the device and save the PNG
dev.off()

###
### Cluster Stats
setwd('C:/Users/ammas/Documents/ECCITE_Seq_Katherine/WNN/Cluster_Stats')

cluster_stats_ADT_sample <- Cluster_Stats_All_Samples(seurat_isotype, group_by_var = 'orig.ident')
write.csv(cluster_stats_ADT_sample,"Cluster_Stats_by_PID.csv", row.names = F)

cluster_stats_ADT_condition <- Cluster_Stats_All_Samples(seurat_isotype, group_by_var = 'condition')
write.csv(cluster_stats_ADT_condition,"Cluster_Stats_by_Condition.csv", row.names = F)

cluster_stats_ADT_stim <- Cluster_Stats_All_Samples(seurat_isotype, group_by_var = 'stim')
write.csv(cluster_stats_ADT_stim,"Cluster_Stats_by_Stim.csv", row.names = F)

cluster_stats_ADT_Azimuth <- Cluster_Stats_All_Samples(seurat_isotype, group_by_var = 'predicted.celltype.l2')
write.csv(cluster_stats_ADT_Azimuth,"Cluster_Stats_by_Azimuth.csv", row.names = F)

