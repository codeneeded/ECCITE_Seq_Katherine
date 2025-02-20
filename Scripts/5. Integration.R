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


setwd("C:/Users/ammas/Documents/ECCITE_Seq_Katherine/Integration")
load.path <- "C:/Users/ammas/Documents/ECCITE_Seq_Katherine/saved_R_data/"

load(paste0(load.path,'Seuratv5_isotype_Assay3.RData'))

################################################# Integrate RNA #############################################
DefaultAssay(seurat_isotype) <- 'RNA'


################################################# Integrate RNA #############################################

### Split Layers for RNA for each batch
seurat_isotype[["RNA"]] <- split(seurat_isotype[["RNA"]], f = seurat_isotype$orig.ident)

# Standard Processing seurat_isotype
seurat_isotype <- NormalizeData(seurat_isotype)
seurat_isotype <- FindVariableFeatures(seurat_isotype)
seurat_isotype <- ScaleData(seurat_isotype)
seurat_isotype <- RunPCA(seurat_isotype)
seurat_isotype <- FindNeighbors(seurat_isotype, dims = 1:30, reduction = "pca")
seurat_isotype <- FindClusters(seurat_isotype, resolution = 2, cluster.name = "unintegrated_clusters")
seurat_isotype <- RunUMAP(seurat_isotype, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")

### seurat_isotype


seurat_isotype <- IntegrateLayers(
  seurat_isotype,
  method= CCAIntegration,
  orig.reduction = "pca",
  assay = 'RNA',
  new.reduction = "integrated.cca.rna"
)

seurat_isotype <- IntegrateLayers(
  seurat_isotype,
  method= FastMNNIntegration,
  assay = 'RNA',
  new.reduction = "integrated.mnn.rna"
)

### EARTH
seurat_isotype <- FindNeighbors(seurat_isotype, reduction = "integrated.cca.rna", dims = 1:30)
seurat_isotype <- FindClusters(seurat_isotype, resolution = 2, cluster.name = "cca_clusters_rna")
seurat_isotype <- RunUMAP(seurat_isotype, reduction = "integrated.cca.rna", dims = 1:30, reduction.name = "umap.cca.rna")
seurat_isotype <- FindNeighbors(seurat_isotype, reduction = "integrated.mnn.rna", dims = 1:30)
seurat_isotype <- FindClusters(seurat_isotype, resolution = 2, cluster.name = "mnn_clusters_rna")
seurat_isotype <- RunUMAP(seurat_isotype, reduction = "integrated.mnn.rna", dims = 1:30, reduction.name = "umap.mnn.rna")

### PLOT

### TARA ALL
p1 <- DimPlot(
  seurat_isotype,
  reduction = "umap.cca.rna",
  group.by = c("predicted.celltype.l2", "cca_clusters_rna"),
  combine = FALSE, label.size = 2
)
p2 <- DimPlot(
  seurat_isotype,
  reduction = "umap.mnn.rna",
  group.by = c("predicted.celltype.l2", "mnn_clusters_rna"),
  combine = FALSE, label.size = 2
)
# Combine the plots with wrap_plots and assign to a variable
combined_plot <- wrap_plots(c(p1, p2), ncol = 2, byrow = FALSE)

# Save the plot using ggsave
ggsave(
  filename = "RNA_Integration.png",  # Specify the file name and format (e.g., .png, .pdf)
  plot = combined_plot,           # The combined plot to save
  width = 18,                     # Set the width of the saved image
  height = 12,                     # Set the height of the saved image
  dpi = 300                       # Set the resolution for the saved image
)

############ ADT Integration #####################

DefaultAssay(seurat_isotype) <-'ADT'

### Convert assay to v5
seurat_isotype[["ADT"]] <- as(object = seurat_isotype[["ADT"]], Class = "Assay5")

### Split Layers for ADT for each batch
seurat_isotype[["ADT"]] <- split(seurat_isotype[["ADT"]], f = seurat_isotype$orig.ident)

# define proteins to use in clustering (non-isptype controls)
prots <- rownames(seurat_isotype[["ADT"]]$data)
isotype_genes <- c('Mouse-IgG1', 'Mouse-IgG2a', 'Mouse-IgG2b', 'Rat-IgG2b', 'Rat-IgG1', 'Rat-IgG2a', 'Hamster-IgG')
prots <- setdiff(prots, isotype_genes)
VariableFeatures(seurat_isotype) <- prots


# Standard Processing seurat_isotype
seurat_isotype <- ScaleData(seurat_isotype)
seurat_isotype <- RunPCA(seurat_isotype,reduction.name = "apca")
seurat_isotype <- FindNeighbors(seurat_isotype, dims = 1:30, reduction = "apca")
seurat_isotype <- FindClusters(seurat_isotype, resolution = 2, cluster.name = "unintegrated_clusters_adt")
seurat_isotype <- RunUMAP(seurat_isotype, dims = 1:30, reduction = "apca", reduction.name = "umap.unintegrated.adt")

### seurat_isotype
seurat_isotype <- IntegrateLayers(
  seurat_isotype,
  orig.reduction = 'apca',
  features = prots,
  method= FastMNNIntegration,
  assay = 'ADT',
  new.reduction = "integrated.mnn.adt"
)

seurat_isotype <- IntegrateLayers(
  seurat_isotype,
  method= CCAIntegration,
  orig.reduction = "apca",
  features = prots,
  assay = 'ADT',
  new.reduction = "integrated.cca.adt"
)

#### UMAP and Clustering #########
### TARA All
seurat_isotype <- FindNeighbors(seurat_isotype, reduction = "integrated.cca.adt", dims = 1:30)
seurat_isotype <- FindClusters(seurat_isotype, resolution = 2, cluster.name = "cca_clusters_adt")
seurat_isotype <- RunUMAP(seurat_isotype, reduction = "integrated.cca.adt", dims = 1:30, reduction.name = "umap.cca.adt")
seurat_isotype <- FindNeighbors(seurat_isotype, reduction = "integrated.mnn.adt", dims = 1:30)
seurat_isotype <- FindClusters(seurat_isotype, resolution = 2, cluster.name = "mnn_clusters_adt")
seurat_isotype <- RunUMAP(seurat_isotype, reduction = "integrated.mnn.adt", dims = 1:30, reduction.name = "umap.mnn.adt")

### PLOT

### TARA ALL
p1 <- DimPlot(
  seurat_isotype,
  reduction = "umap.cca.adt",
  group.by = c("predicted.celltype.l2", "cca_clusters_adt"),
  combine = FALSE, label.size = 2
)
p2 <- DimPlot(
  seurat_isotype,
  reduction = "umap.mnn.adt",
  group.by = c("predicted.celltype.l2", "mnn_clusters_adt"),
  combine = FALSE, label.size = 2
)
# Combine the plots with wrap_plots and assign to a variable
combined_plot <- wrap_plots(c(p1, p2), ncol = 2, byrow = FALSE)

# Save the plot using ggsave
ggsave(
  filename = "seurat_isotype_ADT_Integration.png",  # Specify the file name and format (e.g., .png, .pdf)
  plot = combined_plot,           # The combined plot to save
  width = 18,                     # Set the width of the saved image
  height = 12,                     # Set the height of the saved image
  dpi = 300                       # Set the resolution for the saved image
)

#### Save to load path

save(seurat_isotype, file = paste0(load.path, "seurat_isotype_RNA_ADT.Rdata"))


############# WSNN #######################################

seurat_isotype <- FindMultiModalNeighbors(
  seurat_isotype, reduction.list = list("integrated.mnn.rna", "integrated.cca.adt"), 
  dims.list = list(1:30, 1:14), modality.weight.name = "wnn.weight"
)

################ UMAP ##########################################

### TARA ALL

seurat_isotype <- RunUMAP(seurat_isotype, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")

seurat_isotype <- FindClusters(seurat_isotype, graph.name = "wsnn", algorithm = 2, resolution = 0.5, cluster.name = 'snn.louvianmlr_0.5')
seurat_isotype <- FindClusters(seurat_isotype, graph.name = "wsnn", algorithm = 2, resolution = 1, cluster.name = 'snn.louvianmlr_1')
seurat_isotype <- FindClusters(seurat_isotype, graph.name = "wsnn", algorithm = 2, resolution = 1.5, cluster.name = 'snn.louvianmlr_1.5')

seurat_isotype <- FindClusters(seurat_isotype, graph.name = "wsnn", algorithm = 3, resolution = 0.5, cluster.name = 'snn.slm_0.5')
seurat_isotype <- FindClusters(seurat_isotype, graph.name = "wsnn", algorithm = 3, resolution = 1, cluster.name = 'snn.slm_1')
seurat_isotype <- FindClusters(seurat_isotype, graph.name = "wsnn", algorithm = 3, resolution = 1.5, cluster.name = 'snn.slm_1.5')

### Convert cluster to numeric/factors
# Convert cluster identities to factors with numeric sorting
# Convert cluster identities to factors with numeric sorting
for (cluster_col in c("snn.louvianmlr_0.5", "snn.louvianmlr_1", "snn.louvianmlr_1.5",
                      "snn.slm_0.5", "snn.slm_1", "snn.slm_1.5")) {
  seurat_isotype[[cluster_col]][, 1] <- factor(
    seurat_isotype[[cluster_col]][, 1],
    levels = sort(as.numeric(levels(seurat_isotype[[cluster_col]][, 1])))
  )
}

####### Plot comparing cluster
# Load required library
library(patchwork)

# TARA ALL

# Azimuth plots (first row)
p1 <- DimPlot_scCustom(
  seurat_isotype,
  reduction = "wnn.umap",
  group.by = "predicted.celltype.l1",
  label = TRUE,
  repel = TRUE,
  label.size = 5
) + ggtitle("Azimuth L1")

p2 <- DimPlot_scCustom(
  seurat_isotype,
  reduction = "wnn.umap",
  group.by = "predicted.celltype.l2",
  label = TRUE,
  repel = TRUE,
  label.size = 5
) + ggtitle("Azimuth L2")

p3 <- DimPlot_scCustom(
  seurat_isotype,
  reduction = "wnn.umap",
  group.by = "predicted.celltype.l3",
  label = TRUE,
  repel = TRUE,
  label.size = 5
) + ggtitle("Azimuth L3")

# SLM plots (second row)
p4 <- DimPlot_scCustom(
  seurat_isotype,
  reduction = "wnn.umap",
  group.by = "snn.slm_0.5",
  label = TRUE,
  repel = TRUE,
  label.size = 5
) + ggtitle("snn.slm_0.5")

p5 <- DimPlot_scCustom(
  seurat_isotype,
  reduction = "wnn.umap",
  group.by = "snn.slm_1",
  label = TRUE,
  repel = TRUE,
  label.size = 5
) + ggtitle("snn.slm_1")

p6 <- DimPlot_scCustom(
  seurat_isotype,
  reduction = "wnn.umap",
  group.by = "snn.slm_1.5",
  label = TRUE,
  repel = TRUE,
  label.size = 5
) + ggtitle("snn.slm_1.5")

# LouvainMLR plots (third row)
p7 <- DimPlot_scCustom(
  seurat_isotype,
  reduction = "wnn.umap",
  group.by = "snn.louvianmlr_0.5",
  label = TRUE,
  repel = TRUE,
  label.size = 5
) + ggtitle("snn.louvianmlr_0.5")

p8 <- DimPlot_scCustom(
  seurat_isotype,
  reduction = "wnn.umap",
  group.by = "snn.louvianmlr_1",
  label = TRUE,
  repel = TRUE,
  label.size = 5
) + ggtitle("snn.louvianmlr_1")

p9 <- DimPlot_scCustom(
  seurat_isotype,
  reduction = "wnn.umap",
  group.by = "snn.louvianmlr_1.5",
  label = TRUE,
  repel = TRUE,
  label.size = 5
) + ggtitle("snn.louvianmlr_1.5")

# Combine all plots into a 3x3 grid
combined_plot <- wrap_plots(p1, p2, p3, p4, p5, p6, p7, p8, p9, ncol = 3, nrow = 3)

# Save the combined plot to a file
ggsave(
  filename = "seurat_isotype_WNN_Clusters.png",
  plot = combined_plot,
  width = 24,  # Adjust width as needed
  height = 17,  # Adjust height as needed
  dpi = 300
)

#### Save to load path

save(seurat_isotype, file = paste0(load.path, "seurat_isotype_WNN.Rdata"))



########## Annotation #############

### Cluster Plots

setwd('C:/Users/axi313/Documents/ECCITE_Seq_Katherine/WNN/Cluster_Plots')

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

setwd('C:/Users/axi313/Documents/ECCITE_Seq_Katherine/WNN/Feature_Plot/Stim_Split')

fe <- c('CD69','IL2RA', 'TNFRSF4', 'FOXP3', 'IFNG', 'PDCD1')
fe_p <- c('CD69','IL2RA', 'TNFRSF4', 'IFNGR1', 'PDCD1')

DefaultAssay(seurat_isotype_h) <- 'ADT'
for (i in fe_p) {
  fea.pl <- FeaturePlot(seurat_isotype_h, reduction = 'wnn.umap', split.by = 'stim', features = i)
  ggsave(paste0(i,'_ADT_Featureplot_HIVonly_bystim.png'),dpi=500, width = 13, fea.pl)
}

DefaultAssay(seurat_isotype_h) <- 'RNA'
for (i in fe) {
  fea.pl <- FeaturePlot(seurat_isotype_h, reduction = 'wnn.umap', split.by = 'stim', features = i)
  ggsave(paste0(i,'_RNA_Featureplot_HIVonly_bystim.png'),dpi=500, width = 13, fea.pl)
}

###### Cluster Distribution by donor
setwd('C:/Users/axi313/Documents/ECCITE_Seq_Katherine/WNN/Cluster_Proportions')

# Step 1: Extract relevant information
predicted_celltypes <- seurat_isotype_h@meta.data$predicted.celltype.l2
donor_ids <- seurat_isotype_h@meta.data$orig.ident

# Step 2: Calculate cell counts per cluster per donor
cluster_counts <- table(donor_ids, predicted_celltypes)

# Step 3: Create a table with raw numbers
raw_numbers_table <- as.data.frame.matrix(cluster_counts)

# Print or save the two tables
print("Raw Numbers:")
print(raw_numbers_table)

# Export raw_numbers_table as a CSV file
write.csv(raw_numbers_table, file = "raw_numbers_table.csv", row.names = TRUE)

# Create a stacked bar plot
raw <- seurat_isotype_h@meta.data %>% 
  ggplot(aes(x=orig.ident, fill=predicted.celltype.l2)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")


# Calculate percentage distribution
percentage_distribution <- prop.table(table(seurat_isotype_h@meta.data$orig.ident, seurat_isotype_h@meta.data$predicted.celltype.l2), margin = 1) * 100

# Convert the percentage distribution table to a data frame
percentage_df <- as.data.frame(percentage_distribution)

# Rename the columns for clarity
colnames(percentage_df) <- c("Donor ID", "Cluster", "Percentage")


# Plot the stacked bar plot with % values as fill
percent <- ggplot(percentage_df, aes(x = `Donor ID`, y = Percentage, fill = Cluster)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  ggtitle("Percentage Distribution of Cell Clusters by Donor")

ggsave('Cluster_Raw_byDonor.png',dpi=500, width = 7, raw)
ggsave('Cluster_Prop_byDonor.png',dpi=500, width = 7, percent)

###########################
# load(paste0(load.path,"Seuratv5_WNN_Complete.RData"))

################################## Differential Expression ############################################
setwd('C:/Users/axi313/Documents/ECCITE_Seq_Katherine/WNN/Differential_Expression/Azimuth')


# Assuming 'seurat_isotype_h' is your Seurat object

# Define the clusters
clusters <- levels(as.factor(seurat_isotype_h@meta.data$predicted.celltype.l2))
clusters
# Create an empty list to store differential expression results
de_results <- list()

# Perform differential expression analysis for each cluster
for (cluster in clusters) {
  # Subset the data for the current cluster
  subset_data <- subset(seurat_isotype_h, subset= predicted.celltype.l2 == cluster)
  
  # Check if both conditions have at least 20 cells
  if (sum(subset_data$stim == "Med") >= 20 & sum(subset_data$stim == "HIV_peptide") >= 20) {
    # Perform differential expression analysis comparing HIV peptide stimulation to media stimulation
    de_result <- FindMarkers(subset_data, ident.1 = "HIV_peptide", ident.2 = "Med", group.by = 'stim', test.use = "MAST")
    
    # Add the result to the list
    de_results[[cluster]] <- de_result
  } else {
    # If the condition is not met, skip the cluster and print a message
    cat("Skipping cluster", cluster, "due to insufficient cells for differential expression analysis.\n")
  }
}


# Save the results as CSV files
for (result_name in names(de_results)) {
  write.csv(de_results[[result_name]], file = paste0(result_name, "_HIV_peptide_VS_Media__DGE_results_Azimuth.csv"), row.names = TRUE)
}
