
#
library(ggplot2)
library(scRepertoire)
library(Seurat)
library(ggraph)
library(scater)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(data.table)

################################################# GLobal Variables ###############################################
#Set dir
setwd("C:/Users/axi313/Documents/ECCITE_Seq_Katherine/VDJ")
load.path <- "C:/Users/axi313/Documents/ECCITE_Seq_Katherine/saved_R_data/"
in.path <- "C:/Users/axi313/Documents/ECCITE_Seq_Katherine/Raw_Data/"
t.path <- "/per_sample_outs/vdj_t/filtered_contig_annotations.csv"
b.path <- "/per_sample_outs/vdj_b/filtered_contig_annotations.csv"

#Filenames

f_names <- c("P1042MED", "P1042CD3", "P1042HIV", "P1043MED", "P1043CD3", "P1043HIV", "JLKMED", "JLKCD3")

colorblind_vector <- colorRampPalette(rev(c("#0D0887FF", "#47039FFF", 
                                            "#7301A8FF", "#9C179EFF", "#BD3786FF", "#D8576BFF",
                                            "#ED7953FF","#FA9E3BFF", "#FDC926FF", "#F0F921FF")))

########################################################## Load Files #################################################

for (name in f_names) {
  # Construct the file paths
  t_file <- paste0(in.path, name, t.path)
  b_file <- paste0(in.path, name, b.path)
  
  # Read the files
  t_data <- read.csv(t_file)
  b_data <- read.csv(b_file)
  
  # Create dynamically named variables
  assign(paste0(name, ".TCR"), t_data)
  assign(paste0(name, ".BCR"), b_data)
}


# Create Contig list

f_names.TCR <- paste(f_names, ".TCR", sep="")
f_names.BCR <- paste(f_names, ".BCR", sep="")
contig_list.TCR <- as.list(mget(f_names.TCR))
contig_list.BCR <- as.list(mget(f_names.BCR))


#Combine For downstream Analysis

combined.TCR <- combineTCR(contig_list.TCR,samples = f_names)

combined.BCR <- combineBCR(contig_list.BCR,samples = f_names)


#Add variables

##TCR
combined.TCR <- addVariable(combined.TCR, variable.name = 'Stim',
                            variables = c('Media', 'CD3','HIV_Peptide', 'Media','CD3','HIV_Peptide', 'Media' ,
                                          'CD3'))

combined.TCR <- addVariable(combined.TCR, variable.name = 'Condition',
                            variables = c('HIV+', 'HIV+', 'HIV+','HIV+', 'HIV+','HIV+','HIV-','HIV-'))

##BCR

combined.BCR <- addVariable(combined.BCR, variable.name = 'Stim',
                            variables = c('Media', 'CD3','HIV_Peptide', 'Media','CD3','HIV_Peptide', 'Media' ,
                                          'CD3'))

combined.BCR <- addVariable(combined.BCR, variable.name = 'Condition',
                            variables = c('HIV+', 'HIV+', 'HIV+','HIV+', 'HIV+','HIV+','HIV-','HIV-'))

list.receptors <- c(combined.TCR, combined.BCR)


save(combined.TCR,combined.BCR, file = paste0(load.path, "Seuratv5_TCR_BCR_Split.RData"))

#load(paste0(load.path,"TCR_BCR_Split.RData"))

######################################################## Visualisations #########################################

t.ch <- clonalHomeostasis(combined.TCR, cloneCall = "strict")
t.cp <- clonalProportion(combined.TCR, cloneCall = "strict")

b.ch <- clonalHomeostasis(combined.BCR, cloneCall = "strict")
b.cp <- clonalProportion(combined.BCR, cloneCall = "strict")

ggsave('TCR_Clonal_Homeostasis_strict.png', width = 15, dpi = 500, t.ch)
ggsave('TCR_Clonal_Proportion_strict.png', width = 15, dpi = 500, t.cp)
ggsave('BCR_Clonal_Homeostasis_strict.png', width = 15, dpi = 500, b.ch)
ggsave('BCR_Clonal_Proportion_strict.png', width = 15, dpi = 500, b.cp)


########################################## Merge Seurat #############################################################
load(paste0(load.path,"Seuratv5_WNN_Complete.RData"))

# Access the cell barcodes (Assuming they are in the column names of the data slot)
barcodes <- rownames(seurat_isotype[[]])

# Use gsub to modify the barcodes, removing everything before and including the fourth '_'
modified_barcodes <- gsub(".*_.*_.*_.*_(.*)", "\\1", barcodes)
modified_barcodes <- paste0(seurat_isotype$orig.ident, "_", modified_barcodes)

# Assign the modified barcodes back to the Seurat object
seurat_isotype <- RenameCells(seurat_isotype, new.names = modified_barcodes)


seurat.tcr <- combineExpression(combined.TCR, 
                                seurat_isotype, 
                                cloneCall="strict",
                                group.by = 'sample',
                                cloneSize = c(Single = 1, Small = 5, Medium = 20, Large = 100, Hyperexpanded =
                                                500),
                                proportion = FALSE)

seurat.bcr <- combineExpression(combined.BCR, 
                                seurat_isotype, 
                                cloneCall="strict",
                                group.by = 'sample',
                                cloneSize = c(Single = 1, Small = 5, Medium = 20, Large = 100, Hyperexpanded =
                                                500),
                                proportion = FALSE)



######################################################## Hyper Expansion Plots ##################################

cols <- c("Hyperexpanded (100 < X <= 500)"="#F0F921","Large (20 < X <= 100)" = "#F69441",
          "Medium (5 < X <= 20)" = "#CA4778" ,"Small (1 < X <= 5)" = "#7D06A5",
          "Single (0 < X <= 1)" = "#0D0887")

## first ordering the Clone Size as a factor, this prevents the coloring from being in alphabetical order. 
slot(seurat.tcr, "meta.data")$cloneSize <- factor(slot(seurat.tcr, "meta.data")$cloneSize, 
                                                  levels = c("Hyperexpanded (100 < X <= 500)", 
                                                             "Large (20 < X <= 100)", 
                                                             "Medium (5 < X <= 20)", 
                                                             "Small (1 < X <= 5)", 
                                                             "Single (0 < X <= 1)", NA))
slot(seurat.bcr, "meta.data")$cloneSize <- factor(slot(seurat.bcr, "meta.data")$cloneSize, 
                                                  levels = c("Hyperexpanded (100 < X <= 500)", 
                                                             "Large (20 < X <= 100)", 
                                                             "Medium (5 < X <= 20)", 
                                                             "Small (1 < X <= 5)", 
                                                             "Single (0 < X <= 1)", NA))


### TCR
setwd("C:/Users/axi313/Documents/ECCITE_Seq_Katherine/VDJ/Seurat_Plots")

tcr.condition <- DimPlot(seurat.tcr, group.by = "cloneSize", split.by = 'condition', reduction = 'wnn.umap') +
  scale_color_manual(values=cols) + ggtitle('Expanded TCR Clone Size')
tcr.stim <- DimPlot(seurat.tcr, group.by = "cloneSize", split.by = 'stim', reduction = 'wnn.umap') +
  scale_color_manual(values=cols) + ggtitle('Expanded TCR Clone Size')


ggsave('TCR_Seurat_Condition.png', width = 12, dpi = 500, tcr.condition)
ggsave('TCR_Seurat_Stim.png', width = 12, dpi = 500, tcr.stim)
#ggsave('TCR_Seurat_Sample.png', width = 18, dpi = 500, tcr.sample)


### BCR
bcr.condition <- DimPlot(seurat.bcr, group.by = "cloneSize", split.by = 'condition', reduction = 'wnn.umap') +
  scale_color_manual(values=cols) + ggtitle('Expanded BCR Clone Size')
bcr.stim <- DimPlot(seurat.bcr, group.by = "cloneSize", split.by = 'stim', reduction = 'wnn.umap') +
  scale_color_manual(values=cols) + ggtitle('Expanded BCR Clone Size')
#bcr.sample <- DimPlot(seurat.bcr, group.by = "cloneSize", split.by = 'sample', reduction = 'wnn.umap') +
#  scale_color_manual(values=cols) + ggtitle('Expanded BCR Clone Size')

ggsave('BCR_Seurat_Condition.png', width = 12, dpi = 500, bcr.condition)
ggsave('BCR_Seurat_Stim.png', width = 12, dpi = 500, bcr.stim)
#ggsave('BCR_Seurat_Sample.png', width = 18, dpi = 500, bcr.sample)

### By Sample 
setwd('C:/Users/axi313/Documents/ECCITE_Seq_Katherine/VDJ/Seurat_Plots/By_Sample')

for (i in f_names) {
  plot <- DimPlot(subset(seurat.tcr, subset= `orig.ident`== i), group.by = "cloneSize", reduction = 'wnn.umap') +
    scale_color_manual(values=cols) + ggtitle(paste0(i,' Expanded TCR Clone Size'))
  ggsave(paste0(i,'_TCR_Seurat.png'),width = 9,dpi=500,plot)
}

for (i in f_names) {
  plot <- DimPlot(subset(seurat.bcr, subset= `orig.ident`== i), group.by = "cloneSize", reduction = 'wnn.umap') +
    scale_color_manual(values=cols) + ggtitle(paste0(i,' Expanded BCR Clone Size'))
  ggsave(paste0(i,'_BCR_Seurat.png'),width = 9,dpi=500,plot)
}


############################ Extract output table of most clonal expanded Clone Sizes#############################

### TCR
setwd('C:/Users/axi313/Documents/ECCITE_Seq_Katherine/VDJ/Clone_Table')

seurat.tcr.aa <- seurat.tcr@meta.data %>%
  group_by(CTaa, predicted.celltype.l2, orig.ident) %>%
  count() %>%
  na.omit()

seurat.tcr.gene <- seurat.tcr@meta.data %>%
  group_by(CTgene, predicted.celltype.l2, orig.ident) %>%
  count() %>%
  na.omit()

seurat.tcr.strict <- seurat.tcr@meta.data %>%
  group_by(CTstrict, predicted.celltype.l2, orig.ident) %>%
  count() %>%
  na.omit()


write.csv(seurat.tcr.aa, file = 'TCR_Ranked_Clonal_Expansion_aa_Azimuth.csv', row.names = F)
write.csv(seurat.tcr.gene, file = 'TCR_Ranked_Clonal_Expansion_gene_Azimuth.csv', row.names = F)
write.csv(seurat.tcr.strict, file = 'TCR_Ranked_Clonal_Expansion_strict_Azimuth.csv', row.names = F)

### BCR

seurat.bcr.aa <- seurat.bcr@meta.data %>%
  group_by(CTaa, predicted.celltype.l2, orig.ident) %>%
  count() %>%
  na.omit()

seurat.bcr.gene <- seurat.bcr@meta.data %>%
  group_by(CTgene, predicted.celltype.l2, orig.ident) %>%
  count() %>%
  na.omit()

seurat.bcr.strict <- seurat.bcr@meta.data %>%
  group_by(CTstrict, predicted.celltype.l2, orig.ident) %>%
  count() %>%
  na.omit()


write.csv(seurat.bcr.aa, file = 'BCR_Ranked_Clonal_Expansion_aa_Azimuth.csv', row.names = F)
write.csv(seurat.bcr.gene, file = 'BCR_Ranked_Clonal_Expansion_gene_Azimuth.csv', row.names = F)
write.csv(seurat.bcr.strict, file = 'BCR_Ranked_Clonal_Expansion_strict_Azimuth.csv', row.names = F)

save(seurat.tcr,seurat.bcr, file = paste0(load.path, "Seuratv5_TCR_BCR_VDJ.RData"))
