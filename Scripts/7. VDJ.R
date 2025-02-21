
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
library(scRepertoire)
library(igraph)
library(Cairo)

################################################# GLobal Variables ###############################################
#Set dir
setwd("C:/Users/ammas/Documents/ECCITE_Seq_Katherine/VDJ")
load.path <- "C:/Users/ammas/Documents/ECCITE_Seq_Katherine/saved_R_data/"
in.path <- "C:/Users/ammas/Documents/ECCITE_Seq_Katherine/Raw_Data/"
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


#save(combined.TCR,combined.BCR, file = paste0(load.path, "Seuratv5_TCR_BCR_Split.RData"))

#load(paste0(load.path,"TCR_BCR_Split.RData"))

################################################### Basic Clonal Visualizations ################################################

setwd('C:/Users/ammas/Documents/ECCITE_Seq_Katherine/VDJ/TCR/Clonal_Visualizations'
)

clonalQuant(combined.TCR, 
            cloneCall="strict", 
            chain = "both", 
            scale = FALSE)
ggsave('Unique_Clones_Strict_TRAB_raw.png',width=24,height=12)

### Clonal Abundance
clonalAbundance(combined.TCR, 
                cloneCall = "strict", 
                scale = FALSE)
ggsave('Clonal_Abundance_TRAB_raw.png',width=16,height=12)

### Clonal Length
clonalLength(combined.TCR, 
             cloneCall="aa", 
             chain = "both") 
ggsave('Clonal_Length_TRAB_raw.png',width=16,height=12)

clonalHomeostasis(combined.TCR, 
                  cloneCall = "strict",
                  cloneSize = c(Rare = 0.001, Small = 0.01, Medium = 0.1, Large = 0.3, Hyperexpanded =
                                  1))
ggsave('Clonal_Homeostasis_TRAB_scaled.png',width=24,height=12)


######################################## CD3 Composition ##################################
setwd('C:/Users/ammas/Documents/ECCITE_Seq_Katherine/VDJ/TCR/CD3_Composition')

percentAA(combined.TCR, 
          chain = "TRA", 
          aa.length = 20)
ggsave('Percent_AA_TRA.png',width=26,height=24)

positionalEntropy(combined.TCR, 
                  chain = "both", 
                  aa.length = 20)
ggsave('Positional_Entropy_TRAB.png',width=20,height=15)

vizGenes(combined.TCR, 
         x.axis = "TRAV",
         y.axis = NULL,
         plot = "heatmap",  
         scale = TRUE)
ggsave('Heatmap_TRA_V_gene.png',width=12,height=9)

vizGenes(combined.TCR, 
         x.axis = "TRBV",
         y.axis = NULL,
         plot = "heatmap",  
         scale = TRUE)
ggsave('Heatmap_TRB_V_gene.png',width=12,height=9)

vizGenes(combined.TCR, 
         x.axis = "TRAJ",
         y.axis = NULL,
         plot = "heatmap",  
         scale = TRUE)
ggsave('Heatmap_TRA_J_gene.png',width=12,height=9)

vizGenes(combined.TCR, 
         x.axis = "TRBJ",
         y.axis = NULL,
         plot = "heatmap",  
         scale = TRUE)
ggsave('Heatmap_TRB_J_gene.png',width=12,height=9)

percentKmer(combined.TCR, 
            cloneCall = "aa",
            chain = "TRA", 
            motif.length = 3, 
            top.motifs = 25)
ggsave('Heatmap_TRA_kmer.png',width=12,height=9)

percentKmer(combined.TCR, 
            cloneCall = "aa",
            chain = "TRB", 
            motif.length = 3, 
            top.motifs = 25)
ggsave('Heatmap_TRB_kmer.png',width=12,height=9)

############################### Clonal Diversity and Overlap #####################################################
setwd('C:/Users/ammas/Documents/ECCITE_Seq_Katherine/VDJ/TCR/Clonal_Diversity')


clonalDiversity(combined.TCR, 
                cloneCall = "strict")
ggsave('Clonal_Diversity_Strict.png',width=12,height=9)

clonalOverlap(combined.TCR, 
              cloneCall = "strict", 
              method = "raw")
ggsave('Clonal_Overlap_Strict_raw.png',width=24,height=11)

clonalOverlap(combined.TCR, 
              cloneCall = "gene", 
              method = "raw")
ggsave('Clonal_Overlap_Gene_raw.png',width=24,height=11)


#### 
setwd('C:/Users/ammas/Documents/ECCITE_Seq_Katherine/VDJ/TCR/Clonal_Diversity/By_Sample')

clonalCompare(combined.TCR, 
              top.clones = 20, 
              samples = c("P1042MED", "P1042CD3","P1042HIV"), 
              order.by = c("P1042MED", "P1042CD3","P1042HIV"),
              cloneCall="strict", 
              relabel.clones = T,
              graph = "alluvial")
ggsave('Clonal_Comparison_P1042.png',width=15,height=11)

x <- clonalCompare(combined.TCR, 
              top.clones = 20, 
              samples = c("P1042MED", "P1042CD3","P1042HIV"), 
              order.by = c("P1042MED", "P1042CD3","P1042HIV"),
              cloneCall="strict", 
              relabel.clones = T,
              graph = "alluvial",
              exportTable = T)
write.csv(x,'Clonal_Comparison_P1042.csv',row.names = F)


clonalCompare(combined.TCR, 
              top.clones = 20, 
              samples = c("P1043MED", "P1043CD3","P1043HIV"), 
              order.by = c("P1043MED", "P1043CD3","P1043HIV"),
              cloneCall="strict", 
              relabel.clones = T,
              graph = "alluvial")
ggsave('Clonal_Comparison_P1043.png',width=15,height=11)

x <- clonalCompare(combined.TCR, 
                   top.clones = 20, 
                   samples = c("P1043MED", "P1043CD3","P1043HIV"), 
                   order.by = c("P1043MED", "P1043CD3","P1043HIV"),
                   cloneCall="strict", 
                   relabel.clones = T,
                   graph = "alluvial",
                   exportTable = T)
write.csv(x,'Clonal_Comparison_P1043.csv',row.names = F)

clonalCompare(combined.TCR, 
              top.clones = 20, 
              samples = c("JLKMED", "JLKCD3","JLKHIV"), 
              order.by = c("JLKMED", "JLKCD3","JLKHIV"),
              cloneCall="strict", 
              relabel.clones = T,
              graph = "alluvial")
ggsave('Clonal_Comparison_JLK.png',width=15,height=11)

x <- clonalCompare(combined.TCR, 
                   top.clones = 20, 
                   samples = c("JLKMED", "JLKCD3","JLKHIV"), 
                   order.by = c("JLKMED", "JLKCD3","JLKHIV"),
                   cloneCall="strict", 
                   relabel.clones = T,
                   graph = "alluvial",
                   exportTable = T)
write.csv(x,'Clonal_Comparison_JLK.csv',row.names = F)


########################################## Merge Seurat #############################################################
load(paste0(load.path,"Seuratv5_WNN_Complete.RData"))

# TARA All
# Access the cell barcodes (Assuming they are in the column names of the data slot)
barcodes <- rownames(seurat_isotype[[]])

# Use gsub to modify the barcodes, removing everything before and including the third '_'
modified_barcodes <- gsub(".*_.*_.*_(.*)", "\\1", barcodes)
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


##############

library(viridis)  # Load the viridis package

DimPlot_scCustom(seurat.tcr,split.by = 'orig.ident', reduction = 'wnn.umap', group.by = "cloneSize", split_seurat = T, num_col=3) +
  scale_color_viridis(option = "viridis", discrete = TRUE, na.value = "gray")

ggsave('WNN_TCR_Overlay_byPID.png',dpi=500, height = 11, width = 18)

