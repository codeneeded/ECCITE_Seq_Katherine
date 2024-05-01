###Citeseq Pipeline
#Cite seurat and ds packages
#Load Required Libraries
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(hdf5r)
library(dsb)
library(data.table)
library(ggplot2)

###############################################Path + GLobal Variables ##################################################
##set path to load data
setwd("C:/Users/axi313/Documents/ECCITE_Seq_Katherine/Preprocessing")
in.path <- "C:/Users/axi313/Documents/ECCITE_Seq_Katherine/Raw_Data/"
out.path <- "C:/Users/axi313/Documents/ECCITE_Seq_Katherine/saved_R_data/"
r.str <- "/multi/count/raw_feature_bc_matrix.h5"
f.str <- "/per_sample_outs/count/sample_filtered_feature_bc_matrix.h5"

### 
##Load 10x data as String names

f_names <- c("P1042MED", "P1042CD3", "P1042HIV", "P1043MED", "P1043CD3", "P1043HIV", "JLKMED", "JLKCD3")

####################################### Load Data + Basic Pre-processing ##############################################################
for (i in f_names) {
  #Read in Raw and filtered data from 10x output
  raw <- Read10X_h5(paste0(in.path,i,r.str))
  cells<- Read10X_h5(paste0(in.path,i,f.str))
  raw$`Antibody Capture`@Dimnames[[1]] <-  gsub(pattern = "_TotalSeqC", replacement = "", 
                                                x = raw$`Antibody Capture`@Dimnames[[1]])
  cells$`Antibody Capture`@Dimnames[[1]] <-  gsub(pattern = "_TotalSeqC", replacement = "", 
                                                  x = cells$`Antibody Capture`@Dimnames[[1]])
  # define a vector of cell-containing barcodes and remove them from unfiltered data 
  stained_cells <- colnames(cells$`Gene Expression`)
  background <- setdiff(colnames(raw$`Gene Expression`), stained_cells)
  
  # split the data into separate matrices per assay 
  prot <- raw$`Antibody Capture`
  rna <- raw$`Gene Expression`
  
  # create metadata of droplet QC stats used in standard scRNAseq processing
  rna.size <- log10(Matrix::colSums(rna))
  prot.size <- log10(Matrix::colSums(prot))
  nCount_RNA <- Matrix::colSums(rna) #Molecules per cell
  nCount_ADT <- Matrix::colSums(prot) #Molecules per cell
  nFeature_RNA <- Matrix::colSums(rna > 0) #genes per cell
  nFeature_ADT <- Matrix::colSums(prot > 0) #proteins per cell
  mtgene <- grep(pattern = "^MT-", rownames(rna), value = TRUE)
  mt.prop <- Matrix::colSums(rna[mtgene, ]) / Matrix::colSums(rna)
  md <- as.data.frame(cbind( nCount_RNA, nFeature_RNA, nCount_ADT, nFeature_ADT, rna.size, prot.size,mt.prop))
  
  # add indicator for barcodes Cell Ranger called as cells
  md$drop.class <- ifelse(rownames(md) %in% stained_cells, 'cell', 'background')
  
  # remove barcodes with no evidence of capture in the experiment
  md <- md[md$rna.size > 0 & md$prot.size > 0, ]
  
  ##FIX PROTEIN NAMES
  prot@Dimnames[[1]][prot@Dimnames[[1]] == "human(TCR)"] <- "TCR-AB"
  prot@Dimnames[[1]][prot@Dimnames[[1]] == "Human-TCR "] <- "TCR-vA7.2"
  prot@Dimnames[[1]][prot@Dimnames[[1]] == "Human-TCR"] <- "TCR-vD2"
  prot@Dimnames[[1]][prot@Dimnames[[1]] == "PTPRC"] <- "CD45RA"
  prot@Dimnames[[1]][prot@Dimnames[[1]] == "PTPRC.1"] <- "CD45RO"
  prot@Dimnames[[1]][prot@Dimnames[[1]] == "PTPRC.2"] <- "CD45"
  
  #Assign to final output
  assign(paste0(i,'.md'),md)
  assign(paste0(i,'.prot'),prot)
  assign(paste0(i,'.rna'),rna)
}
###  Droplet Settings

for (i in f_names) {
  #Set metadata object and protein object to prevent constant eval parse calling  
  md <- eval(parse(text= paste0(i,'.md')))
  prot <- eval(parse(text = paste0(i,'.prot')))
  rna <- eval(parse(text = paste0(i,'.rna')))
  
  #Output Plot for detected genes vs the protein library size for cells vs background drops
  png(paste0('Droplet_Thresholds/',i,'_genevsprotlibsize.png'),width = 800, height = 600)
  p <- ggplot(md, aes(x = log10(nFeature_RNA), y = prot.size )) + 
    theme_bw() + 
    geom_bin2d(bins = 300) + 
    scale_fill_viridis_c(option = "C") + 
    facet_wrap(~drop.class)
  print(p)
  dev.off()
}

ggplot(P1042CD3.md, aes(x = log10(nFeature_RNA), y = prot.size)) + 
  stat_summary_2d(aes(z = mt.prop, fill = ..value..), fun = mean, bins = 300) +
  scale_fill_viridis_c(option = "C", na.value = "white") + 
  facet_wrap(~drop.class) +
  theme_bw()

############################################### Droplet Thresholds for DSB Norm###################################################

# P1042MED
md <- P1042MED.md
P1042MED.background_drops = rownames(
  md[ md$prot.size > 1 & 
        md$prot.size < 4 & 
        md$rna.size < 2.5, ]
)

# P1042CD3
md <- P1042CD3.md
P1042CD3.background_drops = rownames(
  md[ md$prot.size > 1 & 
        md$prot.size < 4.5 & 
        md$rna.size < 2.5, ]
)

# P1042HIV
md <- P1042HIV.md
P1042HIV.background_drops = rownames(
  md[ md$prot.size > 1 & 
        md$prot.size < 4.5 & 
        md$rna.size < 2.5, ]
)

# P1043MED
md <- P1043MED.md
P1043MED.background_drops = rownames(
  md[ md$prot.size > 1 & 
        md$prot.size < 4.5 & 
        md$rna.size < 2.5, ]
)

# P1043CD3
md <- P1043CD3.md
P1043CD3.background_drops = rownames(
  md[ md$prot.size > 1 & 
        md$prot.size < 4.5 & 
        md$rna.size < 2.5, ]
)

# P1043HIV
md <- P1043HIV.md
P1043HIV.background_drops = rownames(
  md[ md$prot.size > 1 & 
        md$prot.size < 4.5 & 
        md$rna.size < 2.5, ]
)

# JLKMED
md <- JLKMED.md
JLKMED.background_drops = rownames(
  md[ md$prot.size > 1 & 
        md$prot.size < 4 & 
        md$rna.size < 2.5, ]
)

# JLKCD3
md <- JLKCD3.md
JLKCD3.background_drops = rownames(
  md[ md$prot.size > 1 & 
        md$prot.size < 4 & 
        md$rna.size < 2.5, ]
)

###################################Initial QC + DSB Normalization###################################################


for (i in f_names) {
  #Set metadata object and protein object to prevent constant eval parse calling  
  md <- eval(parse(text= paste0(i,'.md')))
  prot <- eval(parse(text = paste0(i,'.prot')))
  rna <- eval(parse(text = paste0(i,'.rna')))
  background_drops <- eval(parse(text = paste0(i,'.background_drops')))
  
  background.adt.mtx = as.matrix(prot[ , background_drops])
  
  cellmd = md[md$drop.class == 'cell', ]  #Define Cell Metadata on only cells not background  
  
  # filter drops with + / - 3 median absolute deviations from the median library size
  rna.mult = (3*mad(cellmd$rna.size))
  prot.mult = (3*mad(cellmd$prot.size))
  rna.lower = median(cellmd$rna.size) - rna.mult
  rna.upper = median(cellmd$rna.size) + rna.mult
  prot.lower = median(cellmd$prot.size) - prot.mult
  prot.upper = median(cellmd$prot.size) + prot.mult
  
  # filter rows based on droplet qualty control metrics
  qc_cells = rownames(
    cellmd[cellmd$prot.size > prot.lower & 
             cellmd$prot.size < prot.upper & 
             cellmd$rna.size > rna.lower & 
             cellmd$rna.size < rna.upper & 
             cellmd$mt.prop < 0.25, ]
  )
  
  # Output thresholds for quality control metrics as in any standard scRNAseq analysis
  png(paste0(i,'_qc_thresholds.png'),width = 800, height = 600)
  plot_aes = list(theme_bw(), geom_point(shape = 21 , stroke = 0, size = 0.7), scale_fill_viridis_c(option = "C"))
  p1 = ggplot(cellmd, aes(x = rna.size )) + geom_histogram(bins = 50) + theme_bw() + xlab("log10 RNA library size")
  p2 = ggplot(cellmd, aes(x = mt.prop)) + geom_histogram(bins = 50) + theme_bw() + xlab("mitochondrial read proportion")
  p3 = ggplot(cellmd, aes(x = log10(nFeature_RNA), y = rna.size, fill = mt.prop )) + plot_aes
  p4 = ggplot(cellmd, aes(x = nFeature_RNA, y = prot.size, fill = mt.prop )) + plot_aes
  print (p1+p2+p3+p4)
  dev.off()
  cell.adt.raw = as.matrix(prot[ , qc_cells])
  cell.rna.raw = rna[ ,qc_cells]
  cellmd = cellmd[qc_cells, ]
  
  #Proteins without Staining
  pm = sort(apply(cell.adt.raw, 1, max))
  pm2 = apply(background.adt.mtx, 1, max)
  head(pm2)
  
  #Assign it to your final output
  assign(paste0(i,'.cell.adt.raw'),cell.adt.raw)
  assign(paste0(i,'.cell.rna.raw'),cell.rna.raw)
  assign(paste0(i,'.background.adt.mtx'),background.adt.mtx)
  assign(paste0(i,'.cellmd'),cellmd)
  assign(paste0(i,'.pm'),pm)
}

### Look at protien_id without staining deside to remove or not.


for (i in f_names) {
  name <- eval(parse(text = paste0(i,'.pm')))
  x <- as.data.table(name,keep.rownames = T)
  x[,Origin:=..i]
  assign(paste0(i,".dt"),x)
}


# Check if you need to remove proteins without staining
#https://www.rdocumentation.org/packages/dsb/versions/0.3.0
#prot.expres.total <- rbindlist(adt.list)

#In this case we do not

### DSB Normalisation
#Set isotype control

isotype.controls <- c('Mouse-IgG1', 'Mouse-IgG2a','Mouse-IgG2b', 'Rat-IgG2b'
                      ,'Rat-IgG1', 'Rat-IgG2a', 'Hamster-IgG')


#normalize protein data for the cell containing droplets with the dsb method. 

for (i in f_names) {
  cell.adt.raw <- eval(parse(text = paste0(i,'.cell.adt.raw')))
  background.adt.mtx <- eval(parse(text = paste0(i,'.background.adt.mtx')))
  # normalize and denoise with dsb with 
  cells.dsb.norm = DSBNormalizeProtein(
    cell_protein_matrix = cell.adt.raw, 
    empty_drop_matrix = background.adt.mtx, 
    denoise.counts = TRUE, 
    use.isotype.control = TRUE, 
    isotype.control.name.vec = isotype.controls
  )
  cells.dsb.norm <- Matrix(as.matrix(cells.dsb.norm),sparse=TRUE)
  assign(paste0(i,'.cells.dsb.norm'),cells.dsb.norm)
}


########################################Create + Merge Seurat Object (norm dsb+RNA)########################################

# Create Seurat Object

for (i in f_names) {
  cellmd <- eval(parse(text= paste0(i,'.cellmd')))
  cell.adt.raw <- eval(parse(text= paste0(i,'.cell.adt.raw')))
  cells.dsb.norm <- eval(parse(text= paste0(i,'.cells.dsb.norm')))
  cell.rna.raw <- eval(parse(text= paste0(i,'.cell.rna.raw')))
  
  # integrating with Seurat
  stopifnot(isTRUE(all.equal(rownames(cellmd), colnames(cell.adt.raw))))
  stopifnot(isTRUE(all.equal(rownames(cellmd), colnames(cell.rna.raw))))
  
  # create Seurat object note: min.cells is a gene filter, not a cell filter
  s = Seurat::CreateSeuratObject(counts = cell.rna.raw, 
                                 meta.data = cellmd,
                                 assay = "RNA", 
                                 min.cells = 20)
  
  # add dsb normalized matrix "dsb_norm_prot" to the "CITE" assay data slot
  s[["ADT"]] <- CreateAssayObject(data = cells.dsb.norm)
  s$orig.ident <- i
  #Assign
  assign(paste0(i,'.seurat'),s)
}

#Merge Seurat Objects
merged_seurat <- merge(x = P1042MED.seurat, 
                       y = c(P1042CD3.seurat, P1042HIV.seurat, P1043MED.seurat,
                             P1043CD3.seurat, P1043HIV.seurat, JLKMED.seurat, JLKCD3.seurat), 
                       add.cell.id = c("P1042_HIV+_Med", "P1042_HIV+_CD3", "P1042_HIV+_HIV", "P1043_HIV+_Med",
                                       "P1043_HIV+_CD3", "P1043_HIV+_HIV","JLK_HIV-_Med","JLK_HIV-_CD3"))
###################################################### Edit Seurat Metadata ###############################################

## Rename Columns as desired
#colnames(merged_seurat@meta.data) [1] <- "orig.ident"

### Add number of genes per UMI for each cell to metadata
merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)

### Compute percent mitochondrial genes per cell
merged_seurat <- PercentageFeatureSet(merged_seurat, pattern = "^MT-", col.name = "percent_mito")

### Compute percentage of ribosomal genes per cell
merged_seurat <- PercentageFeatureSet(merged_seurat, pattern = "^RP[SL]", col.name = "percent_ribo")

# Percentage hemoglobin genes - includes all genes starting with HB except HBP.
merged_seurat <- PercentageFeatureSet(merged_seurat, "^HB[^(P)]", col.name = "percent_hb")

#Same for platelets
merged_seurat <- PercentageFeatureSet(merged_seurat, "PECAM1|PF4", col.name = "percent_plat")


## Create a metadata Object from seurat to work on

### Create metadata dataframe
metadata <- merged_seurat@meta.data

### Add cell IDs to metadata
metadata$cells <- rownames(metadata)

### Create Stim Column
metadata$PID <- NA
metadata$PID[which(str_detect(metadata$cells, "P1042"))] <- "P1042"
metadata$PID[which(str_detect(metadata$cells, "P1043"))] <- "P1043"
metadata$PID[which(str_detect(metadata$cells, "JLK"))] <- "JLK"
metadata

### Create Stim Column
metadata$stim <- NA
metadata$stim[which(str_detect(metadata$cells, "_Med"))] <- "Med"
metadata$stim[which(str_detect(metadata$cells, "_CD3"))] <- "CD3"
metadata$stim[which(str_detect(metadata$cells, "\\+_HIV"))] <- "HIV_peptide"
metadata

### Create Condition Column
metadata$condition <- NA
metadata$condition[which(str_detect(metadata$cells, "HIV-_"))] <- "HIV-"
metadata$condition[which(str_detect(metadata$cells, "HIV\\+_"))] <- "HIV+"
metadata

### After verifying created metadata is correct, add metadata back to seurat object
merged_seurat@meta.data <- metadata
merged_seurat@meta.data

# Create .RData object to load at any time
save(merged_seurat, file=paste0(out.path,"Seuratv5_CITEseq_dsbnorm_merged_Seurat.RData"))

########################################################################################################################


