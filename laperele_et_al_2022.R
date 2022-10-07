#!/bin/env Rscript

#### License Notice ####

##
# Copyright (c) 2022 Cedars-Sinai Medical Center
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
##

#### Title and Authors ####

##
# Title: iPSC-derived neural progenitor cells secreting GDNF slow disease 
#        progression in rodent models of both ALS and Retinal degeneration  
# Authors:  Alexander H. Laperle*1, Alexandra Moser*1, Veronica J. Garcia 1, 
#           Amanda Wu 1, Aaron Fulton 1, George Lawless 1, Shaughn Bell 1, 
#           Kristina Roxas 1, Roksana Elder 1, Pablo Avalos 1, Bin Lu 1, 
#           Staphany Ramirez 1, Shaomei Wang 1, Clive N. Svendsen**1 
#
# Affiliations: 1 Cedars-Sinai Board of Governors Regenerative Medicine 
#                 Institute, Cedars-Sinai Medical Center, Los Angeles CA
# *These authors contributed equally to this work 
# **Corresponding author email:  Clive.Svendsen@cshs.org 
##

#### Script Information ####

##
# R version 4.2.0
# R Script Title:  2022_laperle_et_al.R
# R Script Author:  Shaughn Bell
# R Script Corresponding Email:  shaughn.bell@cshs.org
##

#### Sample Information #### 

##
# GEO Accession Number:  GSE214210
#
# aligned to homo sapiens hg38 via CellRanger v6.1.2
# used the 10x Genomics supplied reference file "refdata-gex-GRCh38-2020-A.tar.gz"
#
# fastq files used:
#
# 1671_S15_L001_I1_001.fastq.gz
# 1671_S15_L001_R1_001.fastq.gz
# 1671_S15_L001_R2_001.fastq.gz
# 1671_S15_L002_I1_001.fastq.gz
# 1671_S15_L002_R1_001.fastq.gz
# 1671_S15_L002_R2_001.fastq.gz
# 3812_S13_L001_I1_001.fastq.gz
# 3812_S13_L001_R1_001.fastq.gz
# 3812_S13_L001_R2_001.fastq.gz
# 3812_S13_L002_I1_001.fastq.gz
# 3812_S13_L002_R1_001.fastq.gz
# 3812_S13_L002_R2_001.fastq.gz
# 4544_S9_L001_I1_001.fastq.gz
# 4544_S9_L001_R1_001.fastq.gz
# 4544_S9_L001_R2_001.fastq.gz
# 4544_S9_L002_I1_001.fastq.gz
# 4544_S9_L002_R1_001.fastq.gz
# 4544_S9_L002_R2_001.fastq.gz
#
# Cellranger outs used:  
#
# hg38_1671_outs/filtered_feature_bc_matrix/barcodes.tsv.gz
# hg38_1671_outs/filtered_feature_bc_matrix/features.tsv.gz
# hg38_1671_outs/filtered_feature_bc_matrix/matrix.mtx.gz
#
# hg38_3812_outs/filtered_feature_bc_matrix/barcodes.tsv.gz
# hg38_3812_outs/filtered_feature_bc_matrix/features.tsv.gz
# hg38_3812_outs/filtered_feature_bc_matrix/matrix.mtx.gz
#
# hg38_4544_outs/filtered_feature_bc_matrix/barcodes.tsv.gz
# hg38_4544_outs/filtered_feature_bc_matrix/features.tsv.gz
# hg38_4544_outs/filtered_feature_bc_matrix/matrix.mtx.gz
#
# Cellranger outs are not included in the GEO record.
#
# Genesummed matricies of counts:  
# 
# 1_4544_CNS10_GDNF_filtered_genesummed_mtx.csv
# 2_3812_iNPC_GDNF_WT_filtered_genesummed_mtx.csv
# 3_1671_iNPC_GDNF_T12_filtered_genesummed_mtx.csv
#
# Genesummed matrix files are included in the GEO record. 
#
##

#### project information ####

date <- "20220513"
project <- "iNPC_gdnf"
datadir <- "/20220513_iNPC_gdnf" # directory to save files generated
sourcedir <- "/20220419_inpc_outs_copy" # save Cellranger outs here

#### Load required packages ####

library(Seurat)
library(tidyverse)
library(Matrix)
library(ggplot2)
library(cowplot)
library(patchwork)
library(irlba)
library(gridExtra)
library(Vennerable)

#### Set working dir and save session info ####
setwd(datadir)

sessionInfo()

sink(paste0(date,"_",project,"_devtools_sessionInfo.txt"))
devtools::session_info()
sink()

sink(paste0(date,"_",project,"_sessionInfo.txt"))
sessionInfo()
sink()

#### Functions for saving images ####
hirestiff <- function(saveas){
  ggsave(
    saveas,
    plot = last_plot(),
    device = "tiff",
    scale = 1,
    dpi = 600,
    limitsize = TRUE,
    bg = "white"
  )
}

lowrestiff <- function(saveas){
  ggsave(
    saveas,
    plot = last_plot(),
    device = "tiff",
    scale = 1,
    dpi = 100,
    limitsize = TRUE,
    bg = "white"
  )
}

#### Gene Sum ####

  # Some gene symbols refer to more than one ENSG ID.  
  # In order to ensure all ENSG IDs are represented while using the much more
  # user-friendly gene symbol, the cellranger outpus were processed as follows:
  #    1) Cellranger output was loaded into a matrix object
  #    2) All genes with more than one ENSG ID were put into a new sub-matrix
  #    3) Counts from all instances of a given gene symbol were summed
  #    4) The genesummed sub-matrix was rejoined to the main matrix
  #    5) The ENSG ID column was removed
  # Genesummed matrix files are included in the GEO record.  
  # Cellranger outs are not included in the GEO record.

# Create matrix files from the cellranger outs
setwd(sourcedir)

filelist <- list.files(path = sourcedir)
filelist <- as.data.frame(filelist, files = list.files(path = sourcedir))
path_list <- paste0(sourcedir,"/",filelist[,1],"/filtered_feature_bc_matrix")
### IMPORTANT ###
  # Files with numbers in their name may be sorted in a 
  # different order depending on the way files are listed.  
  # Some systems sort by number instead of numeric value.
  # Make sure to use the same sort order as in the path_list object.
path_list

all(file.exists(path_list)) # check that all the files exist

#Read the cellranger output files and save expression matrix csv
for (k in 1:length(path_list)){
  
  # read in and transpose the barcode file
  bc <- as.data.frame(read.delim(file = paste0(path_list[k],"/barcodes.tsv.gz"), header = F))%>%t(.)
  # read in the features file and select just the first two columns
  feat <- as.data.frame(read.delim(file = paste0(path_list[k],"/features.tsv.gz"), header = F))%>%.[,1:2]
  # read in the counts
  matx <- as.data.frame(readMM(file = paste0(path_list[k],"/matrix.mtx.gz")))
  
  colnames(matx) <- bc[1,] #assign barcodes as colnames for the matrix
  colnames(feat) <- c("ensmbl.id", "symbol") #rename the columns of features table
  matx <- cbind(feat, matx) #cbind the features table and the matrix
  
  write.csv(matx, file = paste0(path_list[k],"/filtered_combined_mtx.csv"), row.names = FALSE)
  
}

rm(matx,bc,feat,filelist,k)

# Load the data

data.1671 <- read.csv(file=paste0(path_list[5],"/filtered_combined_mtx.csv"), header = TRUE)
data.3812 <- read.csv(file=paste0(path_list[6],"/filtered_combined_mtx.csv"), header = TRUE)
data.4544 <- read.csv(file=paste0(path_list[7],"/filtered_combined_mtx.csv"), header = TRUE)

# put all the loaded matrices into one list
obj.list <- list(data.1671,data.3812,data.4544)

# iterate the summing over all samples in obj.list
for (j in 1:length(obj.list)){
  
  #create a list of duplicated genes
  gene.duplicates<-obj.list[[j]]$symbol[duplicated(obj.list[[j]]$symbol, incomparables = NA)]%>%as.list(.)%>%unique(.)
  
  #Create empty dataframe to rbind duplicated rows to
  dup.rows<-data.frame()
  
  #Create identical expression matrix to subtract duplicated rows from
  sub.rows<-obj.list[[j]]
  
  #Extract rows containing the duplicated genes. 
  for (i in 1:length(gene.duplicates)){
    dups<-subset(obj.list[[j]], subset = symbol==gene.duplicates[[i]])#subsets duplicated rows
    dup.rows<-rbind(dup.rows, dups)#binds duplicated rows together
    sub.rows<-subset(sub.rows, subset = symbol!=gene.duplicates[[i]])#deleted duplicated rows from expression matrix
  }
  
  #Remove ensembl ID columns
  dup.rows<-dup.rows[,2:length(colnames(dup.rows))]
  sub.rows<-sub.rows[,2:length(colnames(sub.rows))]
  
  #Find the gsum of duplicated rows
  dup.rows<-aggregate(. ~ symbol, data = dup.rows, sum)
  
  #rbind the gsum of duplicated rows with the expression matrix that had these genes removed
  new.mat<-rbind(sub.rows, dup.rows)
  
  if (nrow(new.mat) == length(unique(obj.list[[j]]$symbol))) {
    new.mat <- as.data.frame(new.mat)
    rownames(new.mat) <- new.mat[,1]
    new.mat <- subset(new.mat, select = -c(symbol))
    obj.list[[j]] <- new.mat
    message("Successfully merged!")
  } else {
    message ("There's an error!")
  }
  
}

data.1671 <- obj.list[[1]]
data.3812 <- obj.list[[2]]
data.4544 <- obj.list[[3]]

# save genesummed matrix

write.csv(data.1671, file = paste0(path_list[5],"/filtered_genesummed_mtx.csv"), row.names = FALSE)
write.csv(data.3812, file = paste0(path_list[6],"/filtered_genesummed_mtx.csv"), row.names = FALSE)
write.csv(data.4544, file = paste0(path_list[7],"/filtered_genesummed_mtx.csv"), row.names = FALSE)

# save R objects for easy reloading
save(data.1671,data.3812,data.4544,
     file = paste0(datadir,"/",date,"_",project,"_genesummed.rdata"))

# load(file = paste0(datadir,"/",date,"_",project,"_genesummed.rdata"))

#### Seurat Object ####

# set wd to datadir
setwd(datadir)

# Initialize the Seurat object with the genesummed data.
# Keep all genes expressed in >= 1 cell

s1671 <- CreateSeuratObject(counts = obj.list[[1]], project = "1671", min.cells = 1, min.features = 0)
s3812 <- CreateSeuratObject(counts = obj.list[[2]], project = "3812", min.cells = 1, min.features = 0)
s4544 <- CreateSeuratObject(counts = obj.list[[3]], project = "4544", min.cells = 1, min.features = 0)

rm(data.1671,data.3812,data.4544,
   dup.rows,dups,gene.duplicates,new.mat,
   obj.list,sub.rows,i,j,path_list)

save(s1671,s3812,s4544,
     file = paste0(datadir,"/",date,"_",project,"_seuratobjs.rdata"))

# load(file = paste0(datadir,"/",date,"_",project,"_seuratobjs.rdata"))

#### Add sample metadata ####

s4544[["Vial"]] <-"1024544"
s3812[["Vial"]] <-"1023812"
s1671[["Vial"]] <-"991671"

s4544[["Lot"]] <-"CNS10-GDNF"
s3812[["Lot"]] <-"iNPC-GDNF-WT"
s1671[["Lot"]] <-"iNPC-GDNF-T12"

s4544[["Passage"]] <-"p27"
s3812[["Passage"]] <-"p23"
s1671[["Passage"]] <-"p17"

#### Merge all data into a single seurat object ####

tenx1 <- merge(s4544, y = c(s3812,s1671),
               add.cell.ids = c("4544","3812","1671"), 
               project = project)

tenx1

rm(s1671,s3812,s4544)

##save merged seurat object 
saveRDS(tenx1, file = paste0("./",date,"_",project,"_merged_seurat_prefilter.rds"))
# beep()

# tenx1 <- read_rds(file = paste0("./",date,"_",project,"_merged_seurat_prefilter.rds"))

#### Initial QC Filtering ####

data <- tenx1

# Add mito info to metadata
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
data[["percent.ribo"]] <- PercentageFeatureSet(data, pattern = "^RP[SL]")

##Plots of UMI#, Gene#, and %mito
vp1 <- VlnPlot(data, features = "nCount_RNA", pt.size = 0)+
  ggtitle("nCount_RNA prefiltered")
vp1
pdf(paste0("./",date,"_",project,"_prefilter_","nCount_RNA",".pdf"))
print(vp1)
dev.off()

vp2 <- VlnPlot(data, features = "nFeature_RNA", pt.size = 0) +
  ggtitle("nFeature_RNA prefiltered")
vp2
pdf(paste0("./",date,"_",project,"_prefilter_","nFeature_RNA",".pdf"))
print(vp2)
dev.off()

vp3 <- VlnPlot(data, features = "percent.mt", pt.size = 0) +
  ggtitle("percent.mt prefiltered")
vp3
pdf(paste0("./",date,"_",project,"_prefilter_","percent.mt",".pdf"))
print(vp3)
dev.off()

vp4 <- VlnPlot(data, features = "percent.ribo", pt.size = 0) +
  ggtitle("percent.ribo prefiltered")
vp4
pdf(paste0("./",date,"_",project,"_prefilter_","percent.ribo",".pdf"))
print(vp4)
dev.off()

p1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt") 
p2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 
p1 + p2

pdf(paste0("./",date,"_",project,"_prefilter_","scatterplots",".pdf"))
print(p1 + p2)
dev.off()

# Additional QC and Z-score QC metrics
data[["nUMI.z"]] <- scale(data$nCount_RNA)
data[["nGene.z"]] <- scale(data$nFeature_RNA)
data[["percent.mt.z"]] <- scale(data$percent.mt)
data[["percent.ribo.z"]] <- scale(data$percent.ribo)

##Filter cells based on Z-score
length(data@meta.data$orig.ident)
mean(data@meta.data$percent.mt)
median(data@meta.data$nCount_RNA)
median(data@meta.data$nFeature_RNA)
max(data@meta.data$nCount_RNA)

sink(paste0("./",date,"_",project,"_prefilter_features_counts_pctmito.txt"))
cat("length")
length(data@meta.data$orig.ident)
cat("mean percent mito")
mean(data@meta.data$percent.mt)
cat("median counts")
median(data@meta.data$nCount_RNA)
cat("median features")
median(data@meta.data$nFeature_RNA)
cat("max counts")
max(data@meta.data$nCount_RNA)
sink()

data <- subset(data, subset = percent.mt.z < 3)
data <- subset(data, subset = nGene.z < 3)
data <- subset(data, subset = nUMI.z < 3)
data <- subset(data, subset = percent.ribo.z < 3)

length(data@meta.data$orig.ident)
mean(data@meta.data$percent.mt)
median(data@meta.data$nCount_RNA)
median(data@meta.data$nFeature_RNA)
max(data@meta.data$nCount_RNA)

sink(paste0("./",date,"_",project,"_postfilter_features_counts_pctmito.txt"))
cat("length")
length(data@meta.data$orig.ident)
cat("mean percent mito")
mean(data@meta.data$percent.mt)
cat("median counts")
median(data@meta.data$nCount_RNA)
cat("median features")
median(data@meta.data$nFeature_RNA)
cat("max counts")
max(data@meta.data$nCount_RNA)
sink()

##Plots of UMI#, Gene#, and %mito
vp5 <- VlnPlot(data, features = "nCount_RNA", pt.size = 0, group.by = "orig.ident")  + 
  NoLegend() + ggtitle("nCount_RNA filtered")
vp5
pdf(paste0("./",date,"_",project,"_filtered_","nCount_RNA",".pdf"))
print(vp5)
dev.off()

vp6 <- VlnPlot(data, features = "nFeature_RNA", pt.size = 0, group.by = "orig.ident") + 
  NoLegend() + ggtitle("nFeature_RNA filtered")
vp6
pdf(paste0("./",date,"_",project,"_filtered_","nFeature_RNA",".pdf"))
print(vp6)
dev.off()

vp7 <- VlnPlot(data, features = "percent.mt", pt.size = 0, group.by = "orig.ident") +
  NoLegend() + ggtitle("percent.mt filtered")
vp7
pdf(paste0("./",date,"_",project,"_filtered_","percent.mt",".pdf"))
print(vp7)
dev.off()

vp8 <- VlnPlot(data, features = "percent.ribo", pt.size = 0, group.by = "orig.ident") +
  NoLegend() + ggtitle("percent.ribo filtered")
vp8
pdf(paste0("./",date,"_",project,"_filtered_","percent.ribo",".pdf"))
print(vp8)
dev.off()

p3 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt") 
p4 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 
p3 + p4

pdf(paste0("./",date,"_",project,"_filtered_","scatterplots",".pdf"))
print(p3 + p4)
dev.off()

# Additional QC metrics
qc.metrics <- as_tibble(data[[c("nCount_RNA",
                                "nFeature_RNA",
                                "percent.mt",
                                "percent.ribo")]],
                        rownames="Cell.Barcode")

p5 <- qc.metrics %>%
  arrange(percent.mt) %>%
  ggplot(aes(nCount_RNA,nFeature_RNA,colour=percent.mt)) + 
  geom_point() + 
  scale_color_gradientn(colors=c("black","blue","green2","red","yellow")) +
  ggtitle("QC metrics")
p5
pdf(paste0("./",date,"_",project,"QC_Metrics",".pdf"))
print(p5)
dev.off()

p6 <- qc.metrics %>%
  ggplot(aes(percent.mt)) + 
  geom_histogram(binwidth = 0.5, fill="yellow", colour="black") +
  ggtitle("Distribution of Percentage Mitochondrion")
p6
pdf(paste0("./",date,"_",project,"Dist_mito",".pdf"))
print(p6)
dev.off()

p7 <- qc.metrics %>%
  ggplot(aes(percent.ribo)) + 
  geom_histogram(binwidth = 0.5, fill="yellow", colour="black") +
  ggtitle("Distribution of Percentage Ribosomal")
p7
pdf(paste0("./",date,"_",project,"Dist_ribo",".pdf"))
print(p7)
dev.off()

rm(vp1,vp2,vp3,vp4,vp5,vp6,vp7,vp8,
   p1,p2,p3,p4,p5,p6,p7,qc.metrics)

tenx1 <- data

rm(data)

#### perform integration ####
s.list <- SplitObject(tenx1, split.by = "orig.ident")

# normalize data and find variable features
s.list <- lapply(X = s.list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE)
})

# Next, select features for downstream integration, and run PCA on each 
# object in the list, which is required for running the alternative 
# reciprocal PCA workflow.

features <- SelectIntegrationFeatures(object.list = s.list)
s.list <- lapply(X = s.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

anchors <- FindIntegrationAnchors(object.list = s.list, 
                                  reduction = "rpca",
                                  dims = 1:50)
s.integrated <- IntegrateData(anchorset = anchors, dims = 1:50)
s.integrated <- ScaleData(s.integrated, verbose = FALSE)

#### set factor levels ####

s.integrated$orig.ident <- factor(x = s.integrated$orig.ident, levels = c("4544","3812","1671"))
s.integrated$Lot <- factor(x = s.integrated$Lot, levels = c("CNS10-GDNF","iNPC-GDNF-WT","iNPC-GDNF-T12"))

saveRDS(s.integrated, file = paste0("./",date,"_",project,"_post_integration.rds"))

# s.integrated <- read_rds(file = paste0("./",date,"_",project,"_post_integration.rds"))

#### PCA ####

# this performs PCA on the seurat object
s.integrated <- RunPCA(s.integrated, npcs = 50, verbose = TRUE)

# make PC coordinate object a data frame
xx.coord <- as.data.frame(s.integrated@reductions$pca@cell.embeddings)

# make PC feature loadings object a data frame
xx.gload <- as.data.frame(s.integrated@reductions$pca@feature.loadings)

# calculate eigenvalues for arrays

# generate squares of all sample coordinates
sq.xx.coord <- as.data.frame(xx.coord^2)
# create empty list for eigenvalues first
eig <- c()
# calculate the eigenvalue for each PC in sq.xx.coord by taking the sqrt of the sum of squares
for(i in 1:ncol(sq.xx.coord))
  eig[i] = sqrt(sum(sq.xx.coord[,i]))
# calculate the total variance by adding up all the eigenvalues
sum.eig <- sum(eig)
# calculate the expected contribution of all PCs if they all contribute equally to the total variance
expected.contribution <- sum.eig/(length(xx.coord)-1)
# return the number of principal components with an eigenvalue greater than expected by equal variance
meaningful.PCs <- sum(eig > expected.contribution)

# create empty list for eigenvalue percentage
eig.percent <- c()
# calculate the percentage of the total variance by each PC eigenvalue
for(i in 1:length(eig))
  eig.percent[i] = 100*eig[i]/sum.eig
# sum of all eig.percent should total to 100
sum(eig.percent)
# create empty list for scree values
scree <- c()
# calculate a running total of variance contribution
for(i in 1:length(eig))
  if(i == 1) scree[i] = eig.percent[i] else scree[i] = scree[i-1] + eig.percent[i]

# create data frame for eigenvalue summaries
eigenvalues <- data.frame("PC" = colnames(xx.coord), "eig" = eig, "percent" = eig.percent, "scree" = scree)

# plot scree values
plot(eigenvalues$percent, ylim = c(0,100), type = "S", xlab = "PC", ylab = "Percent of variance",
     main = paste0(date,"_",project," scree plot all samples PCA"))
points(eigenvalues$scree, ylim = c(0,100), type = "p", pch = 16)
lines(eigenvalues$scree)
# add red line to indicate cut-off
cut.off <- 100/(length(eig)-1)
abline(h = cut.off, col = "red")
# add blue line to indicate which PCs are meaningful and kept
abline(v = meaningful.PCs, col = "blue")
text(meaningful.PCs, cut.off, label = paste("cutoff PC",meaningful.PCs),
     adj = c(-0.1, -0.5))

dev.copy(pdf, paste0("./",date,"_",project,"_scree_plot.pdf"))
dev.off()

# meaningful.PCs <- 14

#### Run UMAP and look at UMAP plots ####

s.integrated <- RunUMAP(s.integrated, reduction = "pca", dims = 1:meaningful.PCs, verbose = TRUE)

# update prefixed variable
prefixPC <- paste0("./",date,"_",project,"_",meaningful.PCs,"PCs")

## UMAP plot by sample name ("orig.ident")
DimPlot(s.integrated, reduction = "umap", label = FALSE, pt.size = .25, group.by = "orig.ident")
hirestiff(paste0(prefixPC,"_UMAP_by_sample_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_by_sample_lowres.tiff"))
DimPlot(s.integrated, reduction = "umap", label = FALSE, pt.size = .25, group.by = "Lot")
hirestiff(paste0(prefixPC,"_UMAP_by_lot_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_by_lot_lowres.tiff"))
DimPlot(s.integrated, reduction = "umap", label = FALSE, pt.size = .25, split.by = "Lot") + NoLegend()
hirestiff(paste0(prefixPC,"_UMAP_by_lot_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_by_lot_lowres.tiff"))

#### Cell Cycle Score ####

DefaultAssay(s.integrated) <- "RNA"

s.integrated <- CellCycleScoring(s.integrated,
                               s.features = cc.genes.updated.2019$s.genes,
                               g2m.features = cc.genes.updated.2019$g2m.genes,
                               set.ident = TRUE)

as_tibble(s.integrated[[]]) %>%
  ggplot(aes(Phase)) + geom_bar()

DimPlot(s.integrated, reduction = "umap", label = FALSE, pt.size = .25, group.by = "Phase")
hirestiff(paste0(prefixPC,"_UMAP_by_phase_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_by_phase_lowres.tiff"))

DimPlot(s.integrated, reduction = "umap", label = FALSE, pt.size = .25, split.by = "Lot", group.by = "Phase")
hirestiff(paste0(prefixPC,"_UMAP_by_phase_splitby_orig.ident_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_by_phase_splitby_orig.ident_lowres.tiff"))

DefaultAssay(s.integrated) <- "integrated"

saveRDS(s.integrated, file = paste0(prefixPC,"_seurat_integrated.rds"))

# s.integrated <- read_rds(file = paste0(prefixPC,"_seurat_integrated.rds"))

#### Clustering and Resolution ####

# Determine the K-nearest neighbor graph
s.integrated <- FindNeighbors(object = s.integrated, reduction = "pca", dims = 1:meaningful.PCs)
# s.integrated <- FindNeighbors(object = s.integrated, reduction = "pca", dims = 1:8)

# Determine the clusters                              
s.integrated <- FindClusters(object = s.integrated,
                           resolution = c(0.2))

#update prefix
prefixPCres <- paste0(prefixPC,"_res0.2")

# Plot the UMAP
DimPlot(s.integrated, reduction = "umap", label = TRUE, label.size = 5)  + NoLegend() + 
  ggtitle("integrated_snn_res.0.2")
hirestiff(paste0(prefixPCres,"_UMAP","_by_","cluster","_hires.tiff"))
lowrestiff(paste0(prefixPCres,"_UMAP","_by_","cluster","_lowres.tiff"))

# UMAP of cells in each cluster by orig.ident
DimPlot(s.integrated, reduction = "umap", label = FALSE, split.by = "Lot")  + NoLegend()
hirestiff(paste0(prefixPCres,"_UMAP","_by_","orig.ident","_hires.tiff"))
lowrestiff(paste0(prefixPCres,"_UMAP","_by_","orig.ident","_lowres.tiff"))

# UMAP of cells in each cluster by orig.ident with cluster labels
DimPlot(s.integrated, reduction = "umap", label = TRUE, split.by = "Lot")  + NoLegend()
hirestiff(paste0(prefixPCres,"_UMAP","_by_","orig.ident_with_clusters","_hires.tiff"))
lowrestiff(paste0(prefixPCres,"_UMAP","_by_","orig.ident_with_clusters","_lowres.tiff"))

for (i in 1:length(levels(s.integrated$Type))){
  p <- (DimPlot(s.integrated, 
                reduction = "umap", 
                label = FALSE, 
                pt.size = .25, 
                cells = c(WhichCells(s.integrated, expression = Type == levels(s.integrated$Type)[i]))) + 
          NoLegend() +
          ggtitle(paste0(levels(s.integrated$Type)[i])))
  
  print(p)
  
  hirestiff(paste(prefixPCres,"UMAP","split_by","type","no_cluster_labels",
                  levels(s.integrated$Type)[i],"hires.tiff",sep = "_"))
  lowrestiff(paste(prefixPCres,"UMAP","split_by","type","no_cluster_labels",
                   levels(s.integrated$Type)[i],"lowres.tiff",sep = "_"))
}

# Extract number of cells per cluster per Lot
n_cells <- FetchData(s.integrated, vars = c("ident", "Lot")) %>%
  dplyr::count(ident, Lot) %>%
  tidyr::spread(ident, n)

write.csv(n_cells, file = paste0(prefixPCres,"_cells_per_cluster.csv"))

# save RDS containing reduction and cluster idents
saveRDS(s.integrated, paste0(prefixPCres,"_seurat_after_clustering.rds"))

# s.integrated <- readRDS(file = paste0(prefixPCres,"_seurat_after_clustering.rds"))

#### Normalize RNA slot for visualization ####

# Select the RNA counts slot to be the default assay for visualization purposes
DefaultAssay(s.integrated) <- "RNA"

# Normalize, find variable features, scale data 
s.integrated <- NormalizeData(s.integrated, verbose = FALSE)
s.integrated <- FindVariableFeatures(s.integrated)
all.genes <- rownames(s.integrated)
s.integrated <- ScaleData(s.integrated, features = all.genes)

# save RDS containing RNA normalized data
saveRDS(s.integrated, paste0(prefixPCres,"_seurat_after_RNAnorm.rds"))

# s.integrated <- readRDS(file = paste0(prefixPCres,"_seurat_after_RNAnorm.rds"))

#### Figure 1b ####

# Plot the UMAP
DimPlot(s.integrated, reduction = "umap", label = TRUE, label.size = 5)  + 
  NoLegend() + 
  
hirestiff(paste0(prefixPCres,"_Fig","_1b_","_hires.tiff"))
lowrestiff(paste0(prefixPCres,"_Fig","_1b_","_lowres.tiff"))

#### Figure 1c ####
# for loop that saves UMAP of each Type by cluster without cluster labels

for (i in 1:length(levels(s.integrated$Type))){
  p <- (DimPlot(s.integrated, 
                reduction = "umap", 
                label = FALSE, 
                pt.size = .25, 
                cells = c(WhichCells(s.integrated, expression = Type == levels(s.integrated$Type)[i]))) + 
          NoLegend() +
          ggtitle(paste0(levels(s.integrated$Type)[i])))
  
  print(p)
  
  hirestiff(paste(prefixPCres,"Fig","1c",
                  levels(s.integrated$Type)[i],"hires.tiff",sep = "_"))
  lowrestiff(paste(prefixPCres,"Fig","1c",
                   levels(s.integrated$Type)[i],"lowres.tiff",sep = "_"))
}


#### Figure 1d ####

FeaturePlot(s.integrated, 
            reduction = "umap", 
            features = c("S100B","VIM","CD99","NES"),
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE,
            label.size = 3,
            pt.size = 0.5)

hirestiff(paste(prefixPCres,"Fig","1d","hires.tiff",sep = "_"))
lowrestiff(paste(prefixPCres,"Fig","1d","lowres.tiff",sep = "_"))

#### Figure 1e ####

FeaturePlot(s.integrated, 
            reduction = "umap", 
            features = c("GFAP","AQP4","SLC1A3"),
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE,
            label.size = 3,
            pt.size = 0.5)

hirestiff(paste(prefixPCres,"Fig","1e","hires.tiff",sep = "_"))
lowrestiff(paste(prefixPCres,"Fig","1e","lowres.tiff",sep = "_"))

#### Figure S1a ####

DimPlot(s.integrated, 
        reduction = "umap", 
        label = FALSE, 
        pt.size = .25, 
        group.by = "Phase")

hirestiff(paste0(prefixPCres,"_Fig_S1a_hires.tiff"))
lowrestiff(paste0(prefixPCres,"_Fig_S1a_lowres.tiff"))

#### Figure S1b ####

FeaturePlot(s.integrated, 
            reduction = "umap", 
            features = c("POU5F1","NANOG","SOX2","KLF4"),
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE,
            label.size = 3,
            pt.size = 0.5)

hirestiff(paste(prefixPCres,"Fig","S1b","hires.tiff",sep = "_"))
lowrestiff(paste(prefixPCres,"Fig","S1b","lowres.tiff",sep = "_"))

#### Figure S1c ####

expressiondata <- FetchData(object = s.integrated, 
                            vars = c("orig.ident","Lot",
                                     "POU5F1","NANOG","KLF4"))

expressiondata <- rownames_to_column(expressiondata, var = "cell_id")

cells_pouf1 <- expressiondata$cell_id[expressiondata$POU5F1 > 0]
cells_nanog <- expressiondata$cell_id[expressiondata$NANOG > 0]
cells_klf4 <- expressiondata$cell_id[expressiondata$KLF4 > 0]

genelists <- list(cells_pouf1,cells_nanog,cells_klf4)
names(genelists) <- c("OCT4","NANOG","KLF4")

venngenes <- Venn(genelists)

plot(venngenes, doWeights = FALSE, type = "circles")

pdf(paste(prefixPCres,"Fig","S1c.pdf",sep = "_"))
plot(venngenes, doWeights = FALSE, type = "circles")
dev.off()
#### write metadata ####

write.csv(s.integrated@meta.data, file = paste(prefixPCres,"metadata.csv",sep = "_"))



#### session info ####
 
# R version 4.2.0 (2022-04-22)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Monterey 12.4
# 
# Matrix products: default
# LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] beepr_1.3          irlba_2.3.5        pheatmap_1.0.12    patchwork_1.1.1    viridis_0.6.2      viridisLite_0.4.0 
# [7] cowplot_1.1.1      Matrix_1.4-1       forcats_0.5.1      stringr_1.4.0      dplyr_1.0.9        purrr_0.3.4       
# [13] readr_2.1.2        tidyr_1.2.0        tibble_3.1.7       ggplot2_3.3.6      tidyverse_1.3.1    sp_1.4-7          
# [19] SeuratObject_4.1.0 Seurat_4.1.1      
# 
# loaded via a namespace (and not attached):
#   [1] readxl_1.4.0          backports_1.4.1       plyr_1.8.7            igraph_1.3.1          lazyeval_0.2.2       
# [6] splines_4.2.0         listenv_0.8.0         scattermore_0.8       usethis_2.1.5         digest_0.6.29        
# [11] htmltools_0.5.2       fansi_1.0.3           magrittr_2.0.3        memoise_2.0.1         tensor_1.5           
# [16] cluster_2.1.3         ROCR_1.0-11           tzdb_0.3.0            NCmisc_1.1.6          remotes_2.4.2        
# [21] globals_0.15.0        modelr_0.1.8          matrixStats_0.62.0    spatstat.sparse_2.1-1 prettyunits_1.1.1    
# [26] colorspace_2.0-3      rvest_1.0.2           ggrepel_0.9.1         haven_2.5.0           callr_3.7.0          
# [31] crayon_1.5.1          jsonlite_1.8.0        progressr_0.10.0      spatstat.data_2.2-0   survival_3.3-1       
# [36] zoo_1.8-10            glue_1.6.2            polyclip_1.10-0       gtable_0.3.0          leiden_0.4.2         
# [41] pkgbuild_1.3.1        future.apply_1.9.0    abind_1.4-5           scales_1.2.0          DBI_1.1.2            
# [46] spatstat.random_2.2-0 miniUI_0.1.1.1        Rcpp_1.0.8.3          xtable_1.8-4          reticulate_1.24      
# [51] spatstat.core_2.4-2   htmlwidgets_1.5.4     httr_1.4.3            RColorBrewer_1.1-3    ellipsis_0.3.2       
# [56] ica_1.0-2             pkgconfig_2.0.3       uwot_0.1.11           dbplyr_2.1.1          deldir_1.0-6         
# [61] utf8_1.2.2            tidyselect_1.1.2      rlang_1.0.2           reshape2_1.4.4        later_1.3.0          
# [66] cellranger_1.1.0      munsell_0.5.0         tools_4.2.0           cachem_1.0.6          cli_3.3.0            
# [71] audio_0.1-10          generics_0.1.2        devtools_2.4.3        broom_0.8.0           ggridges_0.5.3       
# [76] fastmap_1.1.0         goftest_1.2-3         processx_3.5.3        fs_1.5.2              fitdistrplus_1.1-8   
# [81] RANN_2.6.1            pbapply_1.5-0         future_1.25.0         nlme_3.1-157          mime_0.12            
# [86] proftools_0.99-3      xml2_1.3.3            rstudioapi_0.13       brio_1.1.3            compiler_4.2.0       
# [91] plotly_4.10.0         png_0.1-7             testthat_3.1.4        spatstat.utils_2.3-1  reprex_2.0.1         
# [96] stringi_1.7.6         ps_1.7.0              desc_1.4.1            rgeos_0.5-9           lattice_0.20-45      
# [101] vctrs_0.4.1           pillar_1.7.0          lifecycle_1.0.1       BiocManager_1.30.17   spatstat.geom_2.4-0  
# [106] lmtest_0.9-40         RcppAnnoy_0.0.19      data.table_1.14.2     httpuv_1.6.5          R6_2.5.1             
# [111] promises_1.2.0.1      KernSmooth_2.23-20    gridExtra_2.3         parallelly_1.31.1     sessioninfo_1.2.2    
# [116] codetools_0.2-18      MASS_7.3-57           assertthat_0.2.1      pkgload_1.2.4         rprojroot_2.0.3      
# [121] withr_2.5.0           sctransform_0.3.3     mgcv_1.8-40           parallel_4.2.0        hms_1.1.1            
# [126] grid_4.2.0            rpart_4.1.16          Rtsne_0.16            shiny_1.7.1           lubridate_1.8.0      
# 

#### devtools session info ####

# ─ Session info ─────────────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.2.0 (2022-04-22)
# os       macOS Monterey 12.4
# system   x86_64, darwin17.0
# ui       RStudio
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       America/Los_Angeles
# date     2022-05-16
# rstudio  2022.02.2+485 Prairie Trillium (desktop)
# pandoc   NA
# 
# ─ Packages ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# package         * version date (UTC) lib source
# abind             1.4-5   2016-07-21 [1] CRAN (R 4.2.0)
# assertthat        0.2.1   2019-03-21 [1] CRAN (R 4.2.0)
# audio             0.1-10  2021-11-25 [1] CRAN (R 4.2.0)
# backports         1.4.1   2021-12-13 [1] CRAN (R 4.2.0)
# beepr           * 1.3     2018-06-04 [1] CRAN (R 4.2.0)
# BiocManager       1.30.17 2022-04-22 [1] CRAN (R 4.2.0)
# brio              1.1.3   2021-11-30 [1] CRAN (R 4.2.0)
# broom             0.8.0   2022-04-13 [1] CRAN (R 4.2.0)
# cachem            1.0.6   2021-08-19 [1] CRAN (R 4.2.0)
# callr             3.7.0   2021-04-20 [1] CRAN (R 4.2.0)
# cellranger        1.1.0   2016-07-27 [1] CRAN (R 4.2.0)
# cli               3.3.0   2022-04-25 [1] CRAN (R 4.2.0)
# cluster           2.1.3   2022-03-28 [1] CRAN (R 4.2.0)
# codetools         0.2-18  2020-11-04 [1] CRAN (R 4.2.0)
# colorspace        2.0-3   2022-02-21 [1] CRAN (R 4.2.0)
# cowplot         * 1.1.1   2020-12-30 [1] CRAN (R 4.2.0)
# crayon            1.5.1   2022-03-26 [1] CRAN (R 4.2.0)
# data.table        1.14.2  2021-09-27 [1] CRAN (R 4.2.0)
# DBI               1.1.2   2021-12-20 [1] CRAN (R 4.2.0)
# dbplyr            2.1.1   2021-04-06 [1] CRAN (R 4.2.0)
# deldir            1.0-6   2021-10-23 [1] CRAN (R 4.2.0)
# desc              1.4.1   2022-03-06 [1] CRAN (R 4.2.0)
# devtools          2.4.3   2021-11-30 [1] CRAN (R 4.2.0)
# digest            0.6.29  2021-12-01 [1] CRAN (R 4.2.0)
# dplyr           * 1.0.9   2022-04-28 [1] CRAN (R 4.2.0)
# ellipsis          0.3.2   2021-04-29 [1] CRAN (R 4.2.0)
# fansi             1.0.3   2022-03-24 [1] CRAN (R 4.2.0)
# fastmap           1.1.0   2021-01-25 [1] CRAN (R 4.2.0)
# fitdistrplus      1.1-8   2022-03-10 [1] CRAN (R 4.2.0)
# forcats         * 0.5.1   2021-01-27 [1] CRAN (R 4.2.0)
# fs                1.5.2   2021-12-08 [1] CRAN (R 4.2.0)
# future            1.25.0  2022-04-24 [1] CRAN (R 4.2.0)
# future.apply      1.9.0   2022-04-25 [1] CRAN (R 4.2.0)
# generics          0.1.2   2022-01-31 [1] CRAN (R 4.2.0)
# ggplot2         * 3.3.6   2022-05-03 [1] CRAN (R 4.2.0)
# ggrepel           0.9.1   2021-01-15 [1] CRAN (R 4.2.0)
# ggridges          0.5.3   2021-01-08 [1] CRAN (R 4.2.0)
# globals           0.15.0  2022-05-09 [1] CRAN (R 4.2.0)
# glue              1.6.2   2022-02-24 [1] CRAN (R 4.2.0)
# goftest           1.2-3   2021-10-07 [1] CRAN (R 4.2.0)
# gridExtra         2.3     2017-09-09 [1] CRAN (R 4.2.0)
# gtable            0.3.0   2019-03-25 [1] CRAN (R 4.2.0)
# haven             2.5.0   2022-04-15 [1] CRAN (R 4.2.0)
# hms               1.1.1   2021-09-26 [1] CRAN (R 4.2.0)
# htmltools         0.5.2   2021-08-25 [1] CRAN (R 4.2.0)
# htmlwidgets       1.5.4   2021-09-08 [1] CRAN (R 4.2.0)
# httpuv            1.6.5   2022-01-05 [1] CRAN (R 4.2.0)
# httr              1.4.3   2022-05-04 [1] CRAN (R 4.2.0)
# ica               1.0-2   2018-05-24 [1] CRAN (R 4.2.0)
# igraph            1.3.1   2022-04-20 [1] CRAN (R 4.2.0)
# irlba           * 2.3.5   2021-12-06 [1] CRAN (R 4.2.0)
# jsonlite          1.8.0   2022-02-22 [1] CRAN (R 4.2.0)
# KernSmooth        2.23-20 2021-05-03 [1] CRAN (R 4.2.0)
# later             1.3.0   2021-08-18 [1] CRAN (R 4.2.0)
# lattice           0.20-45 2021-09-22 [1] CRAN (R 4.2.0)
# lazyeval          0.2.2   2019-03-15 [1] CRAN (R 4.2.0)
# leiden            0.4.2   2022-05-09 [1] CRAN (R 4.2.0)
# lifecycle         1.0.1   2021-09-24 [1] CRAN (R 4.2.0)
# listenv           0.8.0   2019-12-05 [1] CRAN (R 4.2.0)
# lmtest            0.9-40  2022-03-21 [1] CRAN (R 4.2.0)
# lubridate         1.8.0   2021-10-07 [1] CRAN (R 4.2.0)
# magrittr          2.0.3   2022-03-30 [1] CRAN (R 4.2.0)
# MASS              7.3-57  2022-04-22 [1] CRAN (R 4.2.0)
# Matrix          * 1.4-1   2022-03-23 [1] CRAN (R 4.2.0)
# matrixStats       0.62.0  2022-04-19 [1] CRAN (R 4.2.0)
# memoise           2.0.1   2021-11-26 [1] CRAN (R 4.2.0)
# mgcv              1.8-40  2022-03-29 [1] CRAN (R 4.2.0)
# mime              0.12    2021-09-28 [1] CRAN (R 4.2.0)
# miniUI            0.1.1.1 2018-05-18 [1] CRAN (R 4.2.0)
# modelr            0.1.8   2020-05-19 [1] CRAN (R 4.2.0)
# munsell           0.5.0   2018-06-12 [1] CRAN (R 4.2.0)
# NCmisc            1.1.6   2018-11-12 [1] CRAN (R 4.2.0)
# nlme              3.1-157 2022-03-25 [1] CRAN (R 4.2.0)
# parallelly        1.31.1  2022-04-22 [1] CRAN (R 4.2.0)
# patchwork       * 1.1.1   2020-12-17 [1] CRAN (R 4.2.0)
# pbapply           1.5-0   2021-09-16 [1] CRAN (R 4.2.0)
# pheatmap        * 1.0.12  2019-01-04 [1] CRAN (R 4.2.0)
# pillar            1.7.0   2022-02-01 [1] CRAN (R 4.2.0)
# pkgbuild          1.3.1   2021-12-20 [1] CRAN (R 4.2.0)
# pkgconfig         2.0.3   2019-09-22 [1] CRAN (R 4.2.0)
# pkgload           1.2.4   2021-11-30 [1] CRAN (R 4.2.0)
# plotly            4.10.0  2021-10-09 [1] CRAN (R 4.2.0)
# plyr              1.8.7   2022-03-24 [1] CRAN (R 4.2.0)
# png               0.1-7   2013-12-03 [1] CRAN (R 4.2.0)
# polyclip          1.10-0  2019-03-14 [1] CRAN (R 4.2.0)
# prettyunits       1.1.1   2020-01-24 [1] CRAN (R 4.2.0)
# processx          3.5.3   2022-03-25 [1] CRAN (R 4.2.0)
# proftools         0.99-3  2020-07-08 [1] CRAN (R 4.2.0)
# progressr         0.10.0  2021-12-19 [1] CRAN (R 4.2.0)
# promises          1.2.0.1 2021-02-11 [1] CRAN (R 4.2.0)
# ps                1.7.0   2022-04-23 [1] CRAN (R 4.2.0)
# purrr           * 0.3.4   2020-04-17 [1] CRAN (R 4.2.0)
# R6                2.5.1   2021-08-19 [1] CRAN (R 4.2.0)
# RANN              2.6.1   2019-01-08 [1] CRAN (R 4.2.0)
# RColorBrewer      1.1-3   2022-04-03 [1] CRAN (R 4.2.0)
# Rcpp              1.0.8.3 2022-03-17 [1] CRAN (R 4.2.0)
# RcppAnnoy         0.0.19  2021-07-30 [1] CRAN (R 4.2.0)
# readr           * 2.1.2   2022-01-30 [1] CRAN (R 4.2.0)
# readxl            1.4.0   2022-03-28 [1] CRAN (R 4.2.0)
# remotes           2.4.2   2021-11-30 [1] CRAN (R 4.2.0)
# reprex            2.0.1   2021-08-05 [1] CRAN (R 4.2.0)
# reshape2          1.4.4   2020-04-09 [1] CRAN (R 4.2.0)
# reticulate        1.24    2022-01-26 [1] CRAN (R 4.2.0)
# rgeos             0.5-9   2021-12-15 [1] CRAN (R 4.2.0)
# rlang             1.0.2   2022-03-04 [1] CRAN (R 4.2.0)
# ROCR              1.0-11  2020-05-02 [1] CRAN (R 4.2.0)
# rpart             4.1.16  2022-01-24 [1] CRAN (R 4.2.0)
# rprojroot         2.0.3   2022-04-02 [1] CRAN (R 4.2.0)
# rstudioapi        0.13    2020-11-12 [1] CRAN (R 4.2.0)
# Rtsne             0.16    2022-04-17 [1] CRAN (R 4.2.0)
# rvest             1.0.2   2021-10-16 [1] CRAN (R 4.2.0)
# scales            1.2.0   2022-04-13 [1] CRAN (R 4.2.0)
# scattermore       0.8     2022-02-14 [1] CRAN (R 4.2.0)
# sctransform       0.3.3   2022-01-13 [1] CRAN (R 4.2.0)
# sessioninfo       1.2.2   2021-12-06 [1] CRAN (R 4.2.0)
# Seurat          * 4.1.1   2022-05-02 [1] CRAN (R 4.2.0)
# SeuratObject    * 4.1.0   2022-05-01 [1] CRAN (R 4.2.0)
# shiny             1.7.1   2021-10-02 [1] CRAN (R 4.2.0)
# sp              * 1.4-7   2022-04-20 [1] CRAN (R 4.2.0)
# spatstat.core     2.4-2   2022-04-01 [1] CRAN (R 4.2.0)
# spatstat.data     2.2-0   2022-04-18 [1] CRAN (R 4.2.0)
# spatstat.geom     2.4-0   2022-03-29 [1] CRAN (R 4.2.0)
# spatstat.random   2.2-0   2022-03-30 [1] CRAN (R 4.2.0)
# spatstat.sparse   2.1-1   2022-04-18 [1] CRAN (R 4.2.0)
# spatstat.utils    2.3-1   2022-05-06 [1] CRAN (R 4.2.0)
# stringi           1.7.6   2021-11-29 [1] CRAN (R 4.2.0)
# stringr         * 1.4.0   2019-02-10 [1] CRAN (R 4.2.0)
# survival          3.3-1   2022-03-03 [1] CRAN (R 4.2.0)
# tensor            1.5     2012-05-05 [1] CRAN (R 4.2.0)
# testthat          3.1.4   2022-04-26 [1] CRAN (R 4.2.0)
# tibble          * 3.1.7   2022-05-03 [1] CRAN (R 4.2.0)
# tidyr           * 1.2.0   2022-02-01 [1] CRAN (R 4.2.0)
# tidyselect        1.1.2   2022-02-21 [1] CRAN (R 4.2.0)
# tidyverse       * 1.3.1   2021-04-15 [1] CRAN (R 4.2.0)
# tzdb              0.3.0   2022-03-28 [1] CRAN (R 4.2.0)
# usethis           2.1.5   2021-12-09 [1] CRAN (R 4.2.0)
# utf8              1.2.2   2021-07-24 [1] CRAN (R 4.2.0)
# uwot              0.1.11  2021-12-02 [1] CRAN (R 4.2.0)
# vctrs             0.4.1   2022-04-13 [1] CRAN (R 4.2.0)
# viridis         * 0.6.2   2021-10-13 [1] CRAN (R 4.2.0)
# viridisLite     * 0.4.0   2021-04-13 [1] CRAN (R 4.2.0)
# withr             2.5.0   2022-03-03 [1] CRAN (R 4.2.0)
# xml2              1.3.3   2021-11-30 [1] CRAN (R 4.2.0)
# xtable            1.8-4   2019-04-21 [1] CRAN (R 4.2.0)
# zoo               1.8-10  2022-04-15 [1] CRAN (R 4.2.0)
# 
# [1] /Library/Frameworks/R.framework/Versions/4.2/Resources/library
# 
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────

