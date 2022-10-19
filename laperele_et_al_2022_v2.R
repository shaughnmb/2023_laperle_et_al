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
# R Script Title:  laperele_et_al_2022_v2.R
# R Script Author:  Shaughn Bell
# R Script Corresponding Email:  shaughn.bell@cshs.org
#
# Notes: 
#   A) Script makes use of the variables set up under "project information" as 
#      well as additional "prefixes" throughout the script for ease of saving
#      files with a similar path and naming structure.  When reloading data 
#      (see note "B"), you must either use the full path or reload the prefixs
#   B) Script saves intermediate steps at each major manipulation of the seurat
#      object.  These are not required for analysis, and they can be skipped to 
#      save time and disk space.
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
#
# Genesummed matrix files are included in the GEO record. 
#
##

#### project information ####

date <- "20221013"
project <- "CNS10_iNPC_gdnf"
datadir <- "/20221013_iNPC" # directory to save files generated
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

data.3812 <- read.csv(file=paste0(path_list[1],"/filtered_combined_mtx.csv"), header = TRUE)
data.4544 <- read.csv(file=paste0(path_list[2],"/filtered_combined_mtx.csv"), header = TRUE)

# put all the loaded matrices into one list
obj.list <- list(data.3812,data.4544)

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

data.3812 <- obj.list[[1]]
data.4544 <- obj.list[[2]]

# save genesummed matrix

write.csv(data.3812, file = paste0(path_list[1],"/filtered_genesummed_mtx.csv"), row.names = FALSE)
write.csv(data.4544, file = paste0(path_list[2],"/filtered_genesummed_mtx.csv"), row.names = FALSE)

# save R objects for easy reloading
save(data.3812,data.4544,
     file = paste0(datadir,"/",date,"_",project,"_genesummed.rdata"))

# load(file = paste0(datadir,"/",date,"_",project,"_genesummed.rdata"))

#### Seurat Object ####

# set wd to datadir
setwd(datadir)

# Initialize the Seurat object with the genesummed data.
# Keep all genes expressed in >= 1 cell

s3812 <- CreateSeuratObject(counts = obj.list[[1]], project = "3812", min.cells = 1, min.features = 0)
s4544 <- CreateSeuratObject(counts = obj.list[[2]], project = "4544", min.cells = 1, min.features = 0)

rm(data.3812,data.4544,
   dup.rows,dups,gene.duplicates,new.mat,
   obj.list,sub.rows,i,j,path_list)

save(s3812,s4544,
     file = paste0(datadir,"/",date,"_",project,"_seuratobjs.rdata"))

# load(file = paste0(datadir,"/",date,"_",project,"_seuratobjs.rdata"))

#### Add sample metadata ####

s4544[["Vial"]] <-"1024544"
s3812[["Vial"]] <-"1023812"

s4544[["Lot"]] <-"CNS10-GDNF"
s3812[["Lot"]] <-"iNPC-GDNF-WT"

s4544[["Passage"]] <-"p27"
s3812[["Passage"]] <-"p23"

#### Merge all data into a single seurat object ####

tenx1 <- merge(s4544, y = c(s3812),
               add.cell.ids = c("4544","3812"), 
               project = project)

tenx1

rm(s3812,s4544)

##save merged seurat object 
saveRDS(tenx1, file = paste0("./",date,"_",project,"_merged_seurat_prefilter.rds"))

# tenx1 <- read_rds(file = paste0("./",date,"_",project,"_merged_seurat_prefilter.rds"))

#### Initial QC Filtering ####

data <- tenx1 # save merged object in a temp "data" object

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
mean(data@meta.data$percent.ribo)
mean(data@meta.data$nCount_RNA)
median(data@meta.data$nCount_RNA)
mean(data@meta.data$nFeature_RNA)
median(data@meta.data$nFeature_RNA)
max(data@meta.data$nCount_RNA)

sink(paste0("./",date,"_",project,"_prefilter_QC_metrics.txt"))
cat("length")
length(data@meta.data$orig.ident)
cat("mean percent mito")
mean(data@meta.data$percent.mt)
cat("mean percent ribo")
mean(data@meta.data$percent.ribo)
cat("mean counts")
mean(data@meta.data$nCount_RNA)
cat("median counts")
median(data@meta.data$nCount_RNA)
cat("mean features")
mean(data@meta.data$nFeature_RNA)
cat("median features")
median(data@meta.data$nFeature_RNA)
cat("max counts")
max(data@meta.data$nCount_RNA)
sink()

data <- subset(data, subset = percent.mt.z < 3)
data <- subset(data, subset = percent.ribo.z < 3)
data <- subset(data, subset = nGene.z < 3)
data <- subset(data, subset = nUMI.z < 3)

length(data@meta.data$orig.ident)
mean(data@meta.data$percent.mt)
mean(data@meta.data$percent.ribo)
mean(data@meta.data$nCount_RNA)
median(data@meta.data$nCount_RNA)
mean(data@meta.data$nFeature_RNA)
median(data@meta.data$nFeature_RNA)
max(data@meta.data$nCount_RNA)

sink(paste0("./",date,"_",project,"_postfilter_QC_metrics.txt"))
cat("length")
length(data@meta.data$orig.ident)
cat("mean percent mito")
mean(data@meta.data$percent.mt)
cat("mean percent ribo")
mean(data@meta.data$percent.ribo)
cat("mean counts")
mean(data@meta.data$nCount_RNA)
cat("median counts")
median(data@meta.data$nCount_RNA)
cat("mean features")
mean(data@meta.data$nFeature_RNA)
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

rm(s.list,tenx1,anchors)

#### set factor levels ####

s.integrated$orig.ident <- factor(x = s.integrated$orig.ident, levels = c("4544","3812"))
s.integrated$Lot <- factor(x = s.integrated$Lot, levels = c("CNS10-GDNF","iNPC-GDNF-WT"))

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

# Note the calculated # of PCs
# meaningful.PCs <- 13

rm(sq.xx.coord,xx.coord,xx.gload,cut.off,
   eig,eig.percent,expected.contribution,features,
   i,scree,sum.eig,eigenvalues)

#### Run UMAP and look at UMAP plots ####

s.integrated <- RunUMAP(s.integrated, reduction = "pca", dims = 1:meaningful.PCs, verbose = TRUE)

# update prefixed variable
prefixPC <- paste0("./",date,"_",project,"_",meaningful.PCs,"PCs")

## UMAP plot by sample name ("orig.ident")
DimPlot(s.integrated, 
        reduction = "umap", 
        label = FALSE, 
        pt.size = .25, 
        group.by = "orig.ident")
hirestiff(paste0(prefixPC,"_UMAP_by_sample_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_by_sample_lowres.tiff"))

DimPlot(s.integrated, 
        reduction = "umap", 
        label = FALSE, 
        pt.size = .25, 
        group.by = "Lot")
hirestiff(paste0(prefixPC,"_UMAP_by_lot_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_by_lot_lowres.tiff"))

DimPlot(s.integrated, 
        reduction = "umap", 
        label = FALSE, 
        pt.size = .25, 
        split.by = "Lot") + 
  NoLegend()
hirestiff(paste0(prefixPC,"_UMAP_by_lot_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_by_lot_lowres.tiff"))

#### Cell Cycle Score ####

# not all genes will be found in the integrated object, so
# switch back to the RNA slot
DefaultAssay(s.integrated) <- "RNA"

s.integrated <- CellCycleScoring(s.integrated,
                               s.features = cc.genes.updated.2019$s.genes,
                               g2m.features = cc.genes.updated.2019$g2m.genes,
                               set.ident = TRUE)

as_tibble(s.integrated[[]]) %>%
  ggplot(aes(Phase)) + 
  geom_bar()

DimPlot(s.integrated, 
        reduction = "umap", 
        label = FALSE, 
        pt.size = .25, 
        group.by = "Phase")
hirestiff(paste0(prefixPC,"_UMAP_by_phase_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_by_phase_lowres.tiff"))

DimPlot(s.integrated, 
        reduction = "umap", 
        label = FALSE, 
        pt.size = .25, 
        split.by = "Lot", 
        group.by = "Phase")
hirestiff(paste0(prefixPC,"_UMAP_by_phase_splitby_orig.ident_hires.tiff"))
lowrestiff(paste0(prefixPC,"_UMAP_by_phase_splitby_orig.ident_lowres.tiff"))

# switch back to the integrated slot
DefaultAssay(s.integrated) <- "integrated"

saveRDS(s.integrated, file = paste0(prefixPC,"_seurat_integrated.rds"))

# s.integrated <- read_rds(file = paste0(prefixPC,"_seurat_integrated.rds"))

#### Clustering and Resolution ####

# Determine the K-nearest neighbor graph
s.integrated <- FindNeighbors(object = s.integrated, reduction = "pca", dims = 1:meaningful.PCs)

# Determine the clusters                              
s.integrated <- FindClusters(object = s.integrated,
                           resolution = c(0.5))

res <- "_res.0.5"
ident_res <- paste0("integrated_snn",res)

Idents(s.integrated) <- ident_res
DimPlot(s.integrated, reduction = "umap", label = TRUE, label.size = 5) + 
  NoLegend() + 
  ggtitle(paste0(ident_res))

#update prefix
prefixPCres <- paste0(prefixPC,res)

# Plot the UMAP
DimPlot(s.integrated, 
        reduction = "umap", 
        label = TRUE, 
        label.size = 5) + 
  NoLegend() + 
  ggtitle(paste0(ident_res))
hirestiff(paste0(prefixPCres,"_UMAP","_by_","cluster","_hires.tiff"))
lowrestiff(paste0(prefixPCres,"_UMAP","_by_","cluster","_lowres.tiff"))

# UMAP of cells in each cluster by lot
DimPlot(s.integrated, 
        reduction = "umap", 
        label = FALSE, 
        split.by = "Lot") + 
  NoLegend() + 
  ggtitle(paste0(ident_res))
hirestiff(paste0(prefixPCres,"_UMAP","_by_","lot","_hires.tiff"))
lowrestiff(paste0(prefixPCres,"_UMAP","_by_","lot","_lowres.tiff"))

# UMAP of cells in each cluster by Lot with cluster lables
DimPlot(s.integrated, 
        reduction = "umap", 
        label = TRUE, 
        split.by = "Lot") + 
  NoLegend() + 
  ggtitle(paste0(ident_res))
hirestiff(paste0(prefixPCres,"_UMAP","_by_","lot","_hires.tiff"))
lowrestiff(paste0(prefixPCres,"_UMAP","_by_","lot","_lowres.tiff"))

# fix cluster colors
clust_cols <- c("0"="#F8766D","1"="#DB8E00","2"="#AEA200",
                "3"="#64B200","4"="#00BD5C","5"="#00C1A7",
                "6"="#00BADE","7"="#00A6FF","8"="#B385FF",
                "9"="#EF67EB","10"="#FF63B6")

# for loop for saving split plots in the same dimensions as composite plot
# using the fixed cluster color scheme

for (i in 1:length(levels(s.integrated$Lot))){
  p <- (DimPlot(s.integrated, 
                reduction = "umap", 
                label = FALSE, 
                pt.size = .25, 
                cells = c(WhichCells(s.integrated, expression = Lot == levels(s.integrated$Lot)[i])),
                cols = clust_cols) + 
          NoLegend() +
          ggtitle(paste0(levels(s.integrated$Lot)[i])))
  
  print(p)
  
  hirestiff(paste(prefixPCres,"UMAP","split_by","lot","no_cluster_labels",
                  levels(s.integrated$Lot)[i],"hires.tiff",sep = "_"))
  lowrestiff(paste(prefixPCres,"UMAP","split_by","lot","no_cluster_labels",
                   levels(s.integrated$Lot)[i],"lowres.tiff",sep = "_"))
}

# save RDS containing reduction and cluster idents
saveRDS(s.integrated, paste0(prefixPCres,"_seurat_after_clustering.rds"))

# s.integrated <- readRDS(file = paste0(prefixPCres,"_seurat_after_clustering.rds"))

#### Extract number of cells per cluster per Lot ####
n_cells <- FetchData(s.integrated, vars = c("ident", "Lot")) %>%
  dplyr::count(ident, Lot) %>%
  tidyr::spread(ident, n)

write.csv(n_cells, file = paste0(prefixPCres,"_cells_per_cluster.csv"))

rm(n_cells,p,i)

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
DimPlot(s.integrated, 
        reduction = "umap", 
        label = TRUE, 
        label.size = 5)  + 
  NoLegend()
  
hirestiff(paste0(prefixPCres,"_Fig","_1b_","_hires.tiff"))
lowrestiff(paste0(prefixPCres,"_Fig","_1b_","_lowres.tiff"))

#### Figure 1c ####
# for loop that saves UMAP of each Lot by cluster without cluster labels 
# with the same dimensions as composite plot using the fixed cluster color scheme

# fix cluster colors
clust_cols <- c("0"="#F8766D","1"="#DB8E00","2"="#AEA200",
                "3"="#64B200","4"="#00BD5C","5"="#00C1A7",
                "6"="#00BADE","7"="#00A6FF","8"="#B385FF",
                "9"="#EF67EB","10"="#FF63B6")

# for loop
for (i in 1:length(levels(s.integrated$Lot))){
  p <- (DimPlot(s.integrated, 
                reduction = "umap", 
                label = FALSE, 
                pt.size = .25, 
                cells = c(WhichCells(s.integrated, expression = Lot == levels(s.integrated$Lot)[i])),
                cols = clust_cols) + 
          NoLegend() +
          ggtitle(paste0(levels(s.integrated$Lot)[i])))
  
  print(p)
  
  hirestiff(paste(prefixPCres,"Fig","1c",
                  levels(s.integrated$Lot)[i],"hires.tiff",sep = "_"))
  lowrestiff(paste(prefixPCres,"Fig","1c",
                   levels(s.integrated$Lot)[i],"lowres.tiff",sep = "_"))
}

#### Figure 1d ####

FeaturePlot(s.integrated, 
            reduction = "umap", 
            features = c("S100B","VIM","CD99","NES"),
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE,
            label.size = 4,
            pt.size = 0.3)

hirestiff(paste(prefixPCres,"Fig","1d","hires.tiff",sep = "_"))
lowrestiff(paste(prefixPCres,"Fig","1d","lowres.tiff",sep = "_"))

#### Figure 1e ####

FeaturePlot(s.integrated, 
            reduction = "umap", 
            features = c("GFAP","AQP4","SLC1A3"),
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE,
            label.size = 4,
            pt.size = 0.3)

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
            label.size = 4,
            pt.size = 0.3)

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
 
# sessionInfo()
# R version 4.2.1 (2022-06-23)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Monterey 12.6
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
#   [1] Vennerable_3.1.0.9000 gridExtra_2.3         irlba_2.3.5           patchwork_1.1.2       cowplot_1.1.1        
# [6] Matrix_1.5-0          forcats_0.5.2         stringr_1.4.1         dplyr_1.0.10          purrr_0.3.4          
# [11] readr_2.1.2           tidyr_1.2.1           tibble_3.1.8          ggplot2_3.3.6         tidyverse_1.3.2      
# [16] sp_1.5-0              SeuratObject_4.1.1    Seurat_4.1.1         
# 
# loaded via a namespace (and not attached):
#   [1] utf8_1.2.2            reticulate_1.26       tidyselect_1.1.2      htmlwidgets_1.5.4     grid_4.2.1           
# [6] Rtsne_0.16            pROC_1.18.0           devtools_2.4.4        munsell_0.5.0         ragg_1.2.2           
# [11] codetools_0.2-18      ica_1.0-3             interp_1.1-3          future_1.28.0         miniUI_0.1.1.1       
# [16] withr_2.5.0           spatstat.random_2.2-0 colorspace_2.0-3      progressr_0.11.0      knitr_1.40           
# [21] rstudioapi_0.14       stats4_4.2.1          ROCR_1.0-11           tensor_1.5            listenv_0.8.0        
# [26] labeling_0.4.2        polyclip_1.10-0       farver_2.1.1          parallelly_1.32.1     vctrs_0.4.1          
# [31] generics_0.1.3        xfun_0.32             ipred_0.9-13          randomForest_4.7-1.1  R6_2.5.1             
# [36] doParallel_1.0.17     ggbeeswarm_0.6.0      clue_0.3-61           bitops_1.0-7          spatstat.utils_2.3-1 
# [41] cachem_1.0.6          assertthat_0.2.1      promises_1.2.0.1      scales_1.2.1          googlesheets4_1.0.1  
# [46] nnet_7.3-17           beeswarm_0.4.0        rgeos_0.5-9           gtable_0.3.1          globals_0.16.1       
# [51] processx_3.7.0        goftest_1.2-3         timeDate_4021.104     rlang_1.0.5           systemfonts_1.0.4    
# [56] GlobalOptions_0.1.2   splines_4.2.1         lazyeval_0.2.2        gargle_1.2.1          ModelMetrics_1.2.2.2 
# [61] broom_1.0.1           checkmate_2.1.0       spatstat.geom_2.4-0   modelr_0.1.9          BiocManager_1.30.18  
# [66] reshape2_1.4.4        abind_1.4-5           backports_1.4.1       httpuv_1.6.6          Hmisc_4.7-1          
# [71] RBGL_1.72.0           caret_6.0-93          DiagrammeR_1.0.9      tools_4.2.1           lava_1.6.10          
# [76] usethis_2.1.6         ellipsis_0.3.2        spatstat.core_2.4-4   RColorBrewer_1.1-3    proxy_0.4-27         
# [81] BiocGenerics_0.42.0   sessioninfo_1.2.2     ggridges_0.5.3        Rcpp_1.0.9            plyr_1.8.7           
# [86] base64enc_0.1-3       visNetwork_2.1.0      ps_1.7.1              prettyunits_1.1.1     rpart_4.1.16         
# [91] deldir_1.0-6          pbapply_1.5-0         GetoptLong_1.0.5      urlchecker_1.0.1      S4Vectors_0.34.0     
# [96] zoo_1.8-10            haven_2.5.1           nichenetr_1.1.0       ggrepel_0.9.1         cluster_2.1.4        
# [101] fs_1.5.2              magrittr_2.0.3        data.table_1.14.2     scattermore_0.8       circlize_0.4.15      
# [106] reprex_2.0.2          lmtest_0.9-40         RANN_2.6.1            googledrive_2.0.0     fitdistrplus_1.1-8   
# [111] matrixStats_0.62.0    pkgload_1.3.0         hms_1.1.2             mime_0.12             xtable_1.8-4         
# [116] jpeg_0.1-9            readxl_1.4.1          IRanges_2.30.1        proftools_0.99-3      shape_1.4.6          
# [121] compiler_4.2.1        KernSmooth_2.23-20    crayon_1.5.1          htmltools_0.5.3       mgcv_1.8-40          
# [126] later_1.3.0           tzdb_0.3.0            Formula_1.2-4         lubridate_1.8.0       DBI_1.1.3            
# [131] dbplyr_2.2.1          ComplexHeatmap_2.12.1 MASS_7.3-58.1         cli_3.4.0             parallel_4.2.1       
# [136] gower_1.0.0           igraph_1.3.4          pkgconfig_2.0.3       foreign_0.8-82        plotly_4.10.0        
# [141] spatstat.sparse_2.1-1 recipes_1.0.1         xml2_1.3.3            foreach_1.5.2         vipor_0.4.5          
# [146] hardhat_1.2.0         prodlim_2019.11.13    rvest_1.0.3           NCmisc_1.1.6          callr_3.7.2          
# [151] digest_0.6.29         sctransform_0.3.4     RcppAnnoy_0.0.19      graph_1.74.0          spatstat.data_2.2-0  
# [156] cellranger_1.1.0      leiden_0.4.3          htmlTable_2.4.1       uwot_0.1.14           shiny_1.7.2          
# [161] rjson_0.2.21          lifecycle_1.0.2       nlme_3.1-159          jsonlite_1.8.0        limma_3.52.2         
# [166] viridisLite_0.4.1     fansi_1.0.3           pillar_1.8.1          lattice_0.20-45       ggrastr_1.0.1        
# [171] fastmap_1.1.0         httr_1.4.4            pkgbuild_1.3.1        survival_3.4-0        glue_1.6.2           
# [176] remotes_2.4.2         fdrtool_1.2.17        png_0.1-7             iterators_1.0.14      class_7.3-20         
# [181] stringi_1.7.8         profvis_0.3.7         textshaping_0.3.6     caTools_1.18.2        latticeExtra_0.6-30  
# [186] memoise_2.0.1         e1071_1.7-11          future.apply_1.9.1  

#### devtools session info ####

# devtools::session_info()
# ─ Session info ──────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.2.1 (2022-06-23)
# os       macOS Monterey 12.6
# system   x86_64, darwin17.0
# ui       RStudio
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       America/Los_Angeles
# date     2022-10-13
# rstudio  2022.07.1+554 Spotted Wakerobin (desktop)
# pandoc   NA
# 
# ─ Packages ──────────────────────────────────────────────────────────────────────────────────────────────────────────
# package         * version    date (UTC) lib source
# abind             1.4-5      2016-07-21 [1] CRAN (R 4.2.0)
# assertthat        0.2.1      2019-03-21 [1] CRAN (R 4.2.0)
# backports         1.4.1      2021-12-13 [1] CRAN (R 4.2.0)
# base64enc         0.1-3      2015-07-28 [1] CRAN (R 4.2.0)
# beeswarm          0.4.0      2021-06-01 [1] CRAN (R 4.2.0)
# BiocGenerics      0.42.0     2022-04-26 [1] Bioconductor
# BiocManager       1.30.18    2022-05-18 [1] CRAN (R 4.2.0)
# bitops            1.0-7      2021-04-24 [1] CRAN (R 4.2.0)
# broom             1.0.1      2022-08-29 [1] CRAN (R 4.2.0)
# cachem            1.0.6      2021-08-19 [1] CRAN (R 4.2.0)
# callr             3.7.2      2022-08-22 [1] CRAN (R 4.2.0)
# caret             6.0-93     2022-08-09 [1] CRAN (R 4.2.0)
# caTools           1.18.2     2021-03-28 [1] CRAN (R 4.2.0)
# cellranger        1.1.0      2016-07-27 [1] CRAN (R 4.2.0)
# checkmate         2.1.0      2022-04-21 [1] CRAN (R 4.2.0)
# circlize          0.4.15     2022-05-10 [1] CRAN (R 4.2.0)
# class             7.3-20     2022-01-16 [1] CRAN (R 4.2.1)
# cli               3.4.0      2022-09-08 [1] CRAN (R 4.2.0)
# clue              0.3-61     2022-05-30 [1] CRAN (R 4.2.0)
# cluster           2.1.4      2022-08-22 [1] CRAN (R 4.2.0)
# codetools         0.2-18     2020-11-04 [1] CRAN (R 4.2.1)
# colorspace        2.0-3      2022-02-21 [1] CRAN (R 4.2.0)
# ComplexHeatmap    2.12.1     2022-08-09 [1] Bioconductor
# cowplot         * 1.1.1      2020-12-30 [1] CRAN (R 4.2.0)
# crayon            1.5.1      2022-03-26 [1] CRAN (R 4.2.0)
# data.table        1.14.2     2021-09-27 [1] CRAN (R 4.2.0)
# DBI               1.1.3      2022-06-18 [1] CRAN (R 4.2.0)
# dbplyr            2.2.1      2022-06-27 [1] CRAN (R 4.2.0)
# deldir            1.0-6      2021-10-23 [1] CRAN (R 4.2.0)
# devtools          2.4.4      2022-07-20 [1] CRAN (R 4.2.1)
# DiagrammeR        1.0.9      2022-03-05 [1] CRAN (R 4.2.0)
# digest            0.6.29     2021-12-01 [1] CRAN (R 4.2.0)
# doParallel        1.0.17     2022-02-07 [1] CRAN (R 4.2.0)
# dplyr           * 1.0.10     2022-09-01 [1] CRAN (R 4.2.0)
# e1071             1.7-11     2022-06-07 [1] CRAN (R 4.2.0)
# ellipsis          0.3.2      2021-04-29 [1] CRAN (R 4.2.0)
# fansi             1.0.3      2022-03-24 [1] CRAN (R 4.2.0)
# farver            2.1.1      2022-07-06 [1] CRAN (R 4.2.0)
# fastmap           1.1.0      2021-01-25 [1] CRAN (R 4.2.0)
# fdrtool           1.2.17     2021-11-13 [1] CRAN (R 4.2.0)
# fitdistrplus      1.1-8      2022-03-10 [1] CRAN (R 4.2.0)
# forcats         * 0.5.2      2022-08-19 [1] CRAN (R 4.2.0)
# foreach           1.5.2      2022-02-02 [1] CRAN (R 4.2.0)
# foreign           0.8-82     2022-01-16 [1] CRAN (R 4.2.1)
# Formula           1.2-4      2020-10-16 [1] CRAN (R 4.2.0)
# fs                1.5.2      2021-12-08 [1] CRAN (R 4.2.0)
# future            1.28.0     2022-09-02 [1] CRAN (R 4.2.0)
# future.apply      1.9.1      2022-09-07 [1] CRAN (R 4.2.0)
# gargle            1.2.1      2022-09-08 [1] CRAN (R 4.2.0)
# generics          0.1.3      2022-07-05 [1] CRAN (R 4.2.0)
# GetoptLong        1.0.5      2020-12-15 [1] CRAN (R 4.2.0)
# ggbeeswarm        0.6.0      2017-08-07 [1] CRAN (R 4.2.0)
# ggplot2         * 3.3.6      2022-05-03 [1] CRAN (R 4.2.0)
# ggrastr           1.0.1      2021-12-08 [1] CRAN (R 4.2.0)
# ggrepel           0.9.1      2021-01-15 [1] CRAN (R 4.2.0)
# ggridges          0.5.3      2021-01-08 [1] CRAN (R 4.2.0)
# GlobalOptions     0.1.2      2020-06-10 [1] CRAN (R 4.2.0)
# globals           0.16.1     2022-08-28 [1] CRAN (R 4.2.0)
# glue              1.6.2      2022-02-24 [1] CRAN (R 4.2.0)
# goftest           1.2-3      2021-10-07 [1] CRAN (R 4.2.0)
# googledrive       2.0.0      2021-07-08 [1] CRAN (R 4.2.0)
# googlesheets4     1.0.1      2022-08-13 [1] CRAN (R 4.2.0)
# gower             1.0.0      2022-02-03 [1] CRAN (R 4.2.0)
# graph             1.74.0     2022-04-26 [1] Bioconductor
# gridExtra       * 2.3        2017-09-09 [1] CRAN (R 4.2.0)
# gtable            0.3.1      2022-09-01 [1] CRAN (R 4.2.0)
# hardhat           1.2.0      2022-06-30 [1] CRAN (R 4.2.0)
# haven             2.5.1      2022-08-22 [1] CRAN (R 4.2.0)
# Hmisc             4.7-1      2022-08-15 [1] CRAN (R 4.2.0)
# hms               1.1.2      2022-08-19 [1] CRAN (R 4.2.0)
# htmlTable         2.4.1      2022-07-07 [1] CRAN (R 4.2.1)
# htmltools         0.5.3      2022-07-18 [1] CRAN (R 4.2.0)
# htmlwidgets       1.5.4      2021-09-08 [1] CRAN (R 4.2.0)
# httpuv            1.6.6      2022-09-08 [1] CRAN (R 4.2.0)
# httr              1.4.4      2022-08-17 [1] CRAN (R 4.2.0)
# ica               1.0-3      2022-07-08 [1] CRAN (R 4.2.1)
# igraph            1.3.4      2022-07-19 [1] CRAN (R 4.2.0)
# interp            1.1-3      2022-07-13 [1] CRAN (R 4.2.0)
# ipred             0.9-13     2022-06-02 [1] CRAN (R 4.2.0)
# IRanges           2.30.1     2022-08-18 [1] Bioconductor
# irlba           * 2.3.5      2021-12-06 [1] CRAN (R 4.2.0)
# iterators         1.0.14     2022-02-05 [1] CRAN (R 4.2.0)
# jpeg              0.1-9      2021-07-24 [1] CRAN (R 4.2.0)
# jsonlite          1.8.0      2022-02-22 [1] CRAN (R 4.2.0)
# KernSmooth        2.23-20    2021-05-03 [1] CRAN (R 4.2.1)
# knitr             1.40       2022-08-24 [1] CRAN (R 4.2.1)
# labeling          0.4.2      2020-10-20 [1] CRAN (R 4.2.0)
# later             1.3.0      2021-08-18 [1] CRAN (R 4.2.0)
# lattice           0.20-45    2021-09-22 [1] CRAN (R 4.2.1)
# latticeExtra      0.6-30     2022-07-04 [1] CRAN (R 4.2.0)
# lava              1.6.10     2021-09-02 [1] CRAN (R 4.2.0)
# lazyeval          0.2.2      2019-03-15 [1] CRAN (R 4.2.0)
# leiden            0.4.3      2022-09-10 [1] CRAN (R 4.2.0)
# lifecycle         1.0.2      2022-09-09 [1] CRAN (R 4.2.0)
# limma             3.52.2     2022-06-19 [1] Bioconductor
# listenv           0.8.0      2019-12-05 [1] CRAN (R 4.2.0)
# lmtest            0.9-40     2022-03-21 [1] CRAN (R 4.2.0)
# lubridate         1.8.0      2021-10-07 [1] CRAN (R 4.2.0)
# magrittr          2.0.3      2022-03-30 [1] CRAN (R 4.2.0)
# MASS              7.3-58.1   2022-08-03 [1] CRAN (R 4.2.0)
# Matrix          * 1.5-0      2022-09-10 [1] CRAN (R 4.2.0)
# matrixStats       0.62.0     2022-04-19 [1] CRAN (R 4.2.0)
# memoise           2.0.1      2021-11-26 [1] CRAN (R 4.2.0)
# mgcv              1.8-40     2022-03-29 [1] CRAN (R 4.2.1)
# mime              0.12       2021-09-28 [1] CRAN (R 4.2.0)
# miniUI            0.1.1.1    2018-05-18 [1] CRAN (R 4.2.0)
# ModelMetrics      1.2.2.2    2020-03-17 [1] CRAN (R 4.2.0)
# modelr            0.1.9      2022-08-19 [1] CRAN (R 4.2.0)
# munsell           0.5.0      2018-06-12 [1] CRAN (R 4.2.0)
# NCmisc            1.1.6      2018-11-12 [1] CRAN (R 4.2.0)
# nichenetr         1.1.0      2022-05-25 [1] Github (saeyslab/nichenetr@ee3263c)
# nlme              3.1-159    2022-08-09 [1] CRAN (R 4.2.0)
# nnet              7.3-17     2022-01-16 [1] CRAN (R 4.2.1)
# parallelly        1.32.1     2022-07-21 [1] CRAN (R 4.2.1)
# patchwork       * 1.1.2      2022-08-19 [1] CRAN (R 4.2.0)
# pbapply           1.5-0      2021-09-16 [1] CRAN (R 4.2.0)
# pillar            1.8.1      2022-08-19 [1] CRAN (R 4.2.0)
# pkgbuild          1.3.1      2021-12-20 [1] CRAN (R 4.2.0)
# pkgconfig         2.0.3      2019-09-22 [1] CRAN (R 4.2.0)
# pkgload           1.3.0      2022-06-27 [1] CRAN (R 4.2.0)
# plotly            4.10.0     2021-10-09 [1] CRAN (R 4.2.0)
# plyr              1.8.7      2022-03-24 [1] CRAN (R 4.2.0)
# png               0.1-7      2013-12-03 [1] CRAN (R 4.2.0)
# polyclip          1.10-0     2019-03-14 [1] CRAN (R 4.2.0)
# prettyunits       1.1.1      2020-01-24 [1] CRAN (R 4.2.0)
# pROC              1.18.0     2021-09-03 [1] CRAN (R 4.2.0)
# processx          3.7.0      2022-07-07 [1] CRAN (R 4.2.1)
# prodlim           2019.11.13 2019-11-17 [1] CRAN (R 4.2.0)
# proftools         0.99-3     2020-07-08 [1] CRAN (R 4.2.0)
# profvis           0.3.7      2020-11-02 [1] CRAN (R 4.2.0)
# progressr         0.11.0     2022-09-02 [1] CRAN (R 4.2.0)
# promises          1.2.0.1    2021-02-11 [1] CRAN (R 4.2.0)
# proxy             0.4-27     2022-06-09 [1] CRAN (R 4.2.0)
# ps                1.7.1      2022-06-18 [1] CRAN (R 4.2.0)
# purrr           * 0.3.4      2020-04-17 [1] CRAN (R 4.2.0)
# R6                2.5.1      2021-08-19 [1] CRAN (R 4.2.0)
# ragg              1.2.2      2022-02-21 [1] CRAN (R 4.2.0)
# randomForest      4.7-1.1    2022-05-23 [1] CRAN (R 4.2.0)
# RANN              2.6.1      2019-01-08 [1] CRAN (R 4.2.0)
# RBGL              1.72.0     2022-04-26 [1] Bioconductor
# RColorBrewer      1.1-3      2022-04-03 [1] CRAN (R 4.2.0)
# Rcpp              1.0.9      2022-07-08 [1] CRAN (R 4.2.0)
# RcppAnnoy         0.0.19     2021-07-30 [1] CRAN (R 4.2.0)
# readr           * 2.1.2      2022-01-30 [1] CRAN (R 4.2.0)
# readxl            1.4.1      2022-08-17 [1] CRAN (R 4.2.0)
# recipes           1.0.1      2022-07-07 [1] CRAN (R 4.2.1)
# remotes           2.4.2      2021-11-30 [1] CRAN (R 4.2.0)
# reprex            2.0.2      2022-08-17 [1] CRAN (R 4.2.0)
# reshape2          1.4.4      2020-04-09 [1] CRAN (R 4.2.0)
# reticulate        1.26       2022-08-31 [1] CRAN (R 4.2.0)
# rgeos             0.5-9      2021-12-15 [1] CRAN (R 4.2.0)
# rjson             0.2.21     2022-01-09 [1] CRAN (R 4.2.0)
# rlang             1.0.5      2022-08-31 [1] CRAN (R 4.2.0)
# ROCR              1.0-11     2020-05-02 [1] CRAN (R 4.2.0)
# rpart             4.1.16     2022-01-24 [1] CRAN (R 4.2.1)
# rstudioapi        0.14       2022-08-22 [1] CRAN (R 4.2.0)
# Rtsne             0.16       2022-04-17 [1] CRAN (R 4.2.0)
# rvest             1.0.3      2022-08-19 [1] CRAN (R 4.2.0)
# S4Vectors         0.34.0     2022-04-26 [1] Bioconductor
# scales            1.2.1      2022-08-20 [1] CRAN (R 4.2.0)
# scattermore       0.8        2022-02-14 [1] CRAN (R 4.2.0)
# sctransform       0.3.4      2022-08-20 [1] CRAN (R 4.2.0)
# sessioninfo       1.2.2      2021-12-06 [1] CRAN (R 4.2.0)
# Seurat          * 4.1.1      2022-05-02 [1] CRAN (R 4.2.0)
# SeuratObject    * 4.1.1      2022-08-29 [1] CRAN (R 4.2.0)
# shape             1.4.6      2021-05-19 [1] CRAN (R 4.2.0)
# shiny             1.7.2      2022-07-19 [1] CRAN (R 4.2.0)
# sp              * 1.5-0      2022-06-05 [1] CRAN (R 4.2.0)
# spatstat.core     2.4-4      2022-05-18 [1] CRAN (R 4.2.0)
# spatstat.data     2.2-0      2022-04-18 [1] CRAN (R 4.2.0)
# spatstat.geom     2.4-0      2022-03-29 [1] CRAN (R 4.2.0)
# spatstat.random   2.2-0      2022-03-30 [1] CRAN (R 4.2.0)
# spatstat.sparse   2.1-1      2022-04-18 [1] CRAN (R 4.2.0)
# spatstat.utils    2.3-1      2022-05-06 [1] CRAN (R 4.2.0)
# stringi           1.7.8      2022-07-11 [1] CRAN (R 4.2.0)
# stringr         * 1.4.1      2022-08-20 [1] CRAN (R 4.2.0)
# survival          3.4-0      2022-08-09 [1] CRAN (R 4.2.0)
# systemfonts       1.0.4      2022-02-11 [1] CRAN (R 4.2.0)
# tensor            1.5        2012-05-05 [1] CRAN (R 4.2.0)
# textshaping       0.3.6      2021-10-13 [1] CRAN (R 4.2.0)
# tibble          * 3.1.8      2022-07-22 [1] CRAN (R 4.2.0)
# tidyr           * 1.2.1      2022-09-08 [1] CRAN (R 4.2.0)
# tidyselect        1.1.2      2022-02-21 [1] CRAN (R 4.2.0)
# tidyverse       * 1.3.2      2022-07-18 [1] CRAN (R 4.2.0)
# timeDate          4021.104   2022-07-19 [1] CRAN (R 4.2.0)
# tzdb              0.3.0      2022-03-28 [1] CRAN (R 4.2.0)
# urlchecker        1.0.1      2021-11-30 [1] CRAN (R 4.2.0)
# usethis           2.1.6      2022-05-25 [1] CRAN (R 4.2.0)
# utf8              1.2.2      2021-07-24 [1] CRAN (R 4.2.0)
# uwot              0.1.14     2022-08-22 [1] CRAN (R 4.2.0)
# vctrs             0.4.1      2022-04-13 [1] CRAN (R 4.2.0)
# Vennerable      * 3.1.0.9000 2022-05-10 [1] Github (js229/Vennerable@46057c9)
# vipor             0.4.5      2017-03-22 [1] CRAN (R 4.2.0)
# viridisLite       0.4.1      2022-08-22 [1] CRAN (R 4.2.0)
# visNetwork        2.1.0      2021-09-29 [1] CRAN (R 4.2.0)
# withr             2.5.0      2022-03-03 [1] CRAN (R 4.2.0)
# xfun              0.32       2022-08-10 [1] CRAN (R 4.2.0)
# xml2              1.3.3      2021-11-30 [1] CRAN (R 4.2.0)
# xtable            1.8-4      2019-04-21 [1] CRAN (R 4.2.0)
# zoo               1.8-10     2022-04-15 [1] CRAN (R 4.2.0)
# 
# [1] /Library/Frameworks/R.framework/Versions/4.2/Resources/library
# 
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
