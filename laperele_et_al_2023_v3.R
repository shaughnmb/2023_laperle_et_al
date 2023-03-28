#!/bin/env Rscript

#### License Notice ####

##
# Copyright (c) 2023 Cedars-Sinai Medical Center
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
# Title: Human iPSC-derived neural progenitor cells secreting GDNF provide 
#        protection in rodent models of ALS and retinal degeneration 
# Authors:  Alexander H. Laperle*1, Alexandra Moser*1, Pablo Avalos 1, Bin Lu 1, 
#           Amanda Wu 1, Aaron Fulton 1, Stephany Ramirez 1, Veronica J. Garcia 1, 
#           Shaughn Bell 1, Ritchie Ho 1, George Lawless 1, Kristina Roxas 1, 
#           Saba Shahin 1, Oksana Shelest 1, Soshana Svendsen 1, Shaomei Wang 1, 
#           Clive N. Svendsen**1 
#
# Affiliations: 1 Cedars-Sinai Board of Governors Regenerative Medicine 
#                 Institute, Cedars-Sinai Medical Center, Los Angeles CA
# *These authors contributed equally to this work 
# **Corresponding author email:  Clive.Svendsen@cshs.org 
##

#### Script Information ####

##
# R version 4.2.2
# R Script Title:  2023_laperle_et_al_v3.R
# R Script Author:  Shaughn Bell
# R Script Corresponding Email:  Shaughn.Bell@cshs.org
#
# Notes: 
#   A) Script makes use of the variables set up under "project information" as 
#      well as additional "prefixes" throughout the script for ease of saving
#      files with a similar path and naming structure.  When reloading data 
#      (see note "B"), you must either use the full path or reload the prefixes.
#   B) Script saves intermediate steps at each major manipulation of the seurat
#      object via the "saveRDS" function.  If needed, these RDS objects can then 
#      be reloaded to easily restart at one of these save points without needing 
#      to start from scratch.  However, these are not required for analysis, and
#      they can be skipped to save time and disk space.
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

hireseps <- function(saveas){
  ggsave(
    saveas,
    plot = last_plot(),
    device = "eps",
    scale = 1,
    dpi = 600,
    limitsize = TRUE,
    bg = "white",
    width = 4421,
    height = 3221,
    units = c("px")
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

# check that all the files exist
all(file.exists(path_list)) 

#Read the cellranger output files and save expression matrix csv
for (k in 1:length(path_list)) {
  
  # read in and transpose the barcode file
  bc <- as.data.frame(
    read.delim(
      file = paste0(
        path_list[k], "/barcodes.tsv.gz"
      ),
      header = F
    )
  ) %>% t(.)

  # read in the features file and select just the first two columns
  feat <- as.data.frame(
    read.delim(
      file = paste0(
        path_list[k], "/features.tsv.gz"
      ),
      header = F
    )
  ) %>% .[, 1:2]

  # read in the counts
  matx <- as.data.frame(
    readMM(
      file = paste0(
        path_list[k], "/matrix.mtx.gz"
      )
    )
  )

  # assign barcodes as colnames for the matrix
  colnames(matx) <- bc[1, ]

  # rename the columns of features table
  colnames(feat) <- c("ensmbl.id", "symbol")

  # cbind the features table and the matrix
  matx <- cbind(feat, matx)

  write.csv(matx,
    file = paste0(path_list[k], "/filtered_combined_mtx.csv"),
    row.names = FALSE
  )
}

rm(matx, bc, feat, filelist, k)

# check the order of the files in path_list
path_list

# Load the data
# Make sure the index number after path_list matches the correct file

data.3812 <- read.csv(
  file = paste0(
    path_list[1],
    "/filtered_combined_mtx.csv"
  ),
  header = TRUE
)
data.4544 <- read.csv(
  file = paste0(
    path_list[2],
    "/filtered_combined_mtx.csv"
  ),
  header = TRUE
)

# Put all the loaded matrices into one list
obj.list <- list(data.3812, data.4544)

# Iterate the summing over all samples in obj.list
for (j in 1:length(obj.list)) {

    # Create a list of duplicated genes
  gene.duplicates <- obj.list[[j]]$symbol[
    duplicated(obj.list[[j]]$symbol, incomparables = NA)
  ] %>%
    as.list(.) %>%
    unique(.)

  # Create empty dataframe to rbind duplicated rows to
  dup.rows <- data.frame()

  # Create identical expression matrix to subtract duplicated rows from
  sub.rows <- obj.list[[j]]

  # Extract rows containing the duplicated genes.
  for (i in 1:length(gene.duplicates)) {
    # subsets duplicated rows
    dups <- subset(obj.list[[j]], subset = symbol == gene.duplicates[[i]])
    # binds duplicated rows together
    dup.rows <- rbind(dup.rows, dups)
    # deletes duplicated rows from expression matrix
    sub.rows <- subset(sub.rows, subset = symbol != gene.duplicates[[i]])
  }

  # Remove ensembl ID columns
  dup.rows <- dup.rows[, 2:length(colnames(dup.rows))]
  sub.rows <- sub.rows[, 2:length(colnames(sub.rows))]

  # Find the gsum of duplicated rows
  dup.rows <- aggregate(. ~ symbol, data = dup.rows, sum)

  # rbind the gsum of duplicated rows with the 
  # expression matrix that had these genes removed
  new.mat <- rbind(sub.rows, dup.rows)

  if (nrow(new.mat) == length(unique(obj.list[[j]]$symbol))) {
    new.mat <- as.data.frame(new.mat)
    rownames(new.mat) <- new.mat[, 1]
    new.mat <- subset(new.mat, select = -c(symbol))
    obj.list[[j]] <- new.mat
    message("Successfully merged!")
  } else {
    message("There's an error!")
  }
}

data.3812 <- obj.list[[1]]
data.4544 <- obj.list[[2]]

# save genesummed matrix

write.csv(data.3812,
  file = paste0(path_list[1], "/filtered_genesummed_mtx.csv"),
  row.names = FALSE
)
write.csv(data.4544,
  file = paste0(path_list[2], "/filtered_genesummed_mtx.csv"),
  row.names = FALSE
)

# save R objects for easy reloading
save(data.3812, data.4544,
  file = paste0(datadir, "/", date, "_", project, "_genesummed.rdata")
)

# load(file = paste0(datadir,"/",date,"_",project,"_genesummed.rdata"))

#### Seurat Object ####

# set wd to datadir
setwd(datadir)

# Initialize the Seurat object with the genesummed data.
# Keep all genes expressed in >= 1 cell

s3812 <- CreateSeuratObject(
  counts = obj.list[[1]],
  project = "3812",
  min.cells = 1,
  min.features = 0
)

s4544 <- CreateSeuratObject(
  counts = obj.list[[2]],
  project = "4544",
  min.cells = 1,
  min.features = 0
)

rm(
  data.3812, data.4544,
  dup.rows, dups, gene.duplicates, new.mat,
  obj.list, sub.rows, i, j, path_list
)

save(s3812, s4544,
  file = paste0(datadir, "/", date, "_", project, "_seuratobjs.rdata")
)

# load(file = paste0(datadir,"/",date,"_",project,"_seuratobjs.rdata"))

#### Add sample metadata ####

s4544[["Vial"]] <- "1024544"
s3812[["Vial"]] <- "1023812"

s4544[["Lot"]] <- "CNS10-GDNF"
s3812[["Lot"]] <- "iNPC-GDNF-WT"

s4544[["Passage"]] <- "p27"
s3812[["Passage"]] <- "p23"

#### Merge all data into a single seurat object ####

tenx1 <- merge(s4544,
  y = c(s3812),
  add.cell.ids = c("4544", "3812"),
  project = project
)

tenx1

rm(s3812, s4544)

## save merged seurat object
saveRDS(tenx1,
  file = paste0("./", date, "_", project, "_merged_seurat_prefilter.rds")
)

# tenx1 <- read_rds(
#   file = paste0("./", date, "_", project, "_merged_seurat_prefilter.rds")
# )

#### Initial QC Filtering ####

data <- tenx1

# Add mito info to metadata
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
data[["percent.ribo"]] <- PercentageFeatureSet(data, pattern = "^RP[SL]")

## Plots of UMI#, Gene#, and %mito
vp1 <- VlnPlot(data, features = "nCount_RNA", pt.size = 0) +
  ggtitle("nCount_RNA prefiltered")
vp1
pdf(paste0("./", date, "_", project, "_prefilter_", "nCount_RNA", ".pdf"))
print(vp1)
dev.off()

vp2 <- VlnPlot(data, features = "nFeature_RNA", pt.size = 0) +
  ggtitle("nFeature_RNA prefiltered")
vp2
pdf(paste0("./", date, "_", project, "_prefilter_", "nFeature_RNA", ".pdf"))
print(vp2)
dev.off()

vp3 <- VlnPlot(data, features = "percent.mt", pt.size = 0) +
  ggtitle("percent.mt prefiltered")
vp3
pdf(paste0("./", date, "_", project, "_prefilter_", "percent.mt", ".pdf"))
print(vp3)
dev.off()

vp4 <- VlnPlot(data, features = "percent.ribo", pt.size = 0) +
  ggtitle("percent.ribo prefiltered")
vp4
pdf(paste0("./", date, "_", project, "_prefilter_", "percent.ribo", ".pdf"))
print(vp4)
dev.off()

p1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt")
p2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p1 + p2

pdf(paste0("./", date, "_", project, "_prefilter_", "scatterplots", ".pdf"))
print(p1 + p2)
dev.off()

# Additional QC and Z-score QC metrics
data[["nUMI.z"]] <- scale(data$nCount_RNA)
data[["nGene.z"]] <- scale(data$nFeature_RNA)
data[["percent.mt.z"]] <- scale(data$percent.mt)
data[["percent.ribo.z"]] <- scale(data$percent.ribo)

## Filter cells based on Z-score
length(data@meta.data$orig.ident)
mean(data@meta.data$percent.mt)
mean(data@meta.data$percent.ribo)
mean(data@meta.data$nCount_RNA)
median(data@meta.data$nCount_RNA)
mean(data@meta.data$nFeature_RNA)
median(data@meta.data$nFeature_RNA)
max(data@meta.data$nCount_RNA)

sink(paste0("./", date, "_", project, "_prefilter_QC_metrics.txt"))
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

sink(paste0("./", date, "_", project, "_postfilter_QC_metrics.txt"))
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

## Plots of UMI#, Gene#, and %mito
vp5 <- VlnPlot(
  data,
  features = "nCount_RNA", pt.size = 0, group.by = "orig.ident"
) +
  NoLegend() + ggtitle("nCount_RNA filtered")
vp5
pdf(paste0("./", date, "_", project, "_filtered_", "nCount_RNA", ".pdf"))
print(vp5)
dev.off()

vp6 <- VlnPlot(
  data,
  features = "nFeature_RNA", pt.size = 0, group.by = "orig.ident"
) +
  NoLegend() + ggtitle("nFeature_RNA filtered")
vp6
pdf(paste0("./", date, "_", project, "_filtered_", "nFeature_RNA", ".pdf"))
print(vp6)
dev.off()

vp7 <- VlnPlot(
  data,
  features = "percent.mt", pt.size = 0, group.by = "orig.ident"
) +
  NoLegend() + ggtitle("percent.mt filtered")
vp7
pdf(paste0("./", date, "_", project, "_filtered_", "percent.mt", ".pdf"))
print(vp7)
dev.off()

vp8 <- VlnPlot(
  data,
  features = "percent.ribo", pt.size = 0, group.by = "orig.ident"
) +
  NoLegend() + ggtitle("percent.ribo filtered")
vp8
pdf(paste0("./", date, "_", project, "_filtered_", "percent.ribo", ".pdf"))
print(vp8)
dev.off()

p3 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt")
p4 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p3 + p4

pdf(paste0("./", date, "_", project, "_filtered_", "scatterplots", ".pdf"))
print(p3 + p4)
dev.off()

# Additional QC metrics
qc.metrics <- as_tibble(
  data[[c(
    "nCount_RNA",
    "nFeature_RNA",
    "percent.mt",
    "percent.ribo"
  )]],
  rownames = "Cell.Barcode"
)

p5 <- qc.metrics %>%
  arrange(percent.mt) %>%
  ggplot(aes(nCount_RNA, nFeature_RNA, colour = percent.mt)) +
  geom_point() +
  scale_color_gradientn(colors = c(
    "black", "blue", "green2", "red", "yellow"
  )) +
  ggtitle("QC metrics")
p5
pdf(paste0("./", date, "_", project, "QC_Metrics", ".pdf"))
print(p5)
dev.off()

p6 <- qc.metrics %>%
  ggplot(aes(percent.mt)) +
  geom_histogram(binwidth = 0.5, fill = "yellow", colour = "black") +
  ggtitle("Distribution of Percentage Mitochondrion")
p6
pdf(paste0("./", date, "_", project, "Dist_mito", ".pdf"))
print(p6)
dev.off()

p7 <- qc.metrics %>%
  ggplot(aes(percent.ribo)) +
  geom_histogram(binwidth = 0.5, fill = "yellow", colour = "black") +
  ggtitle("Distribution of Percentage Ribosomal")
p7
pdf(paste0("./", date, "_", project, "Dist_ribo", ".pdf"))
print(p7)
dev.off()

rm(
  vp1, vp2, vp3, vp4, vp5, vp6, vp7, vp8,
  p1, p2, p3, p4, p5, p6, p7, qc.metrics
)

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

anchors <- FindIntegrationAnchors(
  object.list = s.list,
  reduction = "rpca",
  dims = 1:50
)

s.integrated <- IntegrateData(anchorset = anchors, dims = 1:50)
s.integrated <- ScaleData(s.integrated, verbose = FALSE)

rm(s.list, tenx1, anchors)

#### set factor levels ####

s.integrated$orig.ident <- factor(
  x = s.integrated$orig.ident,
  levels = c("4544", "3812")
)

s.integrated$Lot <- factor(
  x = s.integrated$Lot,
  levels = c("CNS10-GDNF", "iNPC-GDNF-WT")
)

#### save integrated object ####

saveRDS(s.integrated,
  file = paste0(
    "./",
    date, 
    "_",
    project,
    "_post_integration.rds"
  )
)

# s.integrated <- read_rds(
#   file = paste0(
#     "./",
#     date,
#     "_",
#     project,
#     "_post_integration.rds"
#   )
# )

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
# calculate the eigenvalue for each PC in sq.xx.coord by 
# taking the sqrt of the sum of squares
for (i in 1:ncol(sq.xx.coord)) {
  eig[i] <- sqrt(sum(sq.xx.coord[, i]))
}
# calculate the total variance by adding up all the eigenvalues
sum.eig <- sum(eig)
# calculate the expected contribution of all PCs if 
# they all contribute equally to the total variance
expected.contribution <- sum.eig / (length(xx.coord) - 1)
# return the number of principal components with an 
# eigenvalue greater than expected by equal variance
meaningful.PCs <- sum(eig > expected.contribution)

# create empty list for eigenvalue percentage
eig.percent <- c()
# calculate the percentage of the total variance by each PC eigenvalue
for (i in 1:length(eig)) {
  eig.percent[i] <- 100 * eig[i] / sum.eig
}
# sum of all eig.percent should total to 100
sum(eig.percent)
# create empty list for scree values
scree <- c()
# calculate a running total of variance contribution
for (i in 1:length(eig)) {
  if (i == 1) {
    scree[i] <-
      eig.percent[i]
  } else {
    scree[i] <- scree[i - 1] + eig.percent[i]
  }
}

# create data frame for eigenvalue summaries
eigenvalues <- data.frame(
  "PC" = colnames(xx.coord),
  "eig" = eig,
  "percent" = eig.percent,
  "scree" = scree
)

# plot scree values
plot(eigenvalues$percent,
  ylim = c(0, 100), type = "S", xlab = "PC", ylab = "Percent of variance",
  main = paste0(date, "_", project, " scree plot all samples PCA")
)
points(eigenvalues$scree, ylim = c(0, 100), type = "p", pch = 16)
lines(eigenvalues$scree)
# add red line to indicate cut-off
cut.off <- 100 / (length(eig) - 1)
abline(h = cut.off, col = "red")
# add blue line to indicate which PCs are meaningful and kept
abline(v = meaningful.PCs, col = "blue")
text(meaningful.PCs, cut.off,
  label = paste("cutoff PC", meaningful.PCs),
  adj = c(-0.1, -0.5)
)

dev.copy(pdf, paste0("./", date, "_", project, "_scree_plot.pdf"))
dev.off()

# meaningful.PCs <- 13

rm(
  sq.xx.coord, xx.coord, xx.gload, cut.off,
  eig, eig.percent, expected.contribution, features,
  i, scree, sum.eig, eigenvalues
)

#### Run UMAP and look at UMAP plots ####

s.integrated <- RunUMAP(s.integrated,
  reduction = "pca",
  dims = 1:meaningful.PCs,
  verbose = TRUE
)

# update prefixed variable
prefixPC <- paste0("./", date, "_", project, "_", meaningful.PCs, "PCs")

## UMAP plot by sample name ("orig.ident")
DimPlot(s.integrated,
  reduction = "umap",
  label = FALSE,
  pt.size = .25,
  group.by = "orig.ident"
)
hirestiff(paste0(prefixPC, "_UMAP_by_sample_hires.tiff"))
lowrestiff(paste0(prefixPC, "_UMAP_by_sample_lowres.tiff"))

DimPlot(s.integrated,
  reduction = "umap",
  label = FALSE,
  pt.size = .25,
  group.by = "Lot"
)
hirestiff(paste0(prefixPC, "_UMAP_by_lot_hires.tiff"))
lowrestiff(paste0(prefixPC, "_UMAP_by_lot_lowres.tiff"))

DimPlot(s.integrated,
  reduction = "umap",
  label = FALSE,
  pt.size = .25,
  split.by = "Lot"
) +
  NoLegend()
hirestiff(paste0(prefixPC, "_UMAP_by_lot_hires.tiff"))
lowrestiff(paste0(prefixPC, "_UMAP_by_lot_lowres.tiff"))

#### Cell Cycle Score ####

DefaultAssay(s.integrated) <- "RNA"

s.integrated <- CellCycleScoring(s.integrated,
  s.features = cc.genes.updated.2019$s.genes,
  g2m.features = cc.genes.updated.2019$g2m.genes,
  set.ident = TRUE
)

as_tibble(s.integrated[[]]) %>%
  ggplot(aes(Phase)) +
  geom_bar()

DimPlot(s.integrated,
  reduction = "umap",
  label = FALSE,
  pt.size = .25,
  group.by = "Phase"
)
hirestiff(paste0(prefixPC, "_UMAP_by_phase_hires.tiff"))
lowrestiff(paste0(prefixPC, "_UMAP_by_phase_lowres.tiff"))

DimPlot(s.integrated,
  reduction = "umap",
  label = FALSE,
  pt.size = .25,
  split.by = "Lot",
  group.by = "Phase"
)
hirestiff(paste0(prefixPC, "_UMAP_by_phase_splitby_orig.ident_hires.tiff"))
lowrestiff(paste0(prefixPC, "_UMAP_by_phase_splitby_orig.ident_lowres.tiff"))

DefaultAssay(s.integrated) <- "integrated"

saveRDS(s.integrated, file = paste0(prefixPC, "_seurat_integrated.rds"))

# s.integrated <- read_rds(file = paste0(prefixPC,"_seurat_integrated.rds"))

#### Clustering and Resolution ####

# Determine the K-nearest neighbor graph
s.integrated <- FindNeighbors(
  object = s.integrated,
  reduction = "pca",
  dims = 1:meaningful.PCs
)

# Determine the clusters
s.integrated <- FindClusters(
  object = s.integrated,
  resolution = c(0.5)
)

res <- "_res.0.5"
ident_res <- paste0("integrated_snn", res)

Idents(s.integrated) <- ident_res
DimPlot(s.integrated,
  reduction = "umap",
  label = TRUE,
  label.size = 5
) +
  NoLegend() +
  ggtitle(paste0(ident_res))

# update prefix
prefixPCres <- paste0(prefixPC, res)

# Plot the UMAP
DimPlot(s.integrated,
  reduction = "umap",
  label = TRUE,
  label.size = 5
) +
  NoLegend() +
  ggtitle(paste0(ident_res))
hirestiff(paste0(prefixPCres, "_UMAP", "_by_", "cluster", "_hires.tiff"))
lowrestiff(paste0(prefixPCres, "_UMAP", "_by_", "cluster", "_lowres.tiff"))

# UMAP of cells in each cluster by lot
DimPlot(s.integrated,
  reduction = "umap",
  label = FALSE,
  split.by = "Lot"
) +
  NoLegend() +
  ggtitle(paste0(ident_res))
hirestiff(paste0(prefixPCres, "_UMAP", "_by_", "lot", "_hires.tiff"))
lowrestiff(paste0(prefixPCres, "_UMAP", "_by_", "lot", "_lowres.tiff"))

# UMAP of cells in each cluster by Lot with cluster lables
DimPlot(s.integrated,
  reduction = "umap",
  label = TRUE,
  split.by = "Lot"
) +
  NoLegend() +
  ggtitle(paste0(ident_res))
hirestiff(paste0(prefixPCres, "_UMAP", "_by_", "lot", "_hires.tiff"))
lowrestiff(paste0(prefixPCres, "_UMAP", "_by_", "lot", "_lowres.tiff"))

# fix cluster colors
clust_cols <- c(
  "0" = "#F8766D", "1" = "#DB8E00", "2" = "#AEA200",
  "3" = "#64B200", "4" = "#00BD5C", "5" = "#00C1A7",
  "6" = "#00BADE", "7" = "#00A6FF", "8" = "#B385FF",
  "9" = "#EF67EB", "10" = "#FF63B6"
)

# for loop for saving split plots in the same dimensions as composite plot
# using the fixed cluster color scheme

for (i in 1:length(levels(s.integrated$Lot))) {
  p <- (DimPlot(s.integrated,
    reduction = "umap",
    label = FALSE,
    pt.size = .25,
    cells = c(WhichCells(s.integrated,
      expression = Lot == levels(s.integrated$Lot)[i]
    )),
    cols = clust_cols
  ) +
    NoLegend() +
    ggtitle(paste0(levels(s.integrated$Lot)[i])))

  print(p)

  hirestiff(paste(prefixPCres, "UMAP", "split_by", "lot", "no_cluster_labels",
    levels(s.integrated$Lot)[i], "hires.tiff",
    sep = "_"
  ))
  lowrestiff(paste(prefixPCres, "UMAP", "split_by", "lot", "no_cluster_labels",
    levels(s.integrated$Lot)[i], "lowres.tiff",
    sep = "_"
  ))
}

# Extract number of cells per cluster per Lot
n_cells <- FetchData(s.integrated, vars = c("ident", "Lot")) %>%
  dplyr::count(ident, Lot) %>%
  tidyr::spread(ident, n)

write.csv(n_cells, file = paste0(prefixPCres, "_cells_per_cluster.csv"))

rm(n_cells, p, i)

# save RDS containing reduction and cluster idents
saveRDS(s.integrated, paste0(prefixPCres, "_seurat_after_clustering.rds"))

# s.integrated <- readRDS(
#   file = paste0(prefixPCres, "_seurat_after_clustering.rds")
# )

#### Normalize RNA slot for visualization ####

# Select the RNA counts slot to be the default assay for visualization purposes
DefaultAssay(s.integrated) <- "RNA"

# Normalize, find variable features, scale data
s.integrated <- NormalizeData(s.integrated, verbose = FALSE)
s.integrated <- FindVariableFeatures(s.integrated)
all.genes <- rownames(s.integrated)
s.integrated <- ScaleData(s.integrated, features = all.genes)

# save RDS containing RNA normalized data
saveRDS(s.integrated, paste0(prefixPCres, "_seurat_after_RNAnorm.rds"))

# s.integrated <- readRDS(
#   file = paste0(prefixPCres, "_seurat_after_RNAnorm.rds")
# )

#### Figure 1b ####

# get UMAP coordinates
umap_coords <- s.integrated@reductions[["umap"]]@cell.embeddings

# calculate minimum and maximum values of UMAP1 and UMAP2
# These can be used to fix the axis ranges on subsequent plots using:
# + coord_cartesian(xlim = c(umap1_range),ylim = c(umap2_range))
umap1_range <- range(umap_coords[, 1])
umap2_range <- range(umap_coords[, 2])

# Plot the UMAP
DimPlot(s.integrated,
  reduction = "umap",
  label = FALSE,
  pt.size = 2,
  label.size = 5
) +
  NoLegend() +
  coord_cartesian(
    xlim = c(umap1_range),
    ylim = c(umap2_range)
  )

hireseps(paste0(prefixPCres, "_Fig", "_1b", "_newfunction", "_hires.eps"))

#### Figure 1c ####
# for loop that saves UMAP of each Lot by cluster without cluster labels with 
# the same dimensions as composite plot using the fixed cluster color scheme

# fix cluster colors
clust_cols <- c(
  "0" = "#F8766D", "1" = "#DB8E00", "2" = "#AEA200",
  "3" = "#64B200", "4" = "#00BD5C", "5" = "#00C1A7",
  "6" = "#00BADE", "7" = "#00A6FF", "8" = "#B385FF",
  "9" = "#EF67EB", "10" = "#FF63B6"
)

# for loop
for (i in 1:length(levels(s.integrated$Lot))) {
  p <- (DimPlot(s.integrated,
    reduction = "umap",
    label = FALSE,
    pt.size = 2,
    cells = c(WhichCells(
      s.integrated,
      expression = Lot == levels(s.integrated$Lot)[i]
    )),
    cols = clust_cols
  ) +
    NoLegend() +
    coord_cartesian(
      xlim = c(umap1_range),
      ylim = c(umap2_range)
    ))

  print(p)

  hireseps(paste(prefixPCres, "Fig", "1c", "fixed_axis", "no_title", "new",
    levels(s.integrated$Lot)[i], "hires.eps",
    sep = "_"
  ))
}

#### Figure 1d ####

FeaturePlot(s.integrated,
  reduction = "umap",
  features = c("S100B"),
  order = TRUE,
  min.cutoff = "q10",
  label = FALSE,
  label.size = 4,
  pt.size = 2
)

hireseps(paste(prefixPCres, "Fig", "1d", "S100B", "hires.eps", sep = "_"))

FeaturePlot(s.integrated,
  reduction = "umap",
  features = c("VIM"),
  order = TRUE,
  min.cutoff = "q10",
  label = FALSE,
  label.size = 4,
  pt.size = 2
)

hireseps(paste(prefixPCres, "Fig", "1d", "VIM", "hires.eps", sep = "_"))

FeaturePlot(s.integrated,
  reduction = "umap",
  features = c("CD99"),
  order = TRUE,
  min.cutoff = "q10",
  label = FALSE,
  label.size = 4,
  pt.size = 2
)

hireseps(paste(prefixPCres, "Fig", "1d", "CD99", "hires.eps", sep = "_"))

FeaturePlot(s.integrated,
  reduction = "umap",
  features = c("NES"),
  order = TRUE,
  min.cutoff = "q10",
  label = FALSE,
  label.size = 4,
  pt.size = 2
)

hireseps(paste(prefixPCres, "Fig", "1d", "NES", "hires.eps", sep = "_"))

#### Figure 1e ####

FeaturePlot(s.integrated,
  reduction = "umap",
  features = c("GFAP"),
  order = TRUE,
  min.cutoff = "q10",
  label = FALSE,
  label.size = 4,
  pt.size = 2
)

hireseps(paste(prefixPCres, "Fig", "1e", "GFAP", "hires.eps", sep = "_"))

FeaturePlot(s.integrated,
  reduction = "umap",
  features = c("AQP4"),
  order = TRUE,
  min.cutoff = "q10",
  label = FALSE,
  label.size = 4,
  pt.size = 2
)

hireseps(paste(prefixPCres, "Fig", "1e", "AQP4", "hires.eps", sep = "_"))

FeaturePlot(s.integrated,
  reduction = "umap",
  features = c("SLC1A3"),
  order = TRUE,
  min.cutoff = "q10",
  label = FALSE,
  label.size = 4,
  pt.size = 2
)

hireseps(paste(prefixPCres, "Fig", "1e", "SLC1A3", "hires.eps", sep = "_"))



#### Figure S1b ####

FeaturePlot(s.integrated,
  reduction = "umap",
  features = c("POU5F1"),
  order = TRUE,
  min.cutoff = "q10",
  label = FALSE,
  label.size = 4,
  pt.size = 2
)

hireseps(paste(prefixPCres, "Fig", "S1b", "POU5F1", "hires.eps", sep = "_"))

FeaturePlot(s.integrated,
  reduction = "umap",
  features = c("NANOG"),
  order = TRUE,
  min.cutoff = "q10",
  label = FALSE,
  label.size = 4,
  pt.size = 2
)

hireseps(paste(prefixPCres, "Fig", "S1b", "NANOG", "hires.eps", sep = "_"))

FeaturePlot(s.integrated,
  reduction = "umap",
  features = c("SOX2"),
  order = TRUE,
  min.cutoff = "q10",
  label = FALSE,
  label.size = 4,
  pt.size = 2
)

hireseps(paste(prefixPCres, "Fig", "S1b", "SOX2", "hires.eps", sep = "_"))

FeaturePlot(s.integrated,
  reduction = "umap",
  features = c("KLF4"),
  order = TRUE,
  min.cutoff = "q10",
  label = FALSE,
  label.size = 4,
  pt.size = 2
)

hireseps(paste(prefixPCres, "Fig", "S1b", "KLF4", "hires.eps", sep = "_"))

#### Figure S1c ####

expressiondata <- FetchData(
  object = s.integrated,
  vars = c(
    "orig.ident", "Lot",
    "POU5F1", "NANOG", "KLF4"
  )
)

expressiondata <- rownames_to_column(expressiondata, var = "cell_id")

cells_pouf1 <- expressiondata$cell_id[expressiondata$POU5F1 > 0]
cells_nanog <- expressiondata$cell_id[expressiondata$NANOG > 0]
cells_klf4 <- expressiondata$cell_id[expressiondata$KLF4 > 0]

genelists <- list(cells_pouf1, cells_nanog, cells_klf4)
names(genelists) <- c("OCT4", "NANOG", "KLF4")

venngenes <- Venn(genelists)

plot(venngenes, doWeights = FALSE, type = "circles")

pdf(paste(prefixPCres, "Fig", "S1c.pdf", sep = "_"))
plot(venngenes, doWeights = FALSE, type = "circles")
dev.off()

postscript(paste(prefixPCres, "Fig", "S1c.eps", sep = "_"))
plot(venngenes, doWeights = FALSE, type = "circles")
dev.off()

#### Figure S1d ####

DimPlot(s.integrated,
  reduction = "umap",
  label = FALSE,
  pt.size = 2,
  group.by = "Phase"
)

hireseps(paste0(prefixPCres, "_Fig_S1d_hires.eps"))

#### write metadata ####

write.csv(
  s.integrated@meta.data,
  file = paste(prefixPCres,
    "metadata.csv",
    sep = "_"
  )
)

#### session info ####
 
sessionInfo()
# R version 4.2.2 (2022-10-31)
# Platform: aarch64-apple-darwin20 (64-bit)
# Running under: macOS Ventura 13.2.1
# 
# Matrix products: default
# LAPACK: /Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] Vennerable_3.1.0.9000 gridExtra_2.3         irlba_2.3.5.1         patchwork_1.1.2       cowplot_1.1.1         Matrix_1.5-3         
#   [7] lubridate_1.9.2       forcats_1.0.0         stringr_1.5.0         dplyr_1.1.1           purrr_1.0.1           readr_2.1.4          
#  [13] tidyr_1.3.0           tibble_3.2.1          ggplot2_3.4.1         tidyverse_2.0.0       SeuratObject_4.1.3    Seurat_4.3.0         
# 
# loaded via a namespace (and not attached):
#   [1] utf8_1.2.3             R.utils_2.12.2         spatstat.explore_3.1-0 reticulate_1.28        tidyselect_1.2.0      
#   [6] htmlwidgets_1.6.2      grid_4.2.2             Rtsne_0.16             pROC_1.18.0            devtools_2.4.5        
#  [11] munsell_0.5.0          ragg_1.2.5             codetools_0.2-19       ica_1.0-3              future_1.32.0         
#  [16] miniUI_0.1.1.1         withr_2.5.0            spatstat.random_3.1-4  colorspace_2.1-0       progressr_0.13.0      
#  [21] knitr_1.42             rstudioapi_0.14        stats4_4.2.2           ROCR_1.0-11            tensor_1.5            
#  [26] listenv_0.9.0          labeling_0.4.2         polyclip_1.10-4        farver_2.1.1           parallelly_1.35.0     
#  [31] vctrs_0.6.1            generics_0.1.3         ipred_0.9-14           xfun_0.38              timechange_0.2.0      
#  [36] randomForest_4.7-1.1   R6_2.5.1               doParallel_1.0.17      clue_0.3-64            bitops_1.0-7          
#  [41] spatstat.utils_3.0-2   cachem_1.0.7           promises_1.2.0.1       scales_1.2.1           nnet_7.3-18           
#  [46] gtable_0.3.3           globals_0.16.2         processx_3.8.0         goftest_1.2-3          timeDate_4022.108     
#  [51] rlang_1.1.0            systemfonts_1.0.4      GlobalOptions_0.1.2    splines_4.2.2          lazyeval_0.2.2        
#  [56] ModelMetrics_1.2.2.2   spatstat.geom_3.1-0    checkmate_2.1.0        reshape2_1.4.4         abind_1.4-5           
#  [61] backports_1.4.1        httpuv_1.6.9           Hmisc_5.0-1            RBGL_1.74.0            caret_6.0-94          
#  [66] DiagrammeR_1.0.9       tools_4.2.2            lava_1.7.2.1           usethis_2.1.6          ellipsis_0.3.2        
#  [71] RColorBrewer_1.1-3     proxy_0.4-27           BiocGenerics_0.44.0    sessioninfo_1.2.2      ggridges_0.5.4        
#  [76] Rcpp_1.0.10            plyr_1.8.8             base64enc_0.1-3        visNetwork_2.1.2       ps_1.7.3              
#  [81] prettyunits_1.1.1      rpart_4.1.19           deldir_1.0-6           pbapply_1.7-0          GetoptLong_1.0.5      
#  [86] urlchecker_1.0.1       S4Vectors_0.36.2       zoo_1.8-11             nichenetr_1.1.1        ggrepel_0.9.3         
#  [91] cluster_2.1.4          fs_1.6.1               magrittr_2.0.3         data.table_1.14.8      scattermore_0.8       
#  [96] circlize_0.4.15        lmtest_0.9-40          RANN_2.6.1             R.cache_0.16.0         fitdistrplus_1.1-8    
# [101] matrixStats_0.63.0     pkgload_1.3.2          hms_1.1.3              mime_0.12              evaluate_0.20         
# [106] xtable_1.8-4           IRanges_2.32.0         shape_1.4.6            compiler_4.2.2         KernSmooth_2.23-20    
# [111] crayon_1.5.2           R.oo_1.25.0            htmltools_0.5.5        later_1.3.0            tzdb_0.3.0            
# [116] Formula_1.2-5          DBI_1.1.3              ComplexHeatmap_2.14.0  MASS_7.3-58.3          cli_3.6.1             
# [121] R.methodsS3_1.8.2      parallel_4.2.2         gower_1.0.1            igraph_1.4.1           pkgconfig_2.0.3       
# [126] foreign_0.8-84         sp_1.6-0               plotly_4.10.1          spatstat.sparse_3.0-1  recipes_1.0.5         
# [131] foreach_1.5.2          hardhat_1.2.0          prodlim_2019.11.13     callr_3.7.3            digest_0.6.31         
# [136] sctransform_0.3.5      RcppAnnoy_0.0.20       graph_1.76.0           spatstat.data_3.0-1    rmarkdown_2.21        
# [141] leiden_0.4.3           htmlTable_2.4.1        uwot_0.1.14            shiny_1.7.4            rjson_0.2.21          
# [146] lifecycle_1.0.3        nlme_3.1-162           jsonlite_1.8.4         limma_3.54.2           viridisLite_0.4.1     
# [151] fansi_1.0.4            pillar_1.9.0           lattice_0.20-45        fastmap_1.1.1          httr_1.4.5            
# [156] pkgbuild_1.4.0         survival_3.5-5         glue_1.6.2             remotes_2.4.2          fdrtool_1.2.17        
# [161] png_0.1-8              iterators_1.0.14       class_7.3-21           stringi_1.7.12         profvis_0.3.7         
# [166] textshaping_0.3.6      caTools_1.18.2         memoise_2.0.1          styler_1.9.1           e1071_1.7-13          
# [171] future.apply_1.10.0 

#### devtools session info ####

devtools::session_info()
# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.2.2 (2022-10-31)
# os       macOS Ventura 13.2.1
# system   aarch64, darwin20
# ui       RStudio
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       America/Los_Angeles
# date     2023-03-28
# rstudio  2023.03.0+386 Cherry Blossom (desktop)
# pandoc   2.19.2 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/ (via rmarkdown)
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# package          * version    date (UTC) lib source
# abind              1.4-5      2016-07-21 [1] CRAN (R 4.2.0)
# backports          1.4.1      2021-12-13 [1] CRAN (R 4.2.0)
# base64enc          0.1-3      2015-07-28 [1] CRAN (R 4.2.0)
# BiocGenerics       0.44.0     2022-11-07 [1] Bioconductor
# bitops             1.0-7      2021-04-24 [1] CRAN (R 4.2.0)
# cachem             1.0.7      2023-02-24 [1] CRAN (R 4.2.0)
# callr              3.7.3      2022-11-02 [1] CRAN (R 4.2.0)
# caret              6.0-94     2023-03-21 [1] CRAN (R 4.2.0)
# caTools            1.18.2     2021-03-28 [1] CRAN (R 4.2.0)
# checkmate          2.1.0      2022-04-21 [1] CRAN (R 4.2.0)
# circlize           0.4.15     2022-05-10 [1] CRAN (R 4.2.0)
# class              7.3-21     2023-01-23 [1] CRAN (R 4.2.0)
# cli                3.6.1      2023-03-23 [1] CRAN (R 4.2.0)
# clue               0.3-64     2023-01-31 [1] CRAN (R 4.2.0)
# cluster            2.1.4      2022-08-22 [1] CRAN (R 4.2.2)
# codetools          0.2-19     2023-02-01 [1] CRAN (R 4.2.0)
# colorspace         2.1-0      2023-01-23 [1] CRAN (R 4.2.0)
# ComplexHeatmap     2.14.0     2022-11-07 [1] Bioconductor
# cowplot          * 1.1.1      2020-12-30 [1] CRAN (R 4.2.0)
# crayon             1.5.2      2022-09-29 [1] CRAN (R 4.2.0)
# data.table         1.14.8     2023-02-17 [1] CRAN (R 4.2.0)
# DBI                1.1.3      2022-06-18 [1] CRAN (R 4.2.0)
# deldir             1.0-6      2021-10-23 [1] CRAN (R 4.2.0)
# devtools           2.4.5      2022-10-11 [1] CRAN (R 4.2.0)
# DiagrammeR         1.0.9      2022-03-05 [1] CRAN (R 4.2.0)
# digest             0.6.31     2022-12-11 [1] CRAN (R 4.2.0)
# doParallel         1.0.17     2022-02-07 [1] CRAN (R 4.2.0)
# dplyr            * 1.1.1      2023-03-22 [1] CRAN (R 4.2.0)
# e1071              1.7-13     2023-02-01 [1] CRAN (R 4.2.0)
# ellipsis           0.3.2      2021-04-29 [1] CRAN (R 4.2.0)
# evaluate           0.20       2023-01-17 [1] CRAN (R 4.2.0)
# fansi              1.0.4      2023-01-22 [1] CRAN (R 4.2.0)
# farver             2.1.1      2022-07-06 [1] CRAN (R 4.2.0)
# fastmap            1.1.1      2023-02-24 [1] CRAN (R 4.2.0)
# fdrtool            1.2.17     2021-11-13 [1] CRAN (R 4.2.0)
# fitdistrplus       1.1-8      2022-03-10 [1] CRAN (R 4.2.0)
# forcats          * 1.0.0      2023-01-29 [1] CRAN (R 4.2.0)
# foreach            1.5.2      2022-02-02 [1] CRAN (R 4.2.0)
# foreign            0.8-84     2022-12-06 [1] CRAN (R 4.2.0)
# Formula            1.2-5      2023-02-24 [1] CRAN (R 4.2.0)
# fs                 1.6.1      2023-02-06 [1] CRAN (R 4.2.0)
# future             1.32.0     2023-03-07 [1] CRAN (R 4.2.0)
# future.apply       1.10.0     2022-11-05 [1] CRAN (R 4.2.0)
# generics           0.1.3      2022-07-05 [1] CRAN (R 4.2.0)
# GetoptLong         1.0.5      2020-12-15 [1] CRAN (R 4.2.0)
# ggplot2          * 3.4.1      2023-02-10 [1] CRAN (R 4.2.0)
# ggrepel            0.9.3      2023-02-03 [1] CRAN (R 4.2.0)
# ggridges           0.5.4      2022-09-26 [1] CRAN (R 4.2.0)
# GlobalOptions      0.1.2      2020-06-10 [1] CRAN (R 4.2.0)
# globals            0.16.2     2022-11-21 [1] CRAN (R 4.2.0)
# glue               1.6.2      2022-02-24 [1] CRAN (R 4.2.0)
# goftest            1.2-3      2021-10-07 [1] CRAN (R 4.2.0)
# gower              1.0.1      2022-12-22 [1] CRAN (R 4.2.0)
# graph              1.76.0     2022-11-07 [1] Bioconductor
# gridExtra        * 2.3        2017-09-09 [1] CRAN (R 4.2.0)
# gtable             0.3.3      2023-03-21 [1] CRAN (R 4.2.0)
# hardhat            1.2.0      2022-06-30 [1] CRAN (R 4.2.0)
# Hmisc              5.0-1      2023-03-08 [1] CRAN (R 4.2.0)
# hms                1.1.3      2023-03-21 [1] CRAN (R 4.2.0)
# htmlTable          2.4.1      2022-07-07 [1] CRAN (R 4.2.0)
# htmltools          0.5.5      2023-03-23 [1] CRAN (R 4.2.0)
# htmlwidgets        1.6.2      2023-03-17 [1] CRAN (R 4.2.0)
# httpuv             1.6.9      2023-02-14 [1] CRAN (R 4.2.0)
# httr               1.4.5      2023-02-24 [1] CRAN (R 4.2.0)
# ica                1.0-3      2022-07-08 [1] CRAN (R 4.2.0)
# igraph             1.4.1      2023-02-24 [1] CRAN (R 4.2.0)
# ipred              0.9-14     2023-03-09 [1] CRAN (R 4.2.0)
# IRanges            2.32.0     2022-11-07 [1] Bioconductor
# irlba            * 2.3.5.1    2022-10-03 [1] CRAN (R 4.2.0)
# iterators          1.0.14     2022-02-05 [1] CRAN (R 4.2.0)
# jsonlite           1.8.4      2022-12-06 [1] CRAN (R 4.2.0)
# KernSmooth         2.23-20    2021-05-03 [1] CRAN (R 4.2.2)
# knitr              1.42       2023-01-25 [1] CRAN (R 4.2.0)
# labeling           0.4.2      2020-10-20 [1] CRAN (R 4.2.0)
# later              1.3.0      2021-08-18 [1] CRAN (R 4.2.0)
# lattice            0.20-45    2021-09-22 [1] CRAN (R 4.2.2)
# lava               1.7.2.1    2023-02-27 [1] CRAN (R 4.2.0)
# lazyeval           0.2.2      2019-03-15 [1] CRAN (R 4.2.0)
# leiden             0.4.3      2022-09-10 [1] CRAN (R 4.2.0)
# lifecycle          1.0.3      2022-10-07 [1] CRAN (R 4.2.0)
# limma              3.54.2     2023-03-01 [1] Bioconductor
# listenv            0.9.0      2022-12-16 [1] CRAN (R 4.2.0)
# lmtest             0.9-40     2022-03-21 [1] CRAN (R 4.2.0)
# lubridate        * 1.9.2      2023-02-10 [1] CRAN (R 4.2.0)
# magrittr           2.0.3      2022-03-30 [1] CRAN (R 4.2.0)
# MASS               7.3-58.3   2023-03-07 [1] CRAN (R 4.2.0)
# Matrix           * 1.5-3      2022-11-11 [1] CRAN (R 4.2.0)
# matrixStats        0.63.0     2022-11-18 [1] CRAN (R 4.2.0)
# memoise            2.0.1      2021-11-26 [1] CRAN (R 4.2.0)
# mime               0.12       2021-09-28 [1] CRAN (R 4.2.0)
# miniUI             0.1.1.1    2018-05-18 [1] CRAN (R 4.2.0)
# ModelMetrics       1.2.2.2    2020-03-17 [1] CRAN (R 4.2.0)
# munsell            0.5.0      2018-06-12 [1] CRAN (R 4.2.0)
# nichenetr          1.1.1      2023-03-16 [1] Github (saeyslab/nichenetr@cc3bced)
# nlme               3.1-162    2023-01-31 [1] CRAN (R 4.2.0)
# nnet               7.3-18     2022-09-28 [1] CRAN (R 4.2.2)
# parallelly         1.35.0     2023-03-23 [1] CRAN (R 4.2.0)
# patchwork        * 1.1.2      2022-08-19 [1] CRAN (R 4.2.0)
# pbapply            1.7-0      2023-01-13 [1] CRAN (R 4.2.0)
# pillar             1.9.0      2023-03-22 [1] CRAN (R 4.2.0)
# pkgbuild           1.4.0      2022-11-27 [1] CRAN (R 4.2.0)
# pkgconfig          2.0.3      2019-09-22 [1] CRAN (R 4.2.0)
# pkgload            1.3.2      2022-11-16 [1] CRAN (R 4.2.0)
# plotly             4.10.1     2022-11-07 [1] CRAN (R 4.2.0)
# plyr               1.8.8      2022-11-11 [1] CRAN (R 4.2.0)
# png                0.1-8      2022-11-29 [1] CRAN (R 4.2.0)
# polyclip           1.10-4     2022-10-20 [1] CRAN (R 4.2.0)
# prettyunits        1.1.1      2020-01-24 [1] CRAN (R 4.2.0)
# pROC               1.18.0     2021-09-03 [1] CRAN (R 4.2.0)
# processx           3.8.0      2022-10-26 [1] CRAN (R 4.2.0)
# prodlim            2019.11.13 2019-11-17 [1] CRAN (R 4.2.0)
# profvis            0.3.7      2020-11-02 [1] CRAN (R 4.2.0)
# progressr          0.13.0     2023-01-10 [1] CRAN (R 4.2.0)
# promises           1.2.0.1    2021-02-11 [1] CRAN (R 4.2.0)
# proxy              0.4-27     2022-06-09 [1] CRAN (R 4.2.0)
# ps                 1.7.3      2023-03-21 [1] CRAN (R 4.2.0)
# purrr            * 1.0.1      2023-01-10 [1] CRAN (R 4.2.0)
# R.cache            0.16.0     2022-07-21 [1] CRAN (R 4.2.0)
# R.methodsS3        1.8.2      2022-06-13 [1] CRAN (R 4.2.0)
# R.oo               1.25.0     2022-06-12 [1] CRAN (R 4.2.0)
# R.utils            2.12.2     2022-11-11 [1] CRAN (R 4.2.0)
# R6                 2.5.1      2021-08-19 [1] CRAN (R 4.2.0)
# ragg               1.2.5      2023-01-12 [1] CRAN (R 4.2.0)
# randomForest       4.7-1.1    2022-05-23 [1] CRAN (R 4.2.0)
# RANN               2.6.1      2019-01-08 [1] CRAN (R 4.2.0)
# RBGL               1.74.0     2022-11-07 [1] Bioconductor
# RColorBrewer       1.1-3      2022-04-03 [1] CRAN (R 4.2.0)
# Rcpp               1.0.10     2023-01-22 [1] CRAN (R 4.2.0)
# RcppAnnoy          0.0.20     2022-10-27 [1] CRAN (R 4.2.0)
# readr            * 2.1.4      2023-02-10 [1] CRAN (R 4.2.0)
# recipes            1.0.5      2023-02-20 [1] CRAN (R 4.2.0)
# remotes            2.4.2      2021-11-30 [1] CRAN (R 4.2.0)
# reshape2           1.4.4      2020-04-09 [1] CRAN (R 4.2.0)
# reticulate         1.28       2023-01-27 [1] CRAN (R 4.2.0)
# rjson              0.2.21     2022-01-09 [1] CRAN (R 4.2.0)
# rlang              1.1.0      2023-03-14 [1] CRAN (R 4.2.0)
# rmarkdown          2.21       2023-03-26 [1] CRAN (R 4.2.2)
# ROCR               1.0-11     2020-05-02 [1] CRAN (R 4.2.0)
# rpart              4.1.19     2022-10-21 [1] CRAN (R 4.2.2)
# rstudioapi         0.14       2022-08-22 [1] CRAN (R 4.2.0)
# Rtsne              0.16       2022-04-17 [1] CRAN (R 4.2.0)
# S4Vectors          0.36.2     2023-03-01 [1] Bioconductor
# scales             1.2.1      2022-08-20 [1] CRAN (R 4.2.0)
# scattermore        0.8        2022-02-14 [1] CRAN (R 4.2.0)
# sctransform        0.3.5      2022-09-21 [1] CRAN (R 4.2.0)
# sessioninfo        1.2.2      2021-12-06 [1] CRAN (R 4.2.0)
# Seurat           * 4.3.0      2022-11-18 [1] CRAN (R 4.2.0)
# SeuratObject     * 4.1.3      2022-11-07 [1] CRAN (R 4.2.0)
# shape              1.4.6      2021-05-19 [1] CRAN (R 4.2.0)
# shiny              1.7.4      2022-12-15 [1] CRAN (R 4.2.0)
# sp                 1.6-0      2023-01-19 [1] CRAN (R 4.2.0)
# spatstat.data      3.0-1      2023-03-12 [1] CRAN (R 4.2.0)
# spatstat.explore   3.1-0      2023-03-14 [1] CRAN (R 4.2.0)
# spatstat.geom      3.1-0      2023-03-12 [1] CRAN (R 4.2.0)
# spatstat.random    3.1-4      2023-03-13 [1] CRAN (R 4.2.0)
# spatstat.sparse    3.0-1      2023-03-12 [1] CRAN (R 4.2.0)
# spatstat.utils     3.0-2      2023-03-11 [1] CRAN (R 4.2.0)
# stringi            1.7.12     2023-01-11 [1] CRAN (R 4.2.0)
# stringr          * 1.5.0      2022-12-02 [1] CRAN (R 4.2.0)
# styler             1.9.1      2023-03-04 [1] CRAN (R 4.2.0)
# survival           3.5-5      2023-03-12 [1] CRAN (R 4.2.0)
# systemfonts        1.0.4      2022-02-11 [1] CRAN (R 4.2.0)
# tensor             1.5        2012-05-05 [1] CRAN (R 4.2.0)
# textshaping        0.3.6      2021-10-13 [1] CRAN (R 4.2.0)
# tibble           * 3.2.1      2023-03-20 [1] CRAN (R 4.2.0)
# tidyr            * 1.3.0      2023-01-24 [1] CRAN (R 4.2.0)
# tidyselect         1.2.0      2022-10-10 [1] CRAN (R 4.2.0)
# tidyverse        * 2.0.0      2023-02-22 [1] CRAN (R 4.2.0)
# timechange         0.2.0      2023-01-11 [1] CRAN (R 4.2.0)
# timeDate           4022.108   2023-01-07 [1] CRAN (R 4.2.0)
# tzdb               0.3.0      2022-03-28 [1] CRAN (R 4.2.0)
# urlchecker         1.0.1      2021-11-30 [1] CRAN (R 4.2.0)
# usethis            2.1.6      2022-05-25 [1] CRAN (R 4.2.0)
# utf8               1.2.3      2023-01-31 [1] CRAN (R 4.2.0)
# uwot               0.1.14     2022-08-22 [1] CRAN (R 4.2.0)
# vctrs              0.6.1      2023-03-22 [1] CRAN (R 4.2.0)
# Vennerable       * 3.1.0.9000 2023-03-16 [1] Github (js229/Vennerable@46057c9)
# viridisLite        0.4.1      2022-08-22 [1] CRAN (R 4.2.0)
# visNetwork         2.1.2      2022-09-29 [1] CRAN (R 4.2.0)
# withr              2.5.0      2022-03-03 [1] CRAN (R 4.2.0)
# xfun               0.38       2023-03-24 [1] CRAN (R 4.2.0)
# xtable             1.8-4      2019-04-21 [1] CRAN (R 4.2.0)
# zoo                1.8-11     2022-09-17 [1] CRAN (R 4.2.0)
# 
# [1] /Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/library
# 
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
