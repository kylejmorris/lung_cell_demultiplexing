# load the PBMC dataset from filtered_gene_bc_matrices/hg19

library(dplyr)
library(Seurat)
library(patchwork)
library(readxl)

# Load the PBMC dataset
# PBMC = peripheral blood mononuclear cells (immune cells)
pbmc.data <- Read10X(data.dir = "FunWithKyle/1013_intron/filtered_feature_bc_matrix/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "s1013", min.cells = 3, min.features = 200)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-") # mitochondrial genes are identified as ^MT- in results matrix

# Gene == a piece of RNA that codes for a protein
# Feature_RNA = # genes detected in cell
# Count_RNA = # reads detected in cell

# cutoff any cell in percent.mt above 25%  
pbmc <- subset(pbmc, subset = percent.mt < 25)

# cutoff Feature_RNA <6000  and > 200
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 6000)

# Visualize QC metrics as a violin plot
# Each dot in plot is a cell, and then shows 
p <- VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


# normalize the data
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes across cells
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

p2 <- VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

p3 <- plot(DimPlot(pbmc, reduction = "pca"))

pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)

p4 <- plot(JackStrawPlot(pbmc, dims = 1:20))
#plot(p4)

p5 <- plot(ElbowPlot(pbmc))
plot(p5)


# add 3 columns to object metadata for outcomes of 3 demultiplexing methods. Map the barcode cell from each file with outcome
# append empty column to object metadata in pbmc
#pbmc[['Freemuxlet_SNP']] <- NA
#pbmc[['Souporcell']] <- NA
#pbmc[['Souporcell_SNP']] <- NA

#sort rows of ppmc by barcode
#pbmc <- pbmc[order(pbmc$barcode),]

# read ThressDemulti_S1013 excel into a file object
#obj <- read_excel("FunWithKyle/Demultiplexing_output/ThreeDemulti_S1013.xlsx", col_names = c("barcodes", "Freemuxlet_SNP", "Souporcell", "Souporcell_SNP")) 
#obj.head <- head(obj)

# rename obj column barcode to barcodes
#colFreemuxlet <- obj$Freemuxlet_SNP
#colSouporcell <- obj$Souporcell
#colSouporcell_SNP <- obj$Souporcell_SNP

# sort pbmc data ascending based on barcodes

# for each row in seurat object get barcode 
#print(pbmc@meta.data[,2])
#for (barcode in rownames(pbmc@meta.data)) {
#  # get barcode from row
#  print(barcode)
#
#  # get row from obj where barcode matches
#  row <- obj[obj$barcodes == barcode,]
#
#  print(row)
#
#  ## get value from row
#  a <- row$Freemuxlet_SNP
#  b <- row$Souporcell
#  c <- row$Souporcell_SNP
#  
#  print("Column values of obj with matching barcode")
#  print(typeof(a))
#  print(b)
#  print(c)
#
#  # set value to row in pbmc
#  pbmc.meta.data[i,1] <- "test"
#  #pbmc@meta.data[[i]] <- c(a)
#  #pbmc@meta.data.SetValue(i, "Souporcell", b)
#  #pbmc@meta.data.SetValue(i, "Souporcell_SNP",c)
#}

#print("header after append")
#print(head(pbmc@meta.data, 1000))