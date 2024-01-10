# Example codes for calculating connectivity signature score 
# Using AddModuleScore function in Seurat

library(Seurat)
SeuratObject <- readRDS('YourSeuratObject.rds')

# Connectivity signature (71 DEGs derived from SR101 scRNA-Seq data) 
Signature <- read.csv("ConnectivitySignature.csv",stringsAsFactors = F)

# Connectivity signature score
upG <- as.character(Signature[Signature$GeneGroup == 'Upregulated','Gene'])
downG <- as.character(Signature[Signature$GeneGroup == 'Downregulated','Gene'])
SeuratObject <- AddModuleScore(SeuratObject,list(upG),name='up')
SeuratObject <- AddModuleScore(SeuratObject,list(downG),name='down')
SeuratObject$ScoreWithoutScale <- SeuratObject$up1-SeuratObject$down1

# Scale and center scores across cells and winsorized to -3 and 3
SeuratObject$ConnectivityScore <- scale(SeuratObject$ScoreWithoutScale)
SeuratObject$ConnectivityScore[SeuratObject$ConnectivityScore > 3] <- 3
SeuratObject$ConnectivityScore[SeuratObject$ConnectivityScore < -3] <- -3
