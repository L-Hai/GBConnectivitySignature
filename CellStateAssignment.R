# Example codes for cell state assignment based on cell state signature score
# Using AddModuleScore function in Seurat

library(Seurat)
SeuratObject <- readRDS('YourSeuratObject.rds')

# Cell state signatures obtained from Neftel et. al. 2019
CellStateSignatures <- read.csv('CellStateSignatures.csv',stringsAsFactor=F,na='')
CellStateSignatures <- as.list(CellStateSignatures)
CellStateSignatures <- lapply(CellStateSignatures, function(x) x[!is.na(x)])

DefaultAssay(SeuratObject) <- "RNA"
SeuratObject <- AddModuleScore(SeuratObject, CellStateSignatures, name = 'CellStateSignatures')

# Assign cell state to each cell based on the highest cell state signature score
MetaIndex <- grep('CellStateSignatures', colnames(SeuratObject@meta.data))
colnames(SeuratObject@meta.data)[MetaIndex] <- names(CellStateSignatures)
meta <- SeuratObject@meta.data
CellState.Meta <- meta[,MetaIndex]
CellState.Name <- colnames(CellState.Meta)
CellState.Result <- c()
for(i in 1:nrow(CellState.Meta)){
  CellState.Result <- c(CellState.Result,CellState.Name[which.max(CellState.Meta[i,])])
}
SeuratObject$CellState <- CellState.Result