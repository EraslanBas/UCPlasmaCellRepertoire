source("Main.R")
reticulate::use_python("/home/beraslan/anaconda3/bin/python", required = T)

#library("future")
#library(future.apply)

#plan("multiprocess", workers = 12)
#options(future.globals.maxSize=+Inf)


samps <-  readRDS(paste0(dataDir, "/RDSFiles/allSamplesReplicatesMerged.rds"))


for(i in 1:length(samps)){
  serObj <- samps[[i]]
  print(serObj$orig.ident)
  print(serObj$stim)
  print(dim(serObj@assays$RNA))
  serObj <- serObj[rownames(serObj[["RNA"]]) %ni% IgGenes,]

  allCellNames <- colnames(serObj@assays$RNA)
  sampleIdentifier <- sapply(allCellNames, function(x){strsplit(x,"_")[[1]][2]})
  serObj$sampleName <- sampleIdentifier
  samps[[i]] <- serObj
}

features <- SelectIntegrationFeatures(object.list = samps)

print(str(features))
print(features)

samps <- lapply(X = samps, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE, npcs = 40)
  x <- RunUMAP(x, dims = 1:40)
  x <- RunTSNE(x, dims = 1:40)
})

print(samps)

saveRDS(samps, paste0(dataDir, "/RDSFiles/samps_08112019.rds"))

print("ALoooooo3333333")
anchors <- FindIntegrationAnchors(object.list = samps,
                                  reduction = "cca",
                                  dims = 1:40,
                                  normalization.method="LogNormalize",
                                  reference=c(1,5,7,8,9,11),
                                  scale = F)

print("1111111111")


saveRDS(anchors, paste0(dataDir, "/RDSFiles/seuratIntegratedAnchors.rds"))
samps.integrated <- IntegrateData(anchorset = anchors, dims = 1:40)
print("2222222222")

saveRDS(samps.integrated, paste0(dataDir, "/RDSFiles/seuratIntegratedSamples.rds"))

samps.integrated <- readRDS(paste0(dataDir, "/RDSFiles/seuratIntegratedSamples.rds"))
samps.integrated <- ScaleData(samps.integrated, verbose = FALSE)
samps.integrated <- RunPCA(object = samps.integrated, verbose = FALSE)
samps.integrated <- RunUMAP(samps.integrated, dims = 1:40)
samps.integrated <- RunTSNE(samps.integrated, dims = 1:40)

samps.integrated <- FindNeighbors(samps.integrated, reduction = "pca", dims = 1:40)
samps.integrated <- FindClusters(samps.integrated, resolution = 0.1)
obj.markers <- FindAllMarkers(samps.integrated,
                              min.pct = 0.25,
                              logfc.threshold = 0.25)
saveRDS(obj.markers, paste0(dataDir, "/RDSFiles/seuratIntSampObjMarkers_08112019.rds"))
saveRDS(samps.integrated, paste0(dataDir, "/RDSFiles/seuratIntegratedSamples3_08112019.rds"))


